"""Phase B helper utilities: one-way thermal -> neutronics feedback.

This module helps you move from Phase A to Phase B by updating OpenMC
material temperatures and densities from a prescribed temperature field,
then running neutronics and extracting TBR/heating.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import openmc

from hcpb_openmc_model import make_hcpb_openmc_model
from hcpb_temperature_sweep import (
    HCPBGeom,
    MATERIAL_EXPANSION_MODEL,
    alpha_linear_temp_model,
    density_from_thermal_expansion,
)


@dataclass
class PhaseBResult:
    tbr: float
    tbr_std: float
    heating_total: float
    heating_total_std: float


def apply_material_temperature_feedback(
    model: openmc.Model,
    region_temperatures_c: Mapping[str, float],
    *,
    reference_temperature_c: float,
) -> None:
    """Apply region-wise temperature feedback to OpenMC materials.

    Parameters
    ----------
    model:
        OpenMC model whose materials will be modified in-place.
    region_temperatures_c:
        Mapping from expansion-model key to temperature in Celsius.
        Expected keys (matching MATERIAL_EXPANSION_MODEL):
        - jacket
        - backplate
        - breeder
        - be
        - armour
    reference_temperature_c:
        Reference temperature (°C) used for `rho_ref` acquisition.

    Notes
    -----
    - If `rho_ref_g_cm3` is None for a material, this function approximates
      rho_ref by reversing the expansion law from the current material state
      to `reference_temperature_c`.
    """
    mats_by_name = {m.name: m for m in model.materials}

    name_map = {
        "jacket": "eurofer",
        "backplate": "eurofer",
        "breeder": "Li4SiO4",
        "be": "Be",
        "armour": "tungsten",
    }

    for key, temp_c in region_temperatures_c.items():
        if key not in MATERIAL_EXPANSION_MODEL:
            continue

        # Resolve material object by library name heuristic.
        lib_name = name_map.get(key)
        candidates = [m for m in mats_by_name.values() if lib_name and lib_name.lower() in (m.name or "").lower()]
        if not candidates:
            continue
        mat = candidates[0]

        cfg = MATERIAL_EXPANSION_MODEL[key]
        alpha0 = cfg["alpha0"]
        alpha1 = cfg.get("alpha1", 0.0)
        alpha2 = cfg.get("alpha2", 0.0)
        alpha_unit_scale = cfg.get("alpha_unit_scale", 1.0)

        # Temperature in OpenMC is K
        mat.temperature = temp_c + 273.15

        rho_ref = cfg.get("rho_ref_g_cm3")
        if rho_ref is None:
            # Use current mass density as an anchor; estimate rho_ref at T_ref.
            rho_now = mat.get_mass_density()
            rho_at_tref_from_unit_ref = density_from_thermal_expansion(
                reference_temperature_c,
                rho_ref=1.0,
                alpha0=alpha0,
                alpha1=alpha1,
                alpha2=alpha2,
                alpha_unit_scale=alpha_unit_scale,
            )
            rho_now_from_unit_ref = density_from_thermal_expansion(
                temp_c,
                rho_ref=1.0,
                alpha0=alpha0,
                alpha1=alpha1,
                alpha2=alpha2,
                alpha_unit_scale=alpha_unit_scale,
            )
            rho_ref = rho_now * rho_at_tref_from_unit_ref / rho_now_from_unit_ref

        rho_t = density_from_thermal_expansion(
            temp_c,
            rho_ref=float(rho_ref),
            alpha0=alpha0,
            alpha1=alpha1,
            alpha2=alpha2,
            alpha_unit_scale=alpha_unit_scale,
        )
        mat.set_density("g/cm3", rho_t)

        # optional computed alpha(T), currently for side-effect parity/debug potential
        _ = alpha_linear_temp_model(
            temp_c,
            alpha0=alpha0,
            alpha1=alpha1,
            alpha2=alpha2,
            alpha_unit_scale=alpha_unit_scale,
        )


def run_phase_b_iteration(
    *,
    geom: HCPBGeom,
    gap_front: float,
    breeder_csv: str,
    be_csv: str,
    region_temperatures_c: Mapping[str, float],
    model_reference_temp_c: float = 25.0,
    expansion_reference_temp_c: float = 0.0,
    batches: int = 10,
    particles: int = 200000,
) -> PhaseBResult:
    """Build a model, apply one-way thermal feedback, run OpenMC, return KPIs."""
    model = make_hcpb_openmc_model(
        temp_c=model_reference_temp_c,
        geom=geom,
        gap_front=gap_front,
        breeder_csv=breeder_csv,
        be_csv=be_csv,
        batches=batches,
        particles=particles,
    )

    apply_material_temperature_feedback(
        model,
        region_temperatures_c=region_temperatures_c,
        reference_temperature_c=expansion_reference_temp_c,
    )

    sp_path = model.run(tracks=False)
    with openmc.StatePoint(sp_path) as sp:
        tbr = sp.get_tally(name="TBR")
        heating = sp.get_tally(name="heating_total")
        return PhaseBResult(
            tbr=float(tbr.mean.ravel()[0]),
            tbr_std=float(tbr.std_dev.ravel()[0]),
            heating_total=float(heating.mean.ravel()[0]),
            heating_total_std=float(heating.std_dev.ravel()[0]),
        )
