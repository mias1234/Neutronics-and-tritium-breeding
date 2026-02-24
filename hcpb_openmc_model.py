"""Phase A helper: build your HCPB OpenMC model without running a sweep.

This module is intended for Cardinal OpenMC-only integration:
- constructs geometry/materials/source/tallies
- returns `openmc.Model`
- does not call `model.run()` in the factory
"""

from __future__ import annotations

import openmc

from hcpb_temperature_sweep import (
    HCPBGeom,
    build_hcpb_module_gap,
    build_materials,
    make_pebble_unit_cell,
)


def make_hcpb_openmc_model(
    temp_c: float,
    geom: HCPBGeom,
    gap_front: float,
    breeder_csv: str,
    be_csv: str,
    batches: int = 50,
    particles: int = 200000,
) -> openmc.Model:
    """Build and return an OpenMC model for the HCPB configuration."""
    materials, mats, temp_k, _density_meta = build_materials(temp_c)
    breeder_uc, breeder_pitch = make_pebble_unit_cell(breeder_csv, mats["breeder"], mats["he_purge_mix"])
    be_uc, be_pitch = make_pebble_unit_cell(be_csv, mats["be"], mats["he_purge_mix"])
    geometry, surfs = build_hcpb_module_gap(geom, breeder_uc, breeder_pitch, be_uc, be_pitch, mats, gap_front)

    model = openmc.Model(materials=materials, geometry=geometry, settings=openmc.Settings())
    model.settings.run_mode = "fixed source"
    model.settings.batches = batches
    model.settings.inactive = 0
    model.settings.particles = particles
    model.settings.temperature = {"method": "interpolation", "default": temp_k}

    eps = 1e-3
    x_left = surfs["wx0"].x0
    x_right = min(x for x in [surfs[k].x0 for k in ("armour_front", "fw_front", "xb0", "xp0") if k in surfs] if x > x_left)

    src = openmc.IndependentSource()
    src.space = openmc.stats.Box(
        (x_left + eps, geom.y0 + eps, geom.z0 + eps),
        (x_right - eps, geom.y1 - eps, geom.z1 - eps),
        only_fissionable=False,
    )
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    src.energy = openmc.stats.muir(e0=14.1e6, m_rat=5, kt=10e3)
    model.settings.source = src

    cells = geometry.get_all_cells()
    breeder_cell = next(c for c in cells.values() if c.name == "breeder_middle")

    tbr = openmc.Tally(name="TBR")
    tbr.filters = [openmc.CellFilter([breeder_cell])]
    tbr.nuclides = ["Li6", "Li7"]
    tbr.scores = ["H3-production"]

    all_cells = list(cells.values())
    heating_cell = openmc.Tally(name="heating_by_cell")
    heating_cell.filters = [openmc.CellFilter(all_cells)]
    heating_cell.scores = ["heating"]

    heating_total = openmc.Tally(name="heating_total")
    heating_total.scores = ["heating"]

    model.tallies = openmc.Tallies([tbr, heating_cell, heating_total])
    return model


if __name__ == "__main__":
    geom = HCPBGeom(
        x0=-40,
        x1=50,
        y0=-4.75,
        y1=4.75,
        z0=-10,
        z1=10,
        jacket_thickness=0,
        pre_be_thickness=14,
        armour_thickness=0.2,
        first_wall_thickness=2.5,
        t_be_bottom=4,
        t_plate1=0.5,
        t_breeder=0.5,
        t_plate2=0.5,
        t_be_top=4,
        backplate_thickness=10,
        backplate_mode="inside",
        fw_pipe_ny=4,
        fw_pipe_gap=0.1,
        fw_pipe_wall=0.2,
        fw_pipe_skin=0.2,
    )

    model = make_hcpb_openmc_model(
        temp_c=25,
        geom=geom,
        gap_front=0.2,
        breeder_csv="1.0mm_1.0cm3__polydisperse_5%var_data_4.csv",
        be_csv="1.2mm_radius_0.9cm3_data.csv",
        batches=10,
        particles=200000,
    )
    model.export_to_xml()
