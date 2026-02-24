import os
import glob
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import neutronics_material_maker as nmm
import numpy as np
import openmc
import pandas as pd

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
openmc.config["cross_sections"] = Path.home() / "nuclear_data" / "cross_sections.xml"


REFERENCE_TEMPERATURE_C = 600.0

# User-editable expansion model per solid material.
# alpha(T) = alpha0 + alpha1*dT + alpha2*dT^2, with dT = T_C - REFERENCE_TEMPERATURE_C.
# Set alpha_unit_scale to 1e-6 when alpha coefficients are in ppm/K units.
# If rho_ref_g_cm3 is None, the script uses the library density as rho_ref.
MATERIAL_EXPANSION_MODEL = {
    "jacket": {"rho_ref_g_cm3": None, "alpha0": 12.0, "alpha1": 0.0, "alpha2": 0.0, "alpha_unit_scale": 1e-6},
    "backplate": {"rho_ref_g_cm3": None, "alpha0": 12.0, "alpha1": 0.0, "alpha2": 0.0, "alpha_unit_scale": 1e-6},
    "breeder": {"rho_ref_g_cm3": None, "alpha0": 10.0, "alpha1": 0.0, "alpha2": 0.0, "alpha_unit_scale": 1e-6},
    "be": {"rho_ref_g_cm3": None, "alpha0": 11.3, "alpha1": 0.0, "alpha2": 0.0, "alpha_unit_scale": 1e-6},
    "armour": {"rho_ref_g_cm3": None, "alpha0": 4.5, "alpha1": 0.0, "alpha2": 0.0, "alpha_unit_scale": 1e-6},
}


def alpha_linear_temp_model(temp_c, alpha0, alpha1=0.0, alpha2=0.0, alpha_unit_scale=1.0):
    """Quadratic CTE law alpha(T) in 1/K (name retained for backward compatibility)."""
    dT = temp_c - REFERENCE_TEMPERATURE_C
    alpha = alpha0 + alpha1 * dT + alpha2 * dT**2
    return alpha * alpha_unit_scale


def density_from_thermal_expansion(temp_c, rho_ref, alpha0, alpha1=0.0, alpha2=0.0, alpha_unit_scale=1.0):
    """Return density at `temp_c` using integrated volumetric expansion.

    integral(alpha dT) = (alpha0*dT + 0.5*alpha1*dT^2 + (1/3)*alpha2*dT^3) * alpha_unit_scale
    rho(T) = rho_ref * exp(-3 * integral(alpha dT))
    """
    dT = temp_c - REFERENCE_TEMPERATURE_C
    integ_alpha = (alpha0 * dT + 0.5 * alpha1 * dT**2 + (1.0 / 3.0) * alpha2 * dT**3) * alpha_unit_scale
    return rho_ref * np.exp(-3.0 * integ_alpha)


# -------------------------
# Pebble unit-cell from CSV
# -------------------------
def make_pebble_unit_cell(csv_path, pebble_mat, matrix_mat, shrink=1e-12):
    df = pd.read_csv(csv_path, skiprows=1, header=None, names=["x", "y", "z", "r"])

    rmax = float(df["r"].max())
    half = df[["x", "y", "z"]].abs().max() + rmax
    pitch = tuple((2.0 * half).to_numpy())

    Lx, Ly, Lz = pitch
    xmin, xmax = openmc.XPlane(-Lx / 2), openmc.XPlane(+Lx / 2)
    ymin, ymax = openmc.YPlane(-Ly / 2), openmc.YPlane(+Ly / 2)
    zmin, zmax = openmc.ZPlane(-Lz / 2), openmc.ZPlane(+Lz / 2)
    for s in (xmin, xmax, ymin, ymax, zmin, zmax):
        s.boundary_type = "transmission"

    box = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax

    univ = openmc.Universe(name=f"uc_{pebble_mat.name}")
    matrix_region = box

    for i, row in df.iterrows():
        r = float(row.r) * (1.0 - shrink)
        sph = openmc.Sphere(x0=float(row.x), y0=float(row.y), z0=float(row.z), r=r)
        univ.add_cell(openmc.Cell(name=f"peb_{i}", fill=pebble_mat, region=-sph))
        matrix_region &= +sph

    univ.add_cell(openmc.Cell(name="matrix", fill=matrix_mat, region=matrix_region))
    return univ, pitch


def box_region(x0, x1, y0, y1, z0, z1, boundary="transmission"):
    xp0 = openmc.XPlane(x0, boundary_type=boundary)
    xp1 = openmc.XPlane(x1, boundary_type=boundary)
    yp0 = openmc.YPlane(y0, boundary_type=boundary)
    yp1 = openmc.YPlane(y1, boundary_type=boundary)
    zp0 = openmc.ZPlane(z0, boundary_type=boundary)
    zp1 = openmc.ZPlane(z1, boundary_type=boundary)
    region = +xp0 & -xp1 & +yp0 & -yp1 & +zp0 & -zp1
    return region, (xp0, xp1, yp0, yp1, zp0, zp1)


def fill_box_with_lattice(univ, pitch, x0, x1, y0, y1, z0, z1, outer_univ=None, name="lat"):
    Lx, Ly, Lz = pitch
    nx = int(np.ceil((x1 - x0) / Lx))
    ny = int(np.ceil((y1 - y0) / Ly))
    nz = int(np.ceil((z1 - z0) / Lz))

    lat = openmc.RectLattice(name=name)
    lat.pitch = pitch
    lat.lower_left = (x0, y0, z0)
    lat.universes = [[[univ for _i in range(nx)] for _j in range(ny)] for _k in range(nz)]

    if outer_univ is not None:
        lat.outer = outer_univ

    region, _ = box_region(x0, x1, y0, y1, z0, z1)
    return openmc.Cell(name=name, fill=lat, region=region)


def make_first_wall_cooling_pipes_universe(
    x_fw_front_s,
    x_fw_rear_s,
    yp0_s,
    yp1_s,
    zp0_s,
    zp1_s,
    steel_mat,
    coolant_mat,
    ny=8,
    web=0.2,
    wall=0.1,
    skin=0.1,
    name="fw_pipes_univ",
):
    x_fw0 = float(x_fw_front_s.x0)
    x_fw1 = float(x_fw_rear_s.x0)

    fw_region = +x_fw_front_s & -x_fw_rear_s & +yp0_s & -yp1_s & +zp0_s & -zp1_s

    x_pipe0 = x_fw0 + skin
    x_pipe1 = x_fw1 - skin

    x_cool0 = x_pipe0 + wall
    x_cool1 = x_pipe1 - wall

    xpo0 = openmc.XPlane(x_pipe0, boundary_type="transmission")
    xpo1 = openmc.XPlane(x_pipe1, boundary_type="transmission")
    xpi0 = openmc.XPlane(x_cool0, boundary_type="transmission")
    xpi1 = openmc.XPlane(x_cool1, boundary_type="transmission")

    y0 = float(yp0_s.y0)
    y1 = float(yp1_s.y0)
    Ly = y1 - y0

    edge_web = 0.5 * web
    pipe_outer_y = (Ly - ny * web) / ny

    fw_univ = openmc.Universe(name=name)
    occupied_outer = None

    for i in range(ny):
        yo0 = y0 + edge_web + i * (pipe_outer_y + web)
        yo1 = yo0 + pipe_outer_y

        yo0_s = openmc.YPlane(yo0, boundary_type="transmission")
        yo1_s = openmc.YPlane(yo1, boundary_type="transmission")

        outer_pipe_region = +xpo0 & -xpo1 & +yo0_s & -yo1_s & +zp0_s & -zp1_s

        yi0 = yo0 + wall
        yi1 = yo1 - wall
        yi0_s = openmc.YPlane(yi0, boundary_type="transmission")
        yi1_s = openmc.YPlane(yi1, boundary_type="transmission")

        inner_cool_region = +xpi0 & -xpi1 & +yi0_s & -yi1_s & +zp0_s & -zp1_s
        wall_region = outer_pipe_region & ~inner_cool_region

        fw_univ.add_cell(openmc.Cell(name=f"fw_pipe_wall_{i}", fill=steel_mat, region=wall_region))
        fw_univ.add_cell(openmc.Cell(name=f"fw_pipe_cool_{i}", fill=coolant_mat, region=inner_cool_region))

        occupied_outer = outer_pipe_region if occupied_outer is None else (occupied_outer | outer_pipe_region)

    remainder = fw_region if occupied_outer is None else (fw_region & ~occupied_outer)
    fw_univ.add_cell(openmc.Cell(name="fw_webs_and_skins", fill=steel_mat, region=remainder))

    return fw_univ


@dataclass
class HCPBGeom:
    x0: float
    x1: float
    y0: float
    y1: float
    z0: float
    z1: float
    jacket_thickness: float
    pre_be_thickness: float
    armour_thickness: float
    first_wall_thickness: float
    t_be_bottom: float
    t_plate1: float
    t_breeder: float
    t_plate2: float
    t_be_top: float
    backplate_thickness: float
    backplate_mode: str = "inside"
    fw_pipe_ny: int = 4
    fw_pipe_gap: float = 0.2
    fw_pipe_wall: float = 0.1
    fw_pipe_skin: float = 0.1


def build_hcpb_module_gap(geom, breeder_uc, breeder_pitch, be_uc, be_pitch, mats, gap_front):
    module = openmc.Universe(name="hcpb_module")
    root = openmc.Universe(name="root")

    yp0 = openmc.YPlane(geom.y0, boundary_type="periodic")
    yp1 = openmc.YPlane(geom.y1, boundary_type="periodic")
    yp0.periodic_surface = yp1
    yp1.periodic_surface = yp0
    zp0 = openmc.ZPlane(geom.z0, boundary_type="reflective")
    zp1 = openmc.ZPlane(geom.z1, boundary_type="reflective")

    pre = max(0.0, float(geom.pre_be_thickness))
    tA = max(0.0, float(geom.armour_thickness))
    tFW = max(0.0, float(geom.first_wall_thickness))

    x_pre0 = geom.x0 - pre
    xb0 = openmc.XPlane(geom.x0 - pre, boundary_type="transmission")
    x_fw0 = x_pre0 - tFW
    x_fw_front = openmc.XPlane(x_fw0, boundary_type="transmission")
    x_arm0 = x_fw0 - tA
    x_arm_front = openmc.XPlane(x_arm0, boundary_type="transmission")

    wx0 = openmc.XPlane(x_arm0 - gap_front, boundary_type="vacuum")
    xp0 = openmc.XPlane(geom.x0, boundary_type="transmission")
    xp1 = openmc.XPlane(geom.x1, boundary_type="vacuum")

    void_region = +wx0 & -x_arm_front & +yp0 & -yp1 & +zp0 & -zp1
    armour_region = +x_arm_front & -x_fw_front & +yp0 & -yp1 & +zp0 & -zp1 if tA > 0 else None
    fw_region = +x_fw_front & -xb0 & +yp0 & -yp1 & +zp0 & -zp1 if tFW > 0 else None
    pre_be_region = +xb0 & -xp0 & +yp0 & -yp1 & +zp0 & -zp1 if pre > 0 else None
    module_region = +xp0 & -xp1 & +yp0 & -yp1 & +zp0 & -zp1

    root.add_cell(openmc.Cell(name="void_gap", fill=None, region=void_region))
    if armour_region is not None:
        root.add_cell(openmc.Cell(name="armour_W", fill=mats["armour"], region=armour_region))

    if fw_region is not None:
        fw_pipes_univ = make_first_wall_cooling_pipes_universe(
            x_fw_front_s=x_fw_front,
            x_fw_rear_s=xb0,
            yp0_s=yp0,
            yp1_s=yp1,
            zp0_s=zp0,
            zp1_s=zp1,
            steel_mat=mats["jacket"],
            coolant_mat=mats["he"],
            ny=geom.fw_pipe_ny,
            web=geom.fw_pipe_gap,
            wall=geom.fw_pipe_wall,
            skin=geom.fw_pipe_skin,
            name="first_wall_pipes",
        )
        root.add_cell(openmc.Cell(name="first_wall_pipes", fill=fw_pipes_univ, region=fw_region))

    if pre_be_region is not None:
        he_univ = openmc.Universe(name="he_outer_world")
        he_univ.add_cell(openmc.Cell(fill=mats["he_purge_mix"], region=None))
        pre_be_cell = fill_box_with_lattice(
            be_uc,
            be_pitch,
            x_pre0,
            geom.x0,
            geom.y0,
            geom.y1,
            geom.z0,
            geom.z1,
            outer_univ=he_univ,
            name="pre_be_slab",
        )
        pre_be_cell.region = pre_be_region
        root.add_cell(pre_be_cell)

    root.add_cell(openmc.Cell(name="module", fill=module, region=module_region))

    t = geom.jacket_thickness
    x0i = geom.x0
    z0i, z1i = geom.z0 + t, geom.z1 - t
    y_inner0, y_inner1 = geom.y0 + t, geom.y1 - t

    if geom.backplate_thickness > 0 and geom.backplate_mode == "inside":
        bp = geom.backplate_thickness
        bp_x1 = geom.x1 - t
        bp_x0 = bp_x1 - bp
        x1i = bp_x0
        bp_region, bp_surfs = box_region(bp_x0, bp_x1, y_inner0, y_inner1, z0i, z1i)
        bp_x0_surf, bp_x1_surf = bp_surfs[0], bp_surfs[1]
        module.add_cell(openmc.Cell(name="backplate", fill=mats["backplate"], region=bp_region))
    else:
        x1i = geom.x1 - t
        bp_x0_surf = None

    inner_region, _ = box_region(x0i, x1i, y_inner0, y_inner1, z0i, z1i)
    jacket_region = module_region & ~inner_region
    if geom.jacket_thickness > 0:
        module.add_cell(openmc.Cell(name="jacket", fill=mats["jacket"], region=jacket_region))

    y0b = y_inner0
    y1b = y0b + geom.t_be_bottom
    y2b = y1b + geom.t_plate1
    y3b = y2b + geom.t_breeder
    y4b = y3b + geom.t_plate2
    y5b = y4b + geom.t_be_top
    y_edges = [y0b, y1b, y2b, y3b, y4b, y5b]

    y_planes = [openmc.YPlane(y, boundary_type="transmission") for y in y_edges]
    x0i_s = xp0
    x1i_s = bp_x0_surf if bp_x0_surf is not None else xp1
    z0i_s = openmc.ZPlane(z0i, boundary_type="transmission")
    z1i_s = openmc.ZPlane(z1i, boundary_type="transmission")

    def layer_region(i):
        return +x0i_s & -x1i_s & +y_planes[i] & -y_planes[i + 1] & +z0i_s & -z1i_s

    he_univ = openmc.Universe(name="he_outer")
    he_univ.add_cell(openmc.Cell(fill=mats["he_purge_mix"], region=None))

    c = fill_box_with_lattice(be_uc, be_pitch, geom.x0, x1i_s.x0, y_edges[0], y_edges[1], z0i, z1i, outer_univ=he_univ, name="be_bottom")
    c.region = layer_region(0)
    module.add_cell(c)

    module.add_cell(openmc.Cell(name="plate1_he", fill=mats["he"], region=layer_region(1)))

    c = fill_box_with_lattice(
        breeder_uc, breeder_pitch, geom.x0, x1i_s.x0, y_edges[2], y_edges[3], z0i, z1i, outer_univ=he_univ, name="breeder_middle"
    )
    c.region = layer_region(2)
    module.add_cell(c)

    module.add_cell(openmc.Cell(name="plate2_he", fill=mats["he"], region=layer_region(3)))

    c = fill_box_with_lattice(be_uc, be_pitch, geom.x0, x1i_s.x0, y_edges[4], y_edges[5], z0i, z1i, outer_univ=he_univ, name="be_top")
    c.region = layer_region(4)
    module.add_cell(c)

    surfs = {"wx0": wx0, "armour_front": x_arm_front, "fw_front": x_fw_front, "xb0": xb0, "xp0": xp0, "xp1": xp1, "yp0": yp0, "yp1": yp1, "zp0": zp0, "zp1": zp1}
    if bp_x0_surf is not None:
        surfs["backplate_front"] = bp_x0_surf

    return openmc.Geometry(root), surfs


def build_materials(temp_c):
    temp_k = temp_c + 273.15

    he = nmm.Material.from_library("He", temperature=temp_k, pressure=8e6).openmc_material
    he_purge = nmm.Material.from_library("He", temperature=temp_k, pressure=0.2e6).openmc_material
    water = nmm.Material.from_library("H2O", temperature=temp_k, pressure=0.2e6).openmc_material
    he_purge_mix = nmm.Material.from_mixture(
        name="he_purge_mix", materials=[he_purge, water], fracs=[0.999, 0.001], percent_type="wo"
    ).openmc_material

    jacket = nmm.Material.from_library("eurofer").openmc_material
    backplate = nmm.Material.from_library("eurofer").openmc_material
    breeder = nmm.Material.from_library("Li4SiO4", enrichment=60, enrichment_target="Li6").openmc_material
    be = nmm.Material.from_library("Be").openmc_material
    tungsten = nmm.Material.from_library("tungsten").openmc_material

    # Explicit names so downstream VTK exports can map regions robustly.
    he.name = "coolant_he"
    he_purge.name = "coolant_he_purge"
    he_purge_mix.name = "coolant_he_purge_mix"
    water.name = "coolant_h2o"
    jacket.name = "steel_jacket"
    backplate.name = "steel_backplate"
    breeder.name = "breeder_li4sio4"
    be.name = "be_multiplier"
    tungsten.name = "armour_tungsten"

    for mat in (jacket, backplate, breeder, be, tungsten):
        mat.temperature = temp_k

    # Density variation from thermal expansion.
    # Configure coefficients (and optional reference densities) in
    # MATERIAL_EXPANSION_MODEL at the top of this file.
    density_meta = {}
    solids = {
        "jacket": jacket,
        "backplate": backplate,
        "breeder": breeder,
        "be": be,
        "armour": tungsten,
    }
    for key, mat in solids.items():
        cfg = MATERIAL_EXPANSION_MODEL[key]
        rho_ref = mat.get_mass_density() if cfg["rho_ref_g_cm3"] is None else float(cfg["rho_ref_g_cm3"])
        alpha_t = alpha_linear_temp_model(
            temp_c,
            cfg["alpha0"],
            cfg.get("alpha1", 0.0),
            cfg.get("alpha2", 0.0),
            cfg.get("alpha_unit_scale", 1.0),
        )
        rho_t = density_from_thermal_expansion(
            temp_c,
            rho_ref,
            cfg["alpha0"],
            cfg.get("alpha1", 0.0),
            cfg.get("alpha2", 0.0),
            cfg.get("alpha_unit_scale", 1.0),
        )
        mat.set_density("g/cm3", rho_t)
        density_meta[f"{key}_rho_ref_g_cm3"] = rho_ref
        density_meta[f"{key}_alpha_1_per_K"] = alpha_t
        density_meta[f"{key}_rho_g_cm3"] = mat.get_mass_density()

    materials = openmc.Materials([he, he_purge, jacket, backplate, breeder, be, tungsten, he_purge_mix])
    mats = {
        "he": he,
        "jacket": jacket,
        "backplate": backplate,
        "armour": tungsten,
        "he_purge": he_purge,
        "he_purge_mix": he_purge_mix,
        "breeder": breeder,
        "be": be,
    }
    return materials, mats, temp_k, density_meta


def run_case_and_get_tbr(temp_c, geom, gap_front, breeder_csv, be_csv, batches=10, particles=100000):
    materials, mats, temp_k, density_meta = build_materials(temp_c)

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
    src.space = openmc.stats.Box((x_left + eps, geom.y0 + eps, geom.z0 + eps), (x_right - eps, geom.y1 - eps, geom.z1 - eps), only_fissionable=False)
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
    src.energy = openmc.stats.muir(e0=14.1e6, m_rat=5, kt=10e3)
    model.settings.source = src

    cells = geometry.get_all_cells()
    breeder_cell = next(c for c in cells.values() if c.name == "breeder_middle")

    tbr = openmc.Tally(name="TBR")
    tbr.filters = [openmc.CellFilter([breeder_cell])]
    tbr.nuclides = ["Li6", "Li7"]
    tbr.scores = ["H3-production"]
    model.tallies = openmc.Tallies([tbr])

    for f in glob.glob("statepoint.*.h5"):
        os.remove(f)

    sp_path = model.run(tracks=False)

    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="TBR")
        tbr_val = float(t.mean.ravel()[0])
        tbr_err = float(t.std_dev.ravel()[0])

    return tbr_val, tbr_err, density_meta


def sweep_temperatures_and_plot(
    temps_c,
    geom,
    gap_front,
    breeder_csv="1.0mm_1.0cm3__polydisperse_5%var_data_4.csv",
    be_csv="1.2mm_radius_0.9cm3_data.csv",
    out_csv="tbr_vs_temperature.csv",
    out_png="tbr_vs_temperature.png",
    out_density_png="density_vs_temperature.png",
):
    rows = []
    for tc in temps_c:
        tbr_val, tbr_err, density_meta = run_case_and_get_tbr(tc, geom, gap_front, breeder_csv, be_csv)
        row = {"temperature_C": tc, "tbr": tbr_val, "tbr_std": tbr_err}
        row.update(density_meta)
        rows.append(row)
        print(f"T={tc:.1f} C -> TBR={tbr_val:.6e} +/- {tbr_err:.2e}")

    df = pd.DataFrame(rows).sort_values("temperature_C")
    df.to_csv(out_csv, index=False)

    plt.figure(figsize=(8, 5))
    plt.errorbar(df["temperature_C"], df["tbr"], yerr=df["tbr_std"], marker="o", capsize=3)
    plt.xlabel("Temperature (°C)")
    plt.ylabel("TBR (H3-production/source neutron)")
    plt.title("HCPB TBR vs Temperature")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    density_cols = [c for c in df.columns if c.endswith("_rho_g_cm3")]
    if density_cols:
        plt.figure(figsize=(8, 5))
        for c in density_cols:
            plt.plot(df["temperature_C"], df[c], marker="o", label=c.replace("_rho_g_cm3", ""))
        plt.xlabel("Temperature (°C)")
        plt.ylabel("Density (g/cm³)")
        plt.title("Material density vs temperature (thermal expansion model)")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_density_png, dpi=300)
        plt.close()

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_png}")
    if density_cols:
        print(f"Wrote {out_density_png}")


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

    gap_front = 0.2
    temps_c = [300, 450, 600, 750, 900]
    sweep_temperatures_and_plot(temps_c, geom, gap_front)
