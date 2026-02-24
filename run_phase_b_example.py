"""Example: run one Phase B (one-way thermal->neutronics) iteration."""

from phase_b_temperature_feedback import run_phase_b_iteration
from hcpb_temperature_sweep import HCPBGeom


def main() -> None:
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

    # Prescribed temperatures from thermal side (MOOSE/NekRS) for Phase B.
    region_temperatures_c = {
        "jacket": 420.0,
        "backplate": 380.0,
        "breeder": 650.0,
        "be": 500.0,
        "armour": 700.0,
    }

    result = run_phase_b_iteration(
        geom=geom,
        gap_front=0.2,
        breeder_csv="1.0mm_1.0cm3__polydisperse_5%var_data_4.csv",
        be_csv="1.2mm_radius_0.9cm3_data.csv",
        region_temperatures_c=region_temperatures_c,
        model_reference_temp_c=25.0,
        expansion_reference_temp_c=0.0,
        batches=10,
        particles=200000,
    )

    print(f"Phase B result: TBR = {result.tbr:.6e} +/- {result.tbr_std:.2e}")
    print(f"Phase B result: heating_total = {result.heating_total:.6e} +/- {result.heating_total_std:.2e}")


if __name__ == "__main__":
    main()
