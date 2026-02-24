from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_tbr_vs_temperature(
    csv_path: str = "tbr_vs_temperature.csv",
    out_png: str = "tbr_vs_temperature_from_csv.png",
    title: str = "HCPB TBR vs Temperature",
) -> None:
    """Plot TBR vs temperature using results from hcpb_temperature_sweep.py.

    Expected CSV columns:
      - temperature_C
      - tbr
      - optional tbr_std
    """
    df = pd.read_csv(csv_path)

    required = {"temperature_C", "tbr"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required CSV columns: {sorted(missing)}")

    df = df.sort_values("temperature_C")

    plt.figure(figsize=(8, 5))
    if "tbr_std" in df.columns:
        plt.errorbar(
            df["temperature_C"],
            df["tbr"],
            yerr=df["tbr_std"],
            marker="o",
            capsize=3,
            linewidth=1.5,
        )
    else:
        plt.plot(df["temperature_C"], df["tbr"], marker="o", linewidth=1.5)

    plt.xlabel("Temperature (°C)")
    plt.ylabel("TBR (H3-production/source neutron)")
    plt.title(title)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Loaded: {csv_path}")
    print(f"Wrote:  {out_png}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot TBR vs temperature from sweep CSV output.")
    parser.add_argument("--csv", default="tbr_vs_temperature.csv", help="Input CSV from temperature sweep")
    parser.add_argument("--out", default="tbr_vs_temperature_from_csv.png", help="Output PNG filename")
    parser.add_argument("--title", default="HCPB TBR vs Temperature", help="Plot title")
    args = parser.parse_args()

    plot_tbr_vs_temperature(csv_path=args.csv, out_png=args.out, title=args.title)


if __name__ == "__main__":
    main()
