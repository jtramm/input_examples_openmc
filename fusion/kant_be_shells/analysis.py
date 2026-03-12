#!/usr/bin/env python3
"""
KANT Beryllium Spherical Shell Benchmark - Results Analysis
============================================================

Reads OpenMC statepoint files for shell configurations 1, 2, and 3.
Produces:
  - Total leakage (neutron multiplication factor) for each shell
  - Leakage spectrum plots
  - Runtime statistics
  - results.json summary

Key physics expectation: Be-9(n,2n) causes neutron multiplication, so
the total leakage should EXCEED 1.0 neutrons per source neutron.
Thicker shells produce more multiplication (up to a point).
"""

import json
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openmc


def analyze_statepoint(sp_file, shell_num):
    """Extract leakage data from a statepoint file.

    Parameters
    ----------
    sp_file : str
        Path to the statepoint file.
    shell_num : int
        Shell configuration number (1, 2, or 3).

    Returns
    -------
    dict with keys: shell, thickness, total_leakage, total_leakage_unc,
                    energy_centers, spectrum, spectrum_unc, runtime
    """
    thickness_map = {1: 5.0, 2: 10.0, 3: 17.0}
    thickness = thickness_map[shell_num]

    sp = openmc.StatePoint(sp_file)

    # Total leakage tally
    total_tally = sp.get_tally(name="total_leakage")
    total_mean = total_tally.mean.flatten()[0]
    total_std = total_tally.std_dev.flatten()[0]

    # Spectral leakage tally
    spectrum_tally = sp.get_tally(name="leakage_spectrum")
    spec_mean = spectrum_tally.mean.flatten()
    spec_std = spectrum_tally.std_dev.flatten()

    # Energy bin edges from the energy filter
    energy_filter = spectrum_tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # Shape (N, 2) with [low, high] pairs
    energy_low = energy_bins[:, 0]
    energy_high = energy_bins[:, 1]
    energy_centers = np.sqrt(energy_low * energy_high)  # Geometric mean (eV)

    # Runtime from the statepoint
    runtime = sp.runtime.get("total", 0.0) if hasattr(sp, "runtime") else 0.0

    return {
        "shell": shell_num,
        "thickness": thickness,
        "total_leakage": float(total_mean),
        "total_leakage_unc": float(total_std),
        "energy_centers_eV": energy_centers.tolist(),
        "spectrum": spec_mean.tolist(),
        "spectrum_unc": spec_std.tolist(),
        "runtime_seconds": runtime,
    }


def main():
    results = {}

    # Look for statepoint files for each shell
    for shell_num in [1, 2, 3]:
        shell_dir = f"shell_{shell_num}"
        sp_pattern = os.path.join(shell_dir, "statepoint.*.h5")
        sp_files = sorted(glob.glob(sp_pattern))

        if not sp_files:
            print(f"Shell {shell_num}: no statepoint found in {shell_dir}/")
            continue

        sp_file = sp_files[-1]  # Use the latest statepoint
        print(f"Shell {shell_num}: analyzing {sp_file}")
        results[shell_num] = analyze_statepoint(sp_file, shell_num)

    if not results:
        print("No results found. Run model.py for each shell first.")
        return

    # -----------------------------------------------------------------------
    # Print summary table
    # -----------------------------------------------------------------------
    print("\n" + "=" * 65)
    print("KANT Beryllium Spherical Shell - Leakage Summary")
    print("=" * 65)
    print(f"{'Shell':>6} {'Thickness':>10} {'Leakage':>12} {'Uncertainty':>12} {'Multiplication?':>16}")
    print("-" * 65)

    for shell_num in sorted(results.keys()):
        r = results[shell_num]
        mult = "YES" if r["total_leakage"] > 1.0 else "no"
        print(
            f"{r['shell']:>6d} {r['thickness']:>9.0f} cm "
            f"{r['total_leakage']:>11.4f} {r['total_leakage_unc']:>11.4f} "
            f"{mult:>16}"
        )

    print("=" * 65)
    print("Note: Leakage > 1.0 indicates neutron multiplication via Be-9(n,2n)")

    # -----------------------------------------------------------------------
    # Leakage spectrum plot
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {1: "blue", 2: "green", 3: "red"}
    for shell_num in sorted(results.keys()):
        r = results[shell_num]
        energy_MeV = np.array(r["energy_centers_eV"]) / 1.0e6
        spectrum = np.array(r["spectrum"])

        ax.loglog(
            energy_MeV, spectrum,
            label=f"Shell {shell_num} ({r['thickness']:.0f} cm)",
            color=colors[shell_num], linewidth=1.2
        )

    ax.set_xlabel("Neutron Energy [MeV]", fontsize=12)
    ax.set_ylabel("Surface Current [per source neutron per bin]", fontsize=12)
    ax.set_title("KANT Be Shell - Neutron Leakage Spectra", fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0.05, 15.0)

    plt.tight_layout()
    plt.savefig("leakage_spectra.png", dpi=150, bbox_inches="tight")
    print("\nSaved leakage_spectra.png")

    # -----------------------------------------------------------------------
    # Write results.json
    # -----------------------------------------------------------------------
    output = {
        "benchmark": "KANT Beryllium Spherical Shell",
        "source": "SINBAD / Forschungszentrum Karlsruhe",
        "url": "https://www.oecd-nea.org/science/wprs/shielding/sinbad/kant/fzk-be_a.htm",
        "description": (
            "Neutron leakage spectra through concentric spherical beryllium "
            "shells with a central 14.1 MeV T(d,n) source. Key physics: "
            "Be-9(n,2n) neutron multiplication."
        ),
        "shells": {},
    }

    for shell_num in sorted(results.keys()):
        r = results[shell_num]
        output["shells"][f"shell_{shell_num}"] = {
            "thickness_cm": r["thickness"],
            "inner_radius_cm": 5.0,
            "outer_radius_cm": 5.0 + r["thickness"],
            "total_leakage": r["total_leakage"],
            "total_leakage_uncertainty": r["total_leakage_unc"],
            "neutron_multiplication": r["total_leakage"] > 1.0,
            "runtime_seconds": r["runtime_seconds"],
        }

    with open("results.json", "w") as f:
        json.dump(output, f, indent=2)
    print("Saved results.json")


if __name__ == "__main__":
    main()
