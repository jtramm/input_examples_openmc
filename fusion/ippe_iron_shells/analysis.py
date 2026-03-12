#!/usr/bin/env python3
"""
================================================================================
IPPE Iron Spherical Shell Transmission Benchmark - Post-Processing Analysis
================================================================================

Reads the OpenMC statepoint file produced by model.py and extracts:
  1. Total neutron leakage fraction through the iron shell
  2. Neutron leakage energy spectrum
  3. Runtime statistics
  4. Saves all results to results.json

Usage:
    python analysis.py                          # Default: shell 1
    python analysis.py --shell 3                # Analyze shell 3
    python analysis.py --statepoint statepoint.1.h5  # Specify file
================================================================================
"""

import argparse
import glob
import json
import os
import time

import numpy as np
import openmc

# Import shell data from model.py for reference
from model import SHELL_DATA


def find_statepoint():
    """
    Find the most recent statepoint file in the current directory.

    Returns
    -------
    str
        Path to the statepoint file.
    """
    # Look for statepoint files matching the standard naming pattern
    candidates = sorted(glob.glob("statepoint.*.h5"))
    if not candidates:
        raise FileNotFoundError(
            "No statepoint.*.h5 files found in the current directory. "
            "Run model.py first to generate simulation output."
        )
    # Return the most recent one (highest batch number)
    return candidates[-1]


def analyze(shell_number: int, statepoint_path: str = None):
    """
    Analyze the simulation results for the specified shell configuration.

    Parameters
    ----------
    shell_number : int
        Shell configuration number (1-5), used for labeling.
    statepoint_path : str, optional
        Path to statepoint file. If None, auto-detect.

    Returns
    -------
    dict
        Dictionary of results (also saved to results.json).
    """

    # =========================================================================
    # Load statepoint
    # =========================================================================
    if statepoint_path is None:
        statepoint_path = find_statepoint()

    print(f"Loading statepoint: {statepoint_path}")
    sp = openmc.StatePoint(statepoint_path)

    # Retrieve shell geometry info for labeling
    shell = SHELL_DATA[shell_number]
    r_inner = shell["inner_radius"]
    r_outer = r_inner + shell["thickness"]

    print(f"\nIPPE Iron Shell #{shell_number}: {shell['description']}")
    print(f"  Inner radius:   {r_inner:.1f} cm")
    print(f"  Wall thickness: {shell['thickness']:.1f} cm")
    print(f"  Outer radius:   {r_outer:.1f} cm")

    # =========================================================================
    # Extract total leakage fraction
    # =========================================================================
    # The "Total leakage" tally gives the integral neutron current through
    # the outer surface, normalized per source particle.
    total_tally = sp.get_tally(name="Total leakage")
    total_mean = float(total_tally.mean.flatten()[0])
    total_std_raw = total_tally.std_dev.flatten()[0]
    # With a single batch, std_dev will be nan; replace with 0 for display
    total_std = float(np.nan_to_num(total_std_raw, nan=0.0))

    print(f"\n{'='*60}")
    print(f"TOTAL NEUTRON LEAKAGE FRACTION")
    print(f"{'='*60}")
    print(f"  Leakage fraction:  {total_mean:.6f} +/- {total_std:.6f}")
    print(f"  Relative error:    {total_std/total_mean*100:.2f}%" if total_mean > 0 else "  N/A")
    print(f"  Attenuation:       {1.0 - total_mean:.6f}")

    # =========================================================================
    # Extract leakage spectrum
    # =========================================================================
    # The "Leakage spectrum" tally provides energy-resolved neutron current
    # through the outer surface.
    spectrum_tally = sp.get_tally(name="Leakage spectrum")

    # Get the energy filter to extract bin edges
    energy_filter = spectrum_tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # Shape: (N_bins, 2) with [low, high] edges

    # Extract mean values and uncertainties for each energy bin
    spectrum_mean = spectrum_tally.mean.flatten()     # Per-bin leakage
    # Replace nan uncertainties (single-batch runs) with 0
    spectrum_std = np.nan_to_num(
        spectrum_tally.std_dev.flatten(), nan=0.0
    )

    # Compute bin centers (geometric mean for log-spaced bins) and widths
    energy_low = energy_bins[:, 0]      # Lower bin edges (eV)
    energy_high = energy_bins[:, 1]     # Upper bin edges (eV)
    energy_centers = np.sqrt(energy_low * energy_high)  # Geometric mean (eV)
    energy_widths = energy_high - energy_low             # Bin widths (eV)

    # Compute flux per unit lethargy for spectral plotting
    # d(phi)/d(ln E) = E * d(phi)/dE = spectrum_mean / ln(E_high/E_low)
    lethargy_widths = np.log(energy_high / energy_low)
    spectrum_per_lethargy = spectrum_mean / lethargy_widths

    print(f"\n{'='*60}")
    print(f"LEAKAGE SPECTRUM SUMMARY")
    print(f"{'='*60}")
    print(f"  Energy range:  {energy_low[0]/1e6:.4f} - {energy_high[-1]/1e6:.2f} MeV")
    print(f"  Number of bins: {len(spectrum_mean)}")
    print(f"  Peak bin energy: {energy_centers[np.argmax(spectrum_mean)]/1e6:.3f} MeV")
    print(f"  Integrated spectrum: {np.sum(spectrum_mean):.6f}")

    # Print top 10 energy bins by leakage
    print(f"\n  Top 10 energy bins by leakage:")
    sorted_idx = np.argsort(spectrum_mean)[::-1]
    for i, idx in enumerate(sorted_idx[:10]):
        print(f"    {i+1:2d}. {energy_low[idx]/1e6:.4f} - {energy_high[idx]/1e6:.4f} MeV: "
              f"{spectrum_mean[idx]:.6e} +/- {spectrum_std[idx]:.6e}")

    # =========================================================================
    # Runtime statistics
    # =========================================================================
    print(f"\n{'='*60}")
    print(f"RUNTIME STATISTICS")
    print(f"{'='*60}")

    # Extract runtime info from the statepoint
    runtime = sp.runtime
    print(f"  Total runtime:       {runtime['total']:.2f} s")
    print(f"  Transport time:      {runtime['transport']:.2f} s")
    print(f"  Initialization:      {runtime['total initialization']:.2f} s")
    print(f"  Reading XS:          {runtime['reading cross sections']:.2f} s")

    n_particles = sp.n_particles
    transport_time = runtime['transport']
    if transport_time > 0:
        rate = n_particles / transport_time
        print(f"  Particles simulated: {n_particles:,}")
        print(f"  Transport rate:      {rate:,.0f} particles/sec")

    # =========================================================================
    # Build results dictionary
    # =========================================================================
    results = {
        "benchmark": "IPPE Iron Spherical Shell Transmission",
        "source": "SINBAD Database, OECD-NEA",
        "shell_number": shell_number,
        "geometry": {
            "inner_radius_cm": r_inner,
            "wall_thickness_cm": shell["thickness"],
            "outer_radius_cm": r_outer,
            "description": shell["description"],
        },
        "total_leakage": {
            "mean": total_mean,
            "std_dev": total_std,
            "relative_error_pct": (total_std / total_mean * 100) if total_mean > 0 else None,
        },
        "spectrum": {
            "energy_low_eV": energy_low.tolist(),
            "energy_high_eV": energy_high.tolist(),
            "energy_center_eV": energy_centers.tolist(),
            "leakage_per_source": spectrum_mean.tolist(),
            "uncertainty": spectrum_std.tolist(),
            "leakage_per_lethargy": spectrum_per_lethargy.tolist(),
        },
        "runtime": {
            "total_seconds": runtime["total"],
            "transport_seconds": runtime["transport"],
            "initialization_seconds": runtime["total initialization"],
            "particles": int(n_particles),
            "transport_rate_per_sec": (
                n_particles / transport_time if transport_time > 0 else None
            ),
        },
    }

    # =========================================================================
    # Save to JSON
    # =========================================================================
    output_file = "results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_file}")

    return results


def main():
    """Parse arguments and run the analysis."""

    parser = argparse.ArgumentParser(
        description="Analyze IPPE Iron Shell benchmark simulation results."
    )
    parser.add_argument(
        "--shell", type=int, default=1, choices=[1, 2, 3, 4, 5],
        help="Shell configuration number (1-5). Default: 1"
    )
    parser.add_argument(
        "--statepoint", type=str, default=None,
        help="Path to statepoint file. Default: auto-detect."
    )
    args = parser.parse_args()

    analyze(args.shell, args.statepoint)


if __name__ == "__main__":
    main()
