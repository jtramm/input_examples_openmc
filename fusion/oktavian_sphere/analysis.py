#!/usr/bin/env python3
"""
OKTAVIAN Iron Sphere -- Post-processing and Analysis
=====================================================

Loads the statepoint file produced by OpenMC and extracts:

  1. Total neutron leakage (integral of surface current tally)
  2. Leakage energy spectrum (per-bin surface current values)
  3. Cell flux spectrum in the iron sphere
  4. Runtime and particle tracking rate

Results are printed to the console and saved to results.json.

Prerequisites
-------------
Run model.py --run (or openmc after exporting XML) to produce a statepoint file.
"""

import glob
import json
import sys

import numpy as np

import openmc


def find_statepoint() -> str:
    """Find the most recent statepoint file in the current directory."""
    candidates = sorted(glob.glob("statepoint.*.h5"))
    if not candidates:
        print("ERROR: No statepoint file found. Run the simulation first.")
        sys.exit(1)
    # Return the last one (highest batch number)
    return candidates[-1]


def main():
    sp_path = find_statepoint()
    print(f"Loading statepoint: {sp_path}")
    sp = openmc.StatePoint(sp_path)

    results = {}

    # =====================================================================
    # 1. Neutron leakage spectrum (surface current tally)
    # =====================================================================
    leakage_tally = sp.get_tally(name="Neutron leakage spectrum")

    # Extract the mean and standard deviation arrays
    # Shape: (energy_bins, 1) after squeezing filters
    leakage_mean = leakage_tally.mean.flatten()
    leakage_std = leakage_tally.std_dev.flatten()

    # Get the energy bin edges from the energy filter
    energy_filter = leakage_tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): lower and upper edges

    # Bin centers (geometric mean for log-spaced bins)
    e_low = energy_bins[:, 0]
    e_high = energy_bins[:, 1]
    e_center = np.sqrt(e_low * e_high)

    # Total leakage = sum over all energy bins
    total_leakage = np.sum(leakage_mean)
    total_leakage_std = np.sqrt(np.sum(leakage_std**2))

    print("\n" + "=" * 60)
    print("NEUTRON LEAKAGE (Surface Current on Iron Sphere)")
    print("=" * 60)
    print(f"  Total leakage    : {total_leakage:.6f} +/- {total_leakage_std:.6f}")
    print(f"  (per source particle)")
    print()

    # Print a summary table of the leakage spectrum (selected bins)
    print(f"  {'Energy [MeV]':>14s}  {'Leakage':>12s}  {'Std Dev':>12s}")
    print(f"  {'-'*14}  {'-'*12}  {'-'*12}")
    # Show every 10th bin for brevity
    for i in range(0, len(leakage_mean), 10):
        e_mev = e_center[i] / 1.0e6  # convert eV to MeV
        print(f"  {e_mev:14.4f}  {leakage_mean[i]:12.6e}  {leakage_std[i]:12.6e}")

    # Store leakage spectrum in results
    results["leakage"] = {
        "total": float(total_leakage),
        "total_std": float(total_leakage_std),
        "energy_centers_eV": e_center.tolist(),
        "energy_low_eV": e_low.tolist(),
        "energy_high_eV": e_high.tolist(),
        "mean": leakage_mean.tolist(),
        "std_dev": leakage_std.tolist(),
    }

    # =====================================================================
    # 2. Cell flux spectrum in the iron sphere
    # =====================================================================
    flux_tally = sp.get_tally(name="Iron sphere flux spectrum")

    flux_mean = flux_tally.mean.flatten()
    flux_std = flux_tally.std_dev.flatten()

    total_flux = np.sum(flux_mean)
    total_flux_std = np.sqrt(np.sum(flux_std**2))

    print("\n" + "=" * 60)
    print("IRON SPHERE FLUX (Cell Flux Tally)")
    print("=" * 60)
    print(f"  Integrated flux  : {total_flux:.6e} +/- {total_flux_std:.6e}")
    print(f"  (track-length per source particle)")

    results["flux"] = {
        "total": float(total_flux),
        "total_std": float(total_flux_std),
        "mean": flux_mean.tolist(),
        "std_dev": flux_std.tolist(),
    }

    # =====================================================================
    # 3. Runtime and tracking rate
    # =====================================================================
    print("\n" + "=" * 60)
    print("SIMULATION PERFORMANCE")
    print("=" * 60)

    runtime = sp.runtime
    # runtime is a dict with keys like 'total', 'reading cross sections', etc.
    total_time = runtime.get("total", 0.0)
    transport_time = runtime.get("transport", runtime.get("inactive batches", 0.0))

    n_particles = sp.n_particles * sp.n_batches
    if transport_time > 0:
        tracking_rate = n_particles / transport_time
    else:
        tracking_rate = 0.0

    print(f"  Total runtime       : {total_time:.2f} s")
    print(f"  Transport time      : {transport_time:.2f} s")
    print(f"  Particles simulated : {n_particles:,}")
    if tracking_rate > 0:
        print(f"  Tracking rate       : {tracking_rate:,.0f} particles/s")

    results["performance"] = {
        "total_runtime_s": float(total_time),
        "transport_time_s": float(transport_time),
        "particles": int(n_particles),
        "tracking_rate_per_s": float(tracking_rate),
    }

    # =====================================================================
    # Save results to JSON
    # =====================================================================
    output_file = "results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")


if __name__ == "__main__":
    main()
