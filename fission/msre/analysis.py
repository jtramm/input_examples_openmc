#!/usr/bin/env python3
"""
===============================================================================
MSRE Benchmark Analysis Script
===============================================================================

Reads OpenMC simulation results and compares the calculated k-eff to
reference values from the IRPhEP benchmark and other codes.

Prerequisites:
  Run the OpenMC simulation first:
    python model.py --model homogeneous  (or --model channel)
    openmc --threads <N>

Usage:
  python analysis.py

Outputs:
  - Console summary of k-eff comparison
  - results.json file with structured results
"""

import json
import os

import openmc


# =============================================================================
# Reference k-eff values
# =============================================================================

# These reference values come from the IRPhEP benchmark evaluation and
# the Frontiers in Nuclear Engineering (2024) comparative study.
REFERENCES = {
    "experimental": {
        "keff": 0.99978,
        "uncertainty": 0.00420,
        "source": "IRPhEP MSRE-MSR-EXP-001 (experimental benchmark)",
    },
    "openmc_csg_detailed": {
        "keff": 1.00878,
        "uncertainty": 0.00032,
        "source": "OpenMC CSG detailed model (Frontiers 2024 paper)",
    },
    "serpent": {
        "keff": 1.02132,
        "uncertainty": 0.00003,
        "source": "Serpent 2.1.30 with ENDF/B-VII.1",
    },
}


def analyze_results():
    """Analyze OpenMC statepoint results and compare to reference values.

    This function:
      1. Finds and loads the most recent statepoint file
      2. Extracts the calculated k-eff and its uncertainty
      3. Compares to experimental and computational reference values
      4. Computes the discrepancy in pcm (1 pcm = 1e-5 dk/k)
      5. Writes a structured JSON results file
    """
    # --- Find the statepoint file ---
    # OpenMC writes statepoint files named statepoint.<batch>.h5
    # We want the final one (highest batch number).
    statepoint_files = sorted(
        [f for f in os.listdir(".") if f.startswith("statepoint.") and f.endswith(".h5")]
    )

    if not statepoint_files:
        print("ERROR: No statepoint files found in the current directory.")
        print("Run the OpenMC simulation first:")
        print("  python model.py --model homogeneous")
        print("  openmc --threads <N>")
        return

    latest_sp = statepoint_files[-1]
    print(f"Loading statepoint: {latest_sp}")
    print()

    # --- Load the statepoint and extract k-eff ---
    sp = openmc.StatePoint(latest_sp)

    # k-eff is stored as a tuple: (mean, standard_deviation)
    keff_mean = sp.keff.nominal_value
    keff_std = sp.keff.std_dev

    # Number of active batches used for statistics
    n_batches = sp.n_batches
    n_inactive = sp.n_inactive
    n_active = n_batches - n_inactive
    n_particles = sp.n_particles

    print("=" * 70)
    print("MSRE Benchmark Results")
    print("=" * 70)
    print()
    print(f"Simulation parameters:")
    print(f"  Particles per batch: {n_particles:,}")
    print(f"  Total batches:       {n_batches}")
    print(f"  Inactive batches:    {n_inactive}")
    print(f"  Active batches:      {n_active}")
    print(f"  Total active histories: {n_particles * n_active:,}")
    print()

    # --- Display calculated k-eff ---
    print(f"Calculated k-eff: {keff_mean:.5f} +/- {keff_std:.5f}")
    print()

    # --- Compare to reference values ---
    print("-" * 70)
    print(f"{'Reference':<30s} {'k-eff':>10s} {'Unc.':>10s} {'Delta [pcm]':>14s}")
    print("-" * 70)

    comparisons = {}

    for ref_name, ref_data in REFERENCES.items():
        ref_keff = ref_data["keff"]
        ref_unc = ref_data["uncertainty"]

        # Difference in pcm (parts per cent mille = 1e-5)
        # Positive means our calculation is higher than the reference.
        delta_k = keff_mean - ref_keff
        delta_pcm = delta_k * 1e5  # Convert to pcm

        # Combined uncertainty (root sum of squares)
        combined_unc = (keff_std**2 + ref_unc**2) ** 0.5
        delta_sigma = abs(delta_k) / combined_unc if combined_unc > 0 else float("inf")

        label = ref_name.replace("_", " ").title()
        print(
            f"  {label:<28s} {ref_keff:>10.5f} {ref_unc:>10.5f} {delta_pcm:>+14.0f}"
        )

        comparisons[ref_name] = {
            "reference_keff": ref_keff,
            "reference_uncertainty": ref_unc,
            "delta_keff": round(delta_k, 5),
            "delta_pcm": round(delta_pcm, 1),
            "combined_uncertainty": round(combined_unc, 5),
            "delta_sigma": round(delta_sigma, 2),
            "source": ref_data["source"],
        }

    print("-" * 70)
    print()
    print("Note: Delta [pcm] = (calculated - reference) x 1e5")
    print("      Positive values mean the calculation is higher than the reference.")
    print()

    # --- Qualitative assessment ---
    exp_delta_pcm = comparisons["experimental"]["delta_pcm"]
    exp_delta_sigma = comparisons["experimental"]["delta_sigma"]

    if abs(exp_delta_pcm) < 500:
        print("Assessment: Result is within 500 pcm of the experimental value.")
    elif abs(exp_delta_pcm) < 2000:
        print("Assessment: Result is within 2000 pcm of the experimental value.")
        print("  This is expected for a simplified model (homogenized or unit cell).")
    else:
        print("Assessment: Result differs from the experimental value by more than 2000 pcm.")
        print("  This is likely due to the simplified geometry. A detailed CSG model")
        print("  with explicit graphite stringers would be needed for closer agreement.")

    if exp_delta_sigma < 2.0:
        print(f"  The difference ({exp_delta_sigma:.1f} sigma) is statistically consistent")
        print("  with the experimental benchmark within combined uncertainties.")
    else:
        print(f"  The difference ({exp_delta_sigma:.1f} sigma) indicates a systematic bias")
        print("  from the simplified model geometry.")

    # --- Write results to JSON ---
    results = {
        "benchmark": "MSRE-MSR-EXP-001",
        "description": "Molten Salt Reactor Experiment (ORNL, 1965-1969)",
        "code": "OpenMC",
        "calculated_keff": round(keff_mean, 5),
        "calculated_uncertainty": round(keff_std, 5),
        "simulation": {
            "particles_per_batch": n_particles,
            "total_batches": n_batches,
            "inactive_batches": n_inactive,
            "active_batches": n_active,
            "total_active_histories": n_particles * n_active,
        },
        "comparisons": comparisons,
    }

    results_file = "results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)

    print()
    print(f"Results written to: {results_file}")
    print()


if __name__ == "__main__":
    analyze_results()
