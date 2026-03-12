#!/usr/bin/env python3
"""
=============================================================================
SNAP-10A/2 Space Reactor Benchmark - Results Analysis
=============================================================================

Reads the OpenMC statepoint file from the eigenvalue calculation and
compares the computed k-effective to reference values from ORNL/TM-2005/54.

Reference values for Case 8490a:
  k_exp  = 0.99984   (estimated experimental, -0.02$ subcritical)
  k_calc = 1.00491   (MCNP5, ENDF/B-VI, 2M+ histories)
  sigma  = 0.00069   (MCNP5 computational uncertainty)
  k_calc/k_exp = 1.0051  (normalized result)
  benchmark uncertainty < 0.008

Case description:
  28 SCA-4 fuel elements (positions 1-28 of 37)
  No Lucite rods, no beryllium inserts
  Core vessel water-flooded
  Full water reflection (9 inches in upper tank)

Outputs results.json with key comparison data.
=============================================================================
"""

import json
import glob
import openmc


# =============================================================================
# Reference values from ORNL/TM-2005/54, Table 3
# =============================================================================
# These are the reference values for Case 8490a from the ORNL benchmark report.
REFERENCE_KEFF_EXP = 0.99984       # Estimated experimental k-eff (-0.02$ subcritical)
REFERENCE_KEFF_MCNP5 = 1.00491     # MCNP5 calculated k-eff (ENDF/B-VI)
REFERENCE_SIGMA_MCNP5 = 0.00069    # MCNP5 computational uncertainty (1-sigma)
BENCHMARK_UNCERTAINTY = 0.008       # Conservative benchmark model uncertainty


# =============================================================================
# Load the most recent statepoint file
# =============================================================================
# Find all statepoint files and pick the one with the highest batch number
statepoint_files = sorted(glob.glob("statepoint.*.h5"))
if not statepoint_files:
    print("ERROR: No statepoint file found. Run model.py and then openmc first.")
    raise SystemExit(1)

# Use the last (highest batch number) statepoint
sp_file = statepoint_files[-1]
print(f"Loading statepoint: {sp_file}")
sp = openmc.StatePoint(sp_file)


# =============================================================================
# Extract k-effective results
# =============================================================================
# OpenMC provides several k-eff estimators. The combined estimate uses
# all three (collision, track-length, absorption) for best statistics.
keff_combined = sp.keff           # Combined k-eff (ufloat with uncertainty)
keff_value = keff_combined.nominal_value
keff_sigma = keff_combined.std_dev

# Number of batches and particles
n_batches = sp.n_batches
n_inactive = sp.n_inactive
n_active = n_batches - n_inactive
n_particles = sp.n_particles


# =============================================================================
# Compute comparison metrics
# =============================================================================
# Difference from experimental reference
diff_from_exp = keff_value - REFERENCE_KEFF_EXP
diff_from_exp_pcm = diff_from_exp * 1e5  # Convert to pcm (parts per 100,000)

# Difference from MCNP5 reference
diff_from_mcnp5 = keff_value - REFERENCE_KEFF_MCNP5
diff_from_mcnp5_pcm = diff_from_mcnp5 * 1e5

# Combined uncertainty for comparison with MCNP5
# (root-sum-square of both computational uncertainties)
combined_sigma = (keff_sigma**2 + REFERENCE_SIGMA_MCNP5**2)**0.5

# Number of sigma difference from MCNP5
n_sigma_from_mcnp5 = abs(diff_from_mcnp5) / combined_sigma if combined_sigma > 0 else float('inf')

# Normalized result (OpenMC k-calc / experimental k-exp)
normalized_result = keff_value / REFERENCE_KEFF_EXP

# Is the result within benchmark uncertainty?
within_benchmark_unc = abs(diff_from_exp) < BENCHMARK_UNCERTAINTY


# =============================================================================
# Print results summary
# =============================================================================
print()
print("=" * 70)
print("SNAP-10A/2 SCA-4B Benchmark Results - Case 8490a")
print("=" * 70)
print()
print("Simulation Parameters:")
print(f"  Particles per batch:  {n_particles:,}")
print(f"  Total batches:        {n_batches}")
print(f"  Inactive batches:     {n_inactive}")
print(f"  Active batches:       {n_active}")
print()
print("OpenMC Results:")
print(f"  k-eff (combined):     {keff_value:.5f} +/- {keff_sigma:.5f}")
print()
print("Reference Values (ORNL/TM-2005/54):")
print(f"  k_exp (experimental): {REFERENCE_KEFF_EXP:.5f}")
print(f"  k_calc (MCNP5):       {REFERENCE_KEFF_MCNP5:.5f} +/- {REFERENCE_SIGMA_MCNP5:.5f}")
print(f"  Benchmark unc.:       {BENCHMARK_UNCERTAINTY:.3f}")
print()
print("Comparison:")
print(f"  OpenMC - k_exp:       {diff_from_exp:+.5f}  ({diff_from_exp_pcm:+.1f} pcm)")
print(f"  OpenMC - MCNP5:       {diff_from_mcnp5:+.5f}  ({diff_from_mcnp5_pcm:+.1f} pcm)")
print(f"  Combined sigma:       {combined_sigma:.5f}")
print(f"  |OpenMC - MCNP5| / sigma_combined: {n_sigma_from_mcnp5:.2f}")
print(f"  Normalized (k_calc/k_exp):          {normalized_result:.4f}")
print(f"  Within benchmark uncertainty:       {'YES' if within_benchmark_unc else 'NO'}")
print()

# Assessment
if n_sigma_from_mcnp5 < 2.0:
    agreement = "EXCELLENT"
    detail = "Within 2-sigma of MCNP5 reference"
elif n_sigma_from_mcnp5 < 3.0:
    agreement = "GOOD"
    detail = "Within 3-sigma of MCNP5 reference"
elif within_benchmark_unc:
    agreement = "ACCEPTABLE"
    detail = "Within benchmark model uncertainty (0.008)"
else:
    agreement = "REVIEW NEEDED"
    detail = "Exceeds benchmark model uncertainty"

print(f"Assessment: {agreement} - {detail}")
print("=" * 70)


# =============================================================================
# Write results.json
# =============================================================================
results = {
    "benchmark": {
        "name": "SNAP-10A/2 SCA-4B Case 8490a",
        "source": "ORNL/TM-2005/54",
        "description": (
            "SNAP 10A/2 space reactor benchmark, SCA-4B experimental program, "
            "Case 8490a: 28 SCA-4 fuel elements, water-flooded core, "
            "full water reflection (9 in. upper tank). Simplest Phase I "
            "configuration with no Lucite rods or beryllium inserts."
        ),
        "reference_keff_experimental": REFERENCE_KEFF_EXP,
        "reference_keff_mcnp5": REFERENCE_KEFF_MCNP5,
        "reference_sigma_mcnp5": REFERENCE_SIGMA_MCNP5,
        "benchmark_uncertainty": BENCHMARK_UNCERTAINTY,
    },
    "openmc_results": {
        "keff": float(keff_value),
        "keff_sigma": float(keff_sigma),
        "particles_per_batch": int(n_particles),
        "total_batches": int(n_batches),
        "inactive_batches": int(n_inactive),
        "active_batches": int(n_active),
    },
    "comparison": {
        "diff_from_experimental": float(diff_from_exp),
        "diff_from_experimental_pcm": float(diff_from_exp_pcm),
        "diff_from_mcnp5": float(diff_from_mcnp5),
        "diff_from_mcnp5_pcm": float(diff_from_mcnp5_pcm),
        "combined_sigma": float(combined_sigma),
        "n_sigma_from_mcnp5": float(n_sigma_from_mcnp5),
        "normalized_result": float(normalized_result),
        "within_benchmark_uncertainty": bool(within_benchmark_unc),
        "assessment": agreement,
    },
}

with open("results.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults written to results.json")
