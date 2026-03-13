#!/usr/bin/env python3
"""
Post-processing script for VVER-1000 mock-up benchmark results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, and compares against reference MCS/MCNP6 calculations.

Reference k-eff values (MCS with ENDF/B-VII.1, from Table 3):
  Case 1: 1.00137    Case 2: 1.00307    Case 3: 1.00403
  Case 4: 1.00462    Case 5: 1.00495    Case 6: 1.00532

The experimental k-eff is 1.0 for all cases (critical configurations).
MCS calculations overpredicted by +137 to +532 pcm with ENDF/B-VII.1.

This single-assembly model with reflective BCs will NOT match the full
mock-up keff exactly (the mock-up includes 32 assemblies, external
components, and varying moderator levels). The comparison is indicative
of the assembly-level reactivity and cross-section behavior.

Source:
    F. Setiawan et al., Nuclear Engineering and Technology 53 (2021) 1-18.

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json
import os

import numpy as np


# Reference keff values from MCS with ENDF/B-VII.1 (Table 3 of paper)
# These are for the FULL mock-up model, not a single assembly.
# The experimental keff = 1.0 for all cases.
REFERENCE_CASES = {
    'Case 1': {'keff': 1.00137, 'boron_g_per_kg': 2.85, 'mod_level_cm': 51.34},
    'Case 2': {'keff': 1.00307, 'boron_g_per_kg': 3.63, 'mod_level_cm': 65.91},
    'Case 3': {'keff': 1.00403, 'boron_g_per_kg': 4.06, 'mod_level_cm': 79.11},
    'Case 4': {'keff': 1.00462, 'boron_g_per_kg': 4.44, 'mod_level_cm': 96.71},
    'Case 5': {'keff': 1.00495, 'boron_g_per_kg': 4.53, 'mod_level_cm': 103.37},
    'Case 6': {'keff': 1.00532, 'boron_g_per_kg': 4.68, 'mod_level_cm': 150.0},
}

# For single-assembly infinite lattice, there is no direct experimental
# reference. We note the experimental keff = 1.0 for the full mock-up
# and the MCS overprediction range.
K_EXP = 1.0
MCS_OVERPREDICTION_RANGE_PCM = (137, 559)  # VII.1 and VIII.0 combined


def find_statepoint():
    """Find the most recent statepoint file in the current directory."""
    files = glob.glob('statepoint.*.h5')
    if not files:
        raise FileNotFoundError(
            "No statepoint files found. Run model.py --run first."
        )
    files.sort(key=lambda f: int(f.split('.')[1]))
    return files[-1]


def analyze(statepoint_path=None):
    """
    Load statepoint, extract results, compare to reference.

    Parameters
    ----------
    statepoint_path : str or None
        Path to statepoint HDF5 file. If None, auto-detect.

    Returns
    -------
    dict
        Results dictionary with k-eff, uncertainty, and comparison data.
    """
    import openmc

    if statepoint_path is None:
        statepoint_path = find_statepoint()

    print(f"Loading statepoint: {statepoint_path}")
    sp = openmc.StatePoint(statepoint_path)

    # --- Extract k-effective ---
    keff = sp.keff
    k_mean = keff.nominal_value
    k_std = keff.std_dev

    # --- Compute difference from experiment (keff = 1.0) ---
    delta_k_exp = (k_mean - K_EXP) * 1e5  # pcm
    delta_k_exp_unc = k_std * 1e5

    # --- Runtime information ---
    runtime = sp.runtime
    total_time = runtime.get('total', 0.0) if isinstance(runtime, dict) else 0.0
    if total_time == 0.0:
        try:
            total_time = sp.runtime['total']
        except (TypeError, KeyError):
            total_time = 0.0

    n_particles = sp.n_particles
    n_batches = sp.n_batches
    n_inactive = sp.n_inactive
    total_histories = n_particles * n_batches

    if total_time > 0:
        tracking_rate = total_histories / total_time
    else:
        tracking_rate = 0.0

    # --- Print results ---
    print()
    print("=" * 70)
    print("VVER-1000 Mock-up Benchmark - Single Hexagonal Assembly Results")
    print("Hexagonal lattice, 312 fuel pins, 18 guide tubes, 1 central tube")
    print("=" * 70)
    print(f"  k-eff (OpenMC):         {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  k-eff (experiment):     {K_EXP:.5f} (full mock-up, critical)")
    print(f"  Difference from exp:    {delta_k_exp:+.1f} +/- "
          f"{delta_k_exp_unc:.1f} pcm")
    print(f"  Particles/batch:        {n_particles}")
    print(f"  Active batches:         {n_batches - n_inactive}")
    print(f"  Inactive batches:       {n_inactive}")
    print(f"  Total histories:        {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:          {total_time:.1f} s")
        print(f"  Tracking rate:          {tracking_rate:.0f} particles/s")
    print("=" * 70)

    # --- Context for comparison ---
    print()
    print("Reference context (full 32-assembly mock-up, MCS with ENDF/B-VII.1):")
    print(f"  keff overprediction range: +{MCS_OVERPREDICTION_RANGE_PCM[0]} to "
          f"+{MCS_OVERPREDICTION_RANGE_PCM[1]} pcm vs experiment")
    print()
    print("NOTE: This is a SINGLE ASSEMBLY infinite lattice calculation.")
    print("Direct comparison to the full mock-up keff is not meaningful because:")
    print("  - The mock-up has 32 assemblies with different enrichments")
    print("  - Leakage, reflector, and baffle effects are absent here")
    print("  - The moderator level varies across critical configurations")
    print("  - Control rod insertion (Case 6) is not modeled")
    print()
    print("This model validates the hexagonal lattice construction and provides")
    print("a starting point for more detailed VVER benchmark calculations.")
    print()

    # --- Assess result ---
    if k_mean > 0.9 and k_mean < 1.5:
        verdict = "PHYSICALLY REASONABLE"
    else:
        verdict = "CHECK MODEL - keff outside expected range"

    # --- Save to JSON ---
    results = {
        'benchmark': 'VVER-1000 Mock-up LR-0',
        'description': 'Single hexagonal fuel assembly, 2D infinite lattice',
        'geometry': {
            'lattice_type': 'hexagonal (HexLattice)',
            'n_rings': 11,
            'n_fuel_pins': 312,
            'n_guide_tubes': 18,
            'n_central_tubes': 1,
            'total_positions': 331,
            'pin_pitch_cm': 1.275,
            'assembly_pitch_cm': 23.6,
            'orientation': 'x',
        },
        'reference': {
            'experiment_keff': K_EXP,
            'mcs_overprediction_pcm': list(MCS_OVERPREDICTION_RANGE_PCM),
            'source': ('F. Setiawan et al., Nuclear Engineering and '
                       'Technology 53 (2021) 1-18'),
            'note': 'Reference is for full 32-assembly mock-up, not single assembly'
        },
        'openmc': {
            'k_eff': float(k_mean),
            'uncertainty': float(k_std),
            'particles_per_batch': int(n_particles),
            'active_batches': int(n_batches - n_inactive),
            'inactive_batches': int(n_inactive),
            'total_histories': int(total_histories),
        },
        'comparison': {
            'delta_k_vs_experiment_pcm': float(f"{delta_k_exp:.1f}"),
            'delta_k_uncertainty_pcm': float(f"{delta_k_exp_unc:.1f}"),
            'assessment': verdict,
        }
    }

    if total_time > 0:
        results['openmc']['runtime_seconds'] = float(f"{total_time:.1f}")
        results['openmc']['tracking_rate'] = float(f"{tracking_rate:.0f}")

    output_file = 'results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {output_file}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze VVER-1000 mock-up benchmark OpenMC results'
    )
    parser.add_argument(
        '--statepoint', type=str, default=None,
        help='Path to statepoint file (auto-detect if omitted)'
    )
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
