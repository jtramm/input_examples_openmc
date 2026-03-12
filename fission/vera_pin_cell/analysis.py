#!/usr/bin/env python3
"""
Post-processing script for VERA Problem 1A results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, compares against the reference solution, and saves
results to results.json.

Reference k-eff: 1.187038 +/- 0.000054 (CE KENO-VI, ENDF/B-VII.0)
Source: CASL-U-2012-0131-004

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json
import os


def find_statepoint():
    """Find the most recent statepoint file in the current directory."""
    files = glob.glob('statepoint.*.h5')
    if not files:
        raise FileNotFoundError(
            "No statepoint files found. Run model.py --run first."
        )
    # Sort by batch number (extract integer from filename)
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
        Results dictionary with k-eff, uncertainty, comparison, etc.
    """
    import openmc

    if statepoint_path is None:
        statepoint_path = find_statepoint()

    print(f"Loading statepoint: {statepoint_path}")
    sp = openmc.StatePoint(statepoint_path)

    # --- Extract k-effective ---
    keff = sp.keff
    k_mean = keff.nominal_value       # best estimate of k-eff
    k_std = keff.std_dev              # 1-sigma statistical uncertainty

    # --- Reference value ---
    k_ref = 1.187038        # CE KENO-VI reference
    k_ref_unc = 0.000054    # 1-sigma uncertainty on reference

    # --- Compute difference ---
    # 1 pcm = 1e-5 in k-eff (percent mille)
    delta_k = (k_mean - k_ref) * 1e5  # difference in pcm
    delta_k_unc = (k_std**2 + k_ref_unc**2)**0.5 * 1e5  # combined uncertainty

    # --- Runtime information ---
    # Extract from the statepoint metadata
    runtime = sp.runtime
    total_time = runtime.get('total', 0.0) if isinstance(runtime, dict) else 0.0
    # If runtime is not a dict, try accessing as attributes
    if total_time == 0.0:
        try:
            total_time = sp.runtime['total']
        except (TypeError, KeyError):
            total_time = 0.0

    n_particles = sp.n_particles
    n_batches = sp.n_batches
    n_inactive = sp.n_inactive
    total_histories = n_particles * n_batches

    # Particle tracking rate
    if total_time > 0:
        tracking_rate = total_histories / total_time  # particles/second
    else:
        tracking_rate = 0.0

    # --- Print results ---
    print("\n" + "=" * 60)
    print("VERA Problem 1A - Results Summary")
    print("=" * 60)
    print(f"  k-eff (OpenMC):    {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  k-eff (reference): {k_ref:.6f} +/- {k_ref_unc:.6f}")
    print(f"  Difference:        {delta_k:+.1f} +/- {delta_k_unc:.1f} pcm")
    print(f"  Particles/batch:   {n_particles}")
    print(f"  Active batches:    {n_batches - n_inactive}")
    print(f"  Inactive batches:  {n_inactive}")
    print(f"  Total histories:   {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:     {total_time:.1f} s")
        print(f"  Tracking rate:     {tracking_rate:.0f} particles/s")
    print("=" * 60)

    # Assess agreement
    if abs(delta_k) < 100:
        verdict = "EXCELLENT - within 100 pcm of reference"
    elif abs(delta_k) < 500:
        verdict = "GOOD - within 500 pcm of reference"
    elif abs(delta_k) < 1000:
        verdict = "FAIR - within 1000 pcm of reference"
    else:
        verdict = "POOR - more than 1000 pcm from reference"
    print(f"  Assessment: {verdict}")
    print()

    # Note about cross section library
    print("  NOTE: Differences from the reference may be due to:")
    print("    - Different cross section library (reference used ENDF/B-VII.0)")
    print("    - Temperature treatment (nearest available vs. exact 565K)")
    print("    - Statistical uncertainty")
    print()

    # --- Save to JSON ---
    results = {
        'benchmark': 'VERA Problem 1A',
        'description': 'Single PWR pin cell, hot zero power',
        'reference': {
            'k_eff': k_ref,
            'uncertainty': k_ref_unc,
            'source': 'CE KENO-VI, ENDF/B-VII.0',
            'document': 'CASL-U-2012-0131-004'
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
            'delta_k_pcm': float(f"{delta_k:.1f}"),
            'delta_k_unc_pcm': float(f"{delta_k_unc:.1f}"),
            'assessment': verdict
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
        description='Analyze VERA Problem 1A OpenMC results'
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint file (auto-detect if omitted)')
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
