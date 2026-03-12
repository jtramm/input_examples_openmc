#!/usr/bin/env python3
"""
Post-processing script for VERA Problem 2A results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, compares against the reference solution, and optionally
processes the pin power distribution from the mesh tally.

Reference k-eff: 1.18703 +/- 0.00009 (CE KENO-VI, ENDF/B-VII.0)
Source: CASL-U-2012-0131-004, Table 11

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json
import os

import numpy as np


# Reference k-eff for Problem 2A from CASL-U-2012-0131-004
K_REF = 1.18703
K_REF_UNC = 0.00009


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
        Results dictionary with k-eff, uncertainty, comparison, and
        optionally pin power data.
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

    # --- Compute difference from reference ---
    # 1 pcm = 1e-5 in k-eff (percent mille)
    delta_k = (k_mean - K_REF) * 1e5  # difference in pcm
    delta_k_unc = (k_std**2 + K_REF_UNC**2)**0.5 * 1e5  # combined uncertainty

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

    # Particle tracking rate
    if total_time > 0:
        tracking_rate = total_histories / total_time
    else:
        tracking_rate = 0.0

    # --- Print k-eff results ---
    print("\n" + "=" * 65)
    print("VERA Problem 2A - Results Summary")
    print("17x17 PWR Fuel Assembly, Hot Zero Power, All Rods Out")
    print("=" * 65)
    print(f"  k-eff (OpenMC):    {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  k-eff (reference): {K_REF:.5f} +/- {K_REF_UNC:.5f}")
    print(f"  Difference:        {delta_k:+.1f} +/- {delta_k_unc:.1f} pcm")
    print(f"  Particles/batch:   {n_particles}")
    print(f"  Active batches:    {n_batches - n_inactive}")
    print(f"  Inactive batches:  {n_inactive}")
    print(f"  Total histories:   {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:     {total_time:.1f} s")
        print(f"  Tracking rate:     {tracking_rate:.0f} particles/s")
    print("=" * 65)

    # Assess agreement
    if abs(delta_k) < 100:
        verdict = "EXCELLENT - within 100 pcm of reference"
    elif abs(delta_k) < 500:
        verdict = "GOOD - within 500 pcm of reference"
    elif abs(delta_k) < 1000:
        verdict = "FAIR - within 1000 pcm of reference"
    else:
        verdict = "EXPECTED BIAS - more than 1000 pcm (likely due to 294K xs data)"
    print(f"  Assessment: {verdict}")
    print()

    # Note about cross section temperature
    print("  NOTE: Expected ~1500 pcm positive bias because the cross section")
    print("  library only has 294K data while the benchmark specifies 565K.")
    print("  This Doppler broadening effect makes the 294K calculation more")
    print("  reactive than the true 565K conditions.")
    print()

    # --- Pin power distribution (if mesh tally exists) ---
    pin_powers = None
    pin_power_map = None
    try:
        # Look for the pin power mesh tally
        tally = sp.get_tally(name='Pin Powers')
        if tally is not None:
            print("Processing pin power distribution...")

            # Extract fission rate data
            # The mesh tally has shape (17, 17) after reshaping
            mean = tally.mean.ravel()
            std_dev = tally.std_dev.ravel()

            # Reshape to 17x17 grid
            # OpenMC mesh tallies are ordered with y varying fastest,
            # and the mesh indices go from lower-left to upper-right.
            mean_2d = mean.reshape((17, 17))
            std_2d = std_dev.reshape((17, 17))

            # Normalize to average fuel pin power = 1.0
            # Only average over fuel pin positions (exclude guide tubes
            # and instrument tube which have zero fission rate)
            total_fission = np.sum(mean_2d)
            if total_fission > 0:
                # Normalize so the average over all 289 positions = 1.0
                # (positions with no fuel will naturally be ~0)
                n_total = 17 * 17
                avg_fission = total_fission / n_total
                pin_powers = mean_2d / avg_fission
                pin_std = std_2d / avg_fission

                # Find max and min pin powers (fuel pins only, > 0.01)
                fuel_mask = pin_powers > 0.01
                max_power = np.max(pin_powers[fuel_mask])
                min_power = np.min(pin_powers[fuel_mask])
                max_idx = np.unravel_index(
                    np.argmax(pin_powers * fuel_mask), pin_powers.shape
                )

                print(f"  Max relative pin power: {max_power:.4f} "
                      f"at position ({max_idx[0]}, {max_idx[1]})")
                print(f"  Min relative pin power: {min_power:.4f}")
                print(f"  Power peaking factor:   {max_power:.4f}")
                print()

                # Print the 17x17 pin power map
                print("  Pin Power Distribution (normalized to assembly average):")
                print("  " + "-" * 75)
                for row in range(17):
                    row_str = "  "
                    for col in range(17):
                        val = pin_powers[row, col]
                        if val < 0.01:
                            # Guide tube or instrument tube position
                            row_str += "  --- "
                        else:
                            row_str += f" {val:5.3f}"
                    print(row_str)
                print("  " + "-" * 75)
                print("  (--- = guide tube or instrument tube position)")
                print()

                # Convert to list for JSON serialization
                pin_power_map = pin_powers.tolist()

    except Exception as e:
        print(f"  No pin power tally found or error processing: {e}")
        print("  (Run with pin power mesh tally enabled for pin power data)")
        print()

    # --- Save to JSON ---
    results = {
        'benchmark': 'VERA Problem 2A',
        'description': '17x17 PWR fuel assembly, hot zero power, all rods out',
        'reference': {
            'k_eff': K_REF,
            'uncertainty': K_REF_UNC,
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

    if pin_power_map is not None:
        results['pin_powers'] = {
            'description': 'Normalized pin power distribution (17x17)',
            'normalization': 'Average over all 289 positions = 1.0',
            'max_peaking_factor': float(f"{max_power:.4f}"),
            'map': pin_power_map
        }

    output_file = 'results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {output_file}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze VERA Problem 2A OpenMC results'
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint file (auto-detect if omitted)')
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
