#!/usr/bin/env python3
"""
Post-processing script for NuScale-like SMR benchmark results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, compares against the Serpent reference solution, and
processes the pin power distribution from the mesh tally.

Reference k-eff values from the McSAFER benchmark paper:
  Serpent 2 reference: ~1.1800 (fresh fuel, 3.1% enrichment, no Gd, ARO)
  OpenMC vs Serpent difference: ~355 pcm (reported in the paper)

Note: The exact reference value depends on the cross section library
(ENDF/B-VII.1 vs VIII.0) and the temperature treatment. The ~355 pcm
difference between OpenMC and Serpent is within typical code-to-code
agreement for these types of benchmarks.

Source:
  "Neutronics Benchmark of the NuScale-like SMR in the McSAFER Project"
  J. Bousquet et al., J. Nucl. Eng. 2025, 6(4), 44
  https://www.mdpi.com/2673-4362/6/4/44

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json
import os

import numpy as np


# Reference k-eff for the 3.1% enriched assembly (no Gd, all rods out)
# from the McSAFER benchmark Serpent 2 calculation. This is an approximate
# value; the benchmark paper reports results for multiple codes.
K_REF = 1.1800
K_REF_UNC = 0.00010
K_REF_SOURCE = 'Serpent 2 (McSAFER benchmark, approximate)'


def find_statepoint():
    """Find the most recent statepoint file in the current directory."""
    files = glob.glob('statepoint.*.h5')
    if not files:
        raise FileNotFoundError(
            "No statepoint files found. Run 'python model.py --run' first."
        )
    # Sort by batch number (extract integer from filename)
    files.sort(key=lambda f: int(f.split('.')[1]))
    return files[-1]


def analyze(statepoint_path=None):
    """
    Load statepoint, extract results, compare to Serpent reference.

    Parameters
    ----------
    statepoint_path : str or None
        Path to statepoint HDF5 file. If None, auto-detect.

    Returns
    -------
    dict
        Results dictionary with k-eff, uncertainty, comparison metrics,
        and optionally pin power distribution data.
    """
    import openmc

    if statepoint_path is None:
        statepoint_path = find_statepoint()

    print(f"Loading statepoint: {statepoint_path}")
    sp = openmc.StatePoint(statepoint_path)

    # --- Extract k-effective ---
    # OpenMC returns k-eff as an uncertainties.ufloat object with
    # nominal_value (best estimate) and std_dev (1-sigma uncertainty).
    keff = sp.keff
    k_mean = keff.nominal_value
    k_std = keff.std_dev

    # --- Compute difference from Serpent reference ---
    # Delta-k in pcm (percent mille): 1 pcm = 1e-5 in k-eff units.
    # A positive value means OpenMC is more reactive than the reference.
    delta_k = (k_mean - K_REF) * 1e5  # difference in pcm
    # Combined uncertainty: quadrature sum of OpenMC and reference uncertainties
    delta_k_unc = (k_std**2 + K_REF_UNC**2)**0.5 * 1e5

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

    # Particle tracking rate (histories per second)
    tracking_rate = total_histories / total_time if total_time > 0 else 0.0

    # --- Print results summary ---
    print("\n" + "=" * 70)
    print("NuScale-like SMR Benchmark - Results Summary")
    print("17x17 Fuel Assembly, 3.1% UO2, No Gd, BORON-FREE, All Rods Out")
    print("=" * 70)
    print(f"  k-eff (OpenMC):     {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  k-eff (Serpent ref): {K_REF:.4f} +/- {K_REF_UNC:.5f}")
    print(f"  Difference:          {delta_k:+.1f} +/- {delta_k_unc:.1f} pcm")
    print(f"  Particles/batch:     {n_particles}")
    print(f"  Active batches:      {n_batches - n_inactive}")
    print(f"  Inactive batches:    {n_inactive}")
    print(f"  Total histories:     {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:       {total_time:.1f} s")
        print(f"  Tracking rate:       {tracking_rate:.0f} particles/s")
    print("=" * 70)

    # Assess agreement with reference
    if abs(delta_k) < 200:
        verdict = "EXCELLENT - within 200 pcm of Serpent reference"
    elif abs(delta_k) < 500:
        verdict = "GOOD - within 500 pcm (expected code-to-code difference ~355 pcm)"
    elif abs(delta_k) < 1000:
        verdict = "FAIR - within 1000 pcm of Serpent reference"
    else:
        verdict = "LARGE BIAS - likely due to cross section temperature mismatch"
    print(f"  Assessment: {verdict}")
    print()

    # Note about expected biases
    print("  NOTES:")
    print("  - The McSAFER paper reports ~355 pcm OpenMC vs Serpent difference")
    print("  - Additional bias expected if cross section library only has 294K data")
    print("    while benchmark specifies 900K fuel / 557K moderator temperatures")
    print("  - The boron-free moderator makes the spectrum softer (more thermal)")
    print("    compared to conventional PWR benchmarks with ~1000 ppm boron")
    print()

    # --- Pin power distribution (if mesh tally exists) ---
    pin_powers = None
    pin_power_map = None
    max_power = None
    try:
        tally = sp.get_tally(name='Pin Powers')
        if tally is not None:
            print("Processing pin power distribution...")

            # Extract fission rate data and reshape to 17x17 grid
            mean = tally.mean.ravel()
            std_dev = tally.std_dev.ravel()
            mean_2d = mean.reshape((17, 17))
            std_2d = std_dev.reshape((17, 17))

            # Normalize to average = 1.0 over all 289 positions
            total_fission = np.sum(mean_2d)
            if total_fission > 0:
                n_total = 17 * 17
                avg_fission = total_fission / n_total
                pin_powers = mean_2d / avg_fission
                pin_std = std_2d / avg_fission

                # Statistics over fuel pin positions only (> 0.01 threshold
                # excludes guide tube and instrument tube positions which
                # have zero or near-zero fission rates)
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

                # Print 17x17 pin power map
                print("  Pin Power Distribution (normalized to assembly average):")
                print("  " + "-" * 75)
                for row in range(17):
                    row_str = "  "
                    for col in range(17):
                        val = pin_powers[row, col]
                        if val < 0.01:
                            row_str += "  --- "
                        else:
                            row_str += f" {val:5.3f}"
                    print(row_str)
                print("  " + "-" * 75)
                print("  (--- = guide tube or instrument tube position)")
                print()

                pin_power_map = pin_powers.tolist()

    except Exception as e:
        print(f"  No pin power tally found or error processing: {e}")
        print("  (Run the assembly model with pin power mesh tally for pin data)")
        print()

    # --- Save results to JSON ---
    results = {
        'benchmark': 'McSAFER NuScale-like SMR',
        'description': ('17x17 fuel assembly, 3.1% UO2, no gadolinium, '
                        'boron-free moderator, all rods out'),
        'design_features': {
            'reactor_type': 'NuScale-like Small Modular Reactor',
            'assembly_type': '17x17 (264 fuel + 24 GT + 1 IT)',
            'enrichment_wt_pct': 3.1,
            'moderator_boron_ppm': 0,
            'boron_free': True,
            'burnable_absorber': 'None (simplest assembly type)',
            'pin_pitch_cm': 1.2598,
            'assembly_pitch_cm': 21.5036,
        },
        'reference': {
            'k_eff': K_REF,
            'uncertainty': K_REF_UNC,
            'source': K_REF_SOURCE,
            'paper': 'J. Nucl. Eng. 2025, 6(4), 44',
            'url': 'https://www.mdpi.com/2673-4362/6/4/44',
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
            'assessment': verdict,
            'note': 'McSAFER paper reports ~355 pcm OpenMC vs Serpent difference',
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
        description='Analyze NuScale-like SMR benchmark OpenMC results'
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint file (auto-detect if omitted)')
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
