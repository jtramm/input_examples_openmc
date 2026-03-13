#!/usr/bin/env python3
"""
Post-processing script for the NuScale-like SMR full-core benchmark.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, compares against the Serpent/SCALE reference solution,
and processes the assembly power distribution from the mesh tally.

Reference k-eff values from the McSAFER benchmark:
  HZP, ARO, fresh fuel: k-eff ~ 1.01-1.03 (from Serpent reference)
  SCALE/KENO-VI agrees within 44 pcm of Serpent

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

import numpy as np


# Reference k-eff for the full core (HZP, ARO, fresh fuel).
# The exact value depends on cross section library and whether
# burnable absorbers are modeled. Since this simplified model
# omits Gd2O3, k-eff will be somewhat higher than the full
# benchmark specification.
K_REF = 1.02
K_REF_UNC = 0.0005
K_REF_SOURCE = 'Serpent 2 / SCALE KENO-VI (McSAFER benchmark, approximate)'

# Core loading map for identifying assembly types in the output
CORE_MAP = [
    [0, 0, 3, 3, 3, 0, 0],
    [0, 3, 2, 2, 2, 3, 0],
    [3, 2, 1, 1, 1, 2, 3],
    [3, 2, 1, 1, 1, 2, 3],
    [3, 2, 1, 1, 1, 2, 3],
    [0, 3, 2, 2, 2, 3, 0],
    [0, 0, 3, 3, 3, 0, 0],
]


def find_statepoint():
    """Find the most recent statepoint file in the current directory."""
    files = glob.glob('statepoint.*.h5')
    if not files:
        raise FileNotFoundError(
            "No statepoint files found. Run 'python model.py --run' first."
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
        Results dictionary with k-eff, uncertainty, comparison metrics,
        and optionally assembly power distribution data.
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

    # --- Compute difference from reference ---
    delta_k = (k_mean - K_REF) * 1e5  # pcm
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
    tracking_rate = total_histories / total_time if total_time > 0 else 0.0

    # --- Print results summary ---
    print("\n" + "=" * 70)
    print("McSAFER NuScale-like SMR Full-Core Benchmark - Results")
    print("37 Fuel Assemblies, 3 Enrichment Zones, BORON-FREE, HZP, ARO")
    print("=" * 70)
    print(f"  k-eff (OpenMC):     {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  k-eff (reference):  {K_REF:.4f} +/- {K_REF_UNC:.5f}")
    print(f"  Difference:          {delta_k:+.1f} +/- {delta_k_unc:.1f} pcm")
    print(f"  Particles/batch:     {n_particles}")
    print(f"  Active batches:      {n_batches - n_inactive}")
    print(f"  Inactive batches:    {n_inactive}")
    print(f"  Total histories:     {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:       {total_time:.1f} s")
        print(f"  Tracking rate:       {tracking_rate:.0f} particles/s")
    print("=" * 70)

    # --- Assessment ---
    # Note: since this simplified model omits Gd2O3, k-eff will be
    # higher than the full benchmark specification. Comparison is
    # therefore approximate.
    if abs(delta_k) < 500:
        verdict = "EXCELLENT - within 500 pcm of reference"
    elif abs(delta_k) < 1000:
        verdict = "GOOD - within 1000 pcm of reference"
    elif abs(delta_k) < 2000:
        verdict = "FAIR - within 2000 pcm (expected: model omits Gd2O3)"
    else:
        verdict = "LARGE DIFFERENCE - check model setup and cross sections"
    print(f"  Assessment: {verdict}")
    print()
    print("  NOTES:")
    print("  - This model omits Gd2O3 burnable absorbers for simplicity")
    print("  - Without Gd, k-eff will be HIGHER than the full benchmark spec")
    print("  - The reference k-eff is approximate for this simplified case")
    print("  - SCALE/KENO-VI agrees within 44 pcm of Serpent for the full spec")
    print()

    # --- Assembly power distribution (if mesh tally exists) ---
    assy_powers_map = None
    max_peaking = None
    try:
        tally = sp.get_tally(name='Assembly Powers')
        if tally is not None:
            print("Processing assembly power distribution...")

            mean = tally.mean.ravel()
            mean_2d = mean.reshape((7, 7))

            # Normalize: average over the 37 fuel positions = 1.0
            fuel_mask = np.array(CORE_MAP) > 0
            fuel_total = np.sum(mean_2d[fuel_mask])
            if fuel_total > 0:
                avg_power = fuel_total / 37.0
                assy_powers = mean_2d / avg_power

                max_peaking = np.max(assy_powers[fuel_mask])
                min_peaking = np.min(assy_powers[fuel_mask])

                print(f"  Max assembly peaking: {max_peaking:.4f}")
                print(f"  Min assembly peaking: {min_peaking:.4f}")
                print()

                # Print 7x7 assembly power map
                enrichment_label = {0: '   ', 1: '1.6', 2: '2.4', 3: '3.1'}
                print("  Assembly Power Distribution (normalized, fuel positions only):")
                print("  " + "-" * 55)
                print("  Row   " + "".join(
                    f"  Col {j} " for j in range(7)
                ))
                print("  " + "-" * 55)
                for i in range(7):
                    row_str = f"   {i}  "
                    for j in range(7):
                        if CORE_MAP[i][j] == 0:
                            row_str += "   ---  "
                        else:
                            val = assy_powers[i, j]
                            row_str += f"  {val:5.3f} "
                    print(row_str)
                print("  " + "-" * 55)
                print("  (--- = empty water position)")
                print()

                # Print with enrichment labels
                print("  Assembly Types:")
                print("  " + "-" * 55)
                for i in range(7):
                    row_str = f"   {i}  "
                    for j in range(7):
                        label = enrichment_label[CORE_MAP[i][j]]
                        row_str += f"   {label}  "
                    print(row_str)
                print("  " + "-" * 55)
                print("  (1.6/2.4/3.1 = enrichment in wt% U-235)")
                print()

                assy_powers_map = assy_powers.tolist()

    except Exception as e:
        print(f"  No assembly power tally found or error: {e}")
        print()

    # --- Save results to JSON ---
    results = {
        'benchmark': 'McSAFER NuScale-like SMR Full Core',
        'description': (
            '37 fuel assemblies (17x17), 3 enrichment zones '
            '(1.6/2.4/3.1 wt% U-235), no Gd2O3, boron-free moderator, '
            'HZP, all rods out, fresh fuel'
        ),
        'design_features': {
            'reactor_type': 'NuScale-like Small Modular Reactor',
            'n_assemblies': 37,
            'assembly_type': '17x17 (264 fuel + 24 GT + 1 IT)',
            'enrichment_zones': {
                'type_1': {'enrichment_wt_pct': 1.6, 'count': 9,
                           'location': 'center'},
                'type_2': {'enrichment_wt_pct': 2.4, 'count': 12,
                           'location': 'middle ring'},
                'type_3': {'enrichment_wt_pct': 3.1, 'count': 16,
                           'location': 'outer ring'},
            },
            'moderator_boron_ppm': 0,
            'boron_free': True,
            'active_height_cm': 200.0,
            'power_mwth': 160,
            'pin_pitch_cm': 1.26,
            'assembly_pitch_cm': 21.50,
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
            'note': ('Simplified model without Gd2O3 will have higher k-eff '
                     'than full benchmark specification'),
        }
    }

    if total_time > 0:
        results['openmc']['runtime_seconds'] = float(f"{total_time:.1f}")
        results['openmc']['tracking_rate'] = float(f"{tracking_rate:.0f}")

    if assy_powers_map is not None:
        results['assembly_powers'] = {
            'description': 'Normalized assembly power distribution (7x7)',
            'normalization': 'Average over 37 fuel positions = 1.0',
            'max_peaking_factor': float(f"{max_peaking:.4f}"),
            'map': assy_powers_map,
            'core_map': CORE_MAP,
        }

    output_file = 'results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to {output_file}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze NuScale-like SMR full-core benchmark results'
    )
    parser.add_argument(
        '--statepoint', type=str, default=None,
        help='Path to statepoint file (auto-detect if omitted)'
    )
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
