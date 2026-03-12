#!/usr/bin/env python3
"""
Post-processing script for BEAVRS single fuel assembly results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, and provides context for interpreting the results against
the BEAVRS benchmark measured data.

BEAVRS Reference Data (Real Plant Measurements):
  - HZP critical boron concentration: ~975 ppm (measured)
  - At this boron level, the FULL CORE k-eff = 1.0 (critical)
  - Our single-assembly infinite-lattice result will NOT be 1.0
    because it lacks leakage and inter-assembly heterogeneity effects

Interpreting Results:
  Unlike the VERA benchmark where we compare against a computed reference
  k-eff for the same geometry, BEAVRS provides real plant measurements.
  The measured quantity is the critical boron concentration, not k-eff
  for a single assembly. Therefore:

  1. Our k-eff will be > 1.0 at 975 ppm because:
     - No radial or axial leakage (reflective BCs)
     - No inter-assembly effects
     - All assemblies assumed to be the same enrichment

  2. To compare meaningfully with BEAVRS, one would need a full-core
     model and iterate on boron concentration to find k=1.0.

  3. We can still compare with other codes' single-assembly results
     at the same conditions (enrichment, boron, temperature).

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json
import os

import numpy as np


# BEAVRS measured critical boron concentration at HZP
# This is the primary benchmark quantity - NOT a computed k-eff
BEAVRS_CRITICAL_BORON_PPM = 975
# Note: there is no single-assembly reference k-eff for BEAVRS like
# there is for VERA. We provide context instead of a direct comparison.


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
    Load statepoint, extract results, and provide BEAVRS context.

    Parameters
    ----------
    statepoint_path : str or None
        Path to statepoint HDF5 file. If None, auto-detect.

    Returns
    -------
    dict
        Results dictionary with k-eff, uncertainty, and context.
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

    # --- Print results ---
    print("\n" + "=" * 70)
    print("BEAVRS Benchmark - Single Fuel Assembly Results")
    print("17x17 PWR Assembly, Hot Zero Power, All Rods Out, 975 ppm Boron")
    print("=" * 70)
    print(f"  k-eff (OpenMC):      {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  Particles/batch:     {n_particles}")
    print(f"  Active batches:      {n_batches - n_inactive}")
    print(f"  Inactive batches:    {n_inactive}")
    print(f"  Total histories:     {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:       {total_time:.1f} s")
        print(f"  Tracking rate:       {tracking_rate:.0f} particles/s")
    print("=" * 70)

    # --- Interpretation ---
    print()
    print("  INTERPRETATION:")
    print("  ---------------")
    print(f"  The BEAVRS HZP measured critical boron concentration is "
          f"~{BEAVRS_CRITICAL_BORON_PPM} ppm.")
    print("  At this boron level, the FULL CORE is critical (k=1.0).")
    print()
    print("  Our single-assembly infinite-lattice k-eff is expected to be")
    print("  significantly above 1.0 because:")
    print("    1. No radial leakage (reflective BCs simulate infinite lattice)")
    print("    2. No axial leakage (2D model with reflective axial BCs)")
    print("    3. All assemblies same enrichment (real core has mixed loading)")
    print("    4. The 975 ppm boron was measured for the real heterogeneous core")
    print()

    # Estimate excess reactivity
    rho = (k_mean - 1.0) / k_mean * 1e5  # reactivity in pcm
    print(f"  Excess reactivity:   {rho:+.0f} pcm")
    print(f"  (This excess is expected for an infinite lattice of 3.1% fuel)")
    print()

    # Note about cross section temperature
    print("  NOTE ON TEMPERATURE BIAS:")
    print("  If the cross section library only has 294K data while BEAVRS")
    print("  specifies 565K/600K, there will be a positive bias in k-eff")
    print("  (~1000-1500 pcm) due to missing Doppler broadening of U-238")
    print("  resonances. The Doppler effect increases parasitic absorption")
    print("  at higher temperatures, reducing k-eff.")
    print()

    # --- Pin power distribution ---
    pin_powers = None
    pin_power_map = None
    max_power = None
    try:
        tally = sp.get_tally(name='Pin Powers')
        if tally is not None:
            print("Processing pin power distribution...")

            mean = tally.mean.ravel()
            std_dev = tally.std_dev.ravel()

            # Reshape to 17x17 grid
            mean_2d = mean.reshape((17, 17))
            std_2d = std_dev.reshape((17, 17))

            # Normalize to average over all positions = 1.0
            total_fission = np.sum(mean_2d)
            if total_fission > 0:
                n_total = 17 * 17
                avg_fission = total_fission / n_total
                pin_powers = mean_2d / avg_fission
                pin_std = std_2d / avg_fission

                # Find max and min pin powers (fuel pins only)
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
        print("  (Run with pin power mesh tally enabled for pin power data)")
        print()

    # --- Save to JSON ---
    results = {
        'benchmark': 'BEAVRS Single Fuel Assembly',
        'description': (
            '17x17 PWR fuel assembly from the BEAVRS benchmark, '
            'hot zero power, all rods out, 975 ppm soluble boron'
        ),
        'reference': {
            'type': 'measured plant data',
            'critical_boron_ppm': BEAVRS_CRITICAL_BORON_PPM,
            'note': (
                'BEAVRS provides measured critical boron concentration for '
                'the full core, not a single-assembly reference k-eff. '
                'Direct k-eff comparison is not applicable for this model.'
            ),
            'source': 'MIT Computational Reactor Physics Group',
            'document': 'BEAVRS v1.0.1 specification'
        },
        'openmc': {
            'k_eff': float(k_mean),
            'uncertainty': float(k_std),
            'excess_reactivity_pcm': float(f"{rho:.0f}"),
            'particles_per_batch': int(n_particles),
            'active_batches': int(n_batches - n_inactive),
            'inactive_batches': int(n_inactive),
            'total_histories': int(total_histories),
        },
        'model_parameters': {
            'enrichment_wpct': 3.1,
            'boron_ppm': 975,
            'fuel_temperature_K': 600,
            'moderator_temperature_K': 565,
            'moderator_density_gcc': 0.740,
            'boundary_conditions': 'reflective (infinite lattice)',
            'geometry': '2D (reflective axial BCs)'
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
        description='Analyze BEAVRS single fuel assembly OpenMC results'
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint file (auto-detect if omitted)')
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
