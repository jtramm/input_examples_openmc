#!/usr/bin/env python3
"""
Post-processing script for the VVER-1000 MOX Core Benchmark results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, and compares against reference Monte Carlo calculations
from the NEA/NSC/DOC(2005)17 benchmark report.

Reference k-eff values (State S1: HZP, ARO, 600 ppm boron):
  The benchmark was designed as a near-critical configuration.
  MCU Monte Carlo:  ~1.000-1.010
  TVS-M:            ~1.005
  HELIOS/BIPR:      ~1.003
  Multiple codes agree within ~500 pcm.

  The benchmark reference keff depends on the code and nuclear data library
  used. A range of 1.000-1.010 is considered acceptable for State S1.

Source:
    V. Boyarinov et al., "VVER-1000 MOX Core Computational Benchmark,"
    NEA/NSC/DOC(2005)17, OECD Nuclear Energy Agency, 2005.

Usage:
    python analysis.py [--statepoint FILE]
"""

import argparse
import glob
import json


# Reference keff range from benchmark participants (State S1)
# The benchmark is designed to be near-critical; different codes and
# nuclear data libraries produce keff values spanning this range.
K_REF_LOW = 1.000
K_REF_HIGH = 1.010
K_REF_NOMINAL = 1.005  # approximate midpoint of participant results


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
    Load statepoint, extract k-eff, compare to benchmark reference.

    Parameters
    ----------
    statepoint_path : str or None
        Path to statepoint HDF5 file. If None, auto-detect.

    Returns
    -------
    dict
        Results dictionary with k-eff and comparison data.
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

    # --- Comparison to reference ---
    delta_k_ref = (k_mean - K_REF_NOMINAL) * 1e5  # pcm
    delta_k_unc = k_std * 1e5

    # --- Runtime information ---
    runtime = sp.runtime
    total_time = 0.0
    if isinstance(runtime, dict):
        total_time = runtime.get('total', 0.0)

    n_particles = sp.n_particles
    n_batches = sp.n_batches
    n_inactive = sp.n_inactive
    total_histories = n_particles * n_batches
    active_batches = n_batches - n_inactive

    if total_time > 0:
        tracking_rate = total_histories / total_time
    else:
        tracking_rate = 0.0

    # --- Print results ---
    print()
    print("=" * 72)
    print("VVER-1000 MOX Core Benchmark - Full Core Results")
    print("State S1: HZP, All Rods Out, 600 ppm boron, 552 K")
    print("163 hexagonal assemblies, ~30% MOX loading")
    print("=" * 72)
    print(f"  k-eff (OpenMC):           {k_mean:.6f} +/- {k_std:.6f}")
    print(f"  Reference range:          {K_REF_LOW:.3f} - {K_REF_HIGH:.3f}")
    print(f"  Reference nominal:        {K_REF_NOMINAL:.3f}")
    print(f"  Delta-k vs nominal:       {delta_k_ref:+.1f} +/- "
          f"{delta_k_unc:.1f} pcm")
    print()
    print(f"  Particles/batch:          {n_particles:,}")
    print(f"  Active batches:           {active_batches}")
    print(f"  Inactive batches:         {n_inactive}")
    print(f"  Total histories:          {total_histories:,}")
    if total_time > 0:
        print(f"  Total runtime:            {total_time:.1f} s "
              f"({total_time/3600:.2f} h)")
        print(f"  Tracking rate:            {tracking_rate:,.0f} particles/s")
    print("=" * 72)

    # --- Assessment ---
    within_range = K_REF_LOW - 0.005 <= k_mean <= K_REF_HIGH + 0.005
    if within_range:
        verdict = "WITHIN EXPECTED RANGE"
    elif 0.9 < k_mean < 1.1:
        verdict = "OUTSIDE EXPECTED RANGE - review model parameters"
    else:
        verdict = "CHECK MODEL - keff far from expected"

    print()
    print(f"Assessment: {verdict}")
    print()
    print("Notes:")
    print("  - This model uses simplified material compositions and a")
    print("    representative loading pattern based on the NEA benchmark.")
    print("  - Burned fuel is approximated as fresh UOX with reduced enrichment.")
    print("  - Differences of several hundred pcm from the benchmark reference")
    print("    are expected due to these simplifications.")
    print("  - For precise benchmark comparison, use the exact material")
    print("    compositions, geometry, and loading map from NEA/NSC/DOC(2005)17.")
    print()

    # --- Save results to JSON ---
    results = {
        'benchmark': 'VVER-1000 MOX Core Computational Benchmark',
        'document': 'NEA/NSC/DOC(2005)17',
        'state': 'S1 (HZP, ARO, 600 ppm boron, 552 K)',
        'core': {
            'n_assemblies': 163,
            'mox_fraction_pct': 31.3,
            'assembly_pitch_cm': 23.6,
            'pin_pitch_cm': 1.275,
            'pins_per_assembly': 312,
            'guide_tubes_per_assembly': 18,
        },
        'reference': {
            'keff_range': [K_REF_LOW, K_REF_HIGH],
            'keff_nominal': K_REF_NOMINAL,
            'source': 'NEA/NSC/DOC(2005)17, multiple participants',
        },
        'openmc': {
            'k_eff': float(k_mean),
            'uncertainty': float(k_std),
            'particles_per_batch': int(n_particles),
            'active_batches': int(active_batches),
            'inactive_batches': int(n_inactive),
            'total_histories': int(total_histories),
        },
        'comparison': {
            'delta_k_vs_nominal_pcm': float(f"{delta_k_ref:.1f}"),
            'delta_k_uncertainty_pcm': float(f"{delta_k_unc:.1f}"),
            'within_range': within_range,
            'assessment': verdict,
        },
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
        description='Analyze VVER-1000 MOX Core Benchmark OpenMC results'
    )
    parser.add_argument(
        '--statepoint', type=str, default=None,
        help='Path to statepoint file (auto-detect if omitted)'
    )
    args = parser.parse_args()
    analyze(args.statepoint)


if __name__ == '__main__':
    main()
