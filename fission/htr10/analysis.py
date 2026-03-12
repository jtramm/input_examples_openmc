#!/usr/bin/env python3
"""
Post-processing script for HTR-10 Pebble Bed Reactor results.

Loads the OpenMC statepoint file, extracts k-effective and runtime
information, compares against published reference solutions, and saves
results to results.json.

Reference k-eff values for HTR-10 initial criticality:
  - Experimental critical height: ~123 cm (measured by adding pebbles
    until criticality was reached)
  - Various code predictions for k-eff at ~126 cm loading height:
      MCNP:    1.0270 +/- 0.0020  (homogenized pebble model)
      VSOP:    1.0364             (diffusion + cell homogenization)
      Serpent: ~1.01-1.03         (varies with modeling detail)

  NOTE: The single-pebble "pin" model gives k-infinity of the fuel lattice,
  not the whole-core k-eff. k-inf for a TRISO pebble lattice is typically
  in the range 1.4-1.7 depending on enrichment and modeling approach.

Sources:
  [1] IAEA-TECDOC-1382 (2003)
  [2] Wu et al., Nuclear Engineering and Design 218 (2002) 25-32
  [3] Terry et al., HTR10-GCR-RESR-001, IRPhEP (2007)

Usage:
    python analysis.py [--statepoint FILE] [--model {pin,core}]
"""

import argparse
import glob
import json
import os


# Reference values for comparison
# For the pin (infinite lattice) model, k-inf is expected around 1.45-1.65
# depending on cross sections, homogenization, and enrichment treatment.
# For the core model, k-eff should be near 1.0 at the critical height.
REFERENCES = {
    'pin': {
        'k_eff': 1.55,       # approximate k-inf for HTR-10 pebble lattice
        'uncertainty': 0.05,  # wide band since this depends strongly on XS lib
        'source': 'Typical TRISO pebble k-inf (varies with modeling approach)',
        'document': 'IAEA-TECDOC-1382'
    },
    'core': {
        'k_eff': 1.027,
        'uncertainty': 0.002,
        'source': 'MCNP reference calculation, homogenized pebble model',
        'document': 'IAEA-TECDOC-1382; Terry et al. HTR10-GCR-RESR-001'
    }
}


def find_statepoint():
    """Find the most recent statepoint file in the current directory."""
    files = glob.glob('statepoint.*.h5')
    if not files:
        raise FileNotFoundError(
            "No statepoint files found. Run model.py --run first."
        )
    files.sort(key=lambda f: int(f.split('.')[1]))
    return files[-1]


def analyze(statepoint_path=None, model_type='pin'):
    """
    Load statepoint, extract results, compare to reference.

    Parameters
    ----------
    statepoint_path : str or None
        Path to statepoint HDF5 file. If None, auto-detect.
    model_type : str
        'pin' or 'core', determines which reference values to use.

    Returns
    -------
    dict
        Results dictionary.
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

    # --- Reference value ---
    ref = REFERENCES[model_type]
    k_ref = ref['k_eff']
    k_ref_unc = ref['uncertainty']

    # --- Compute difference ---
    delta_k = (k_mean - k_ref) * 1e5  # pcm
    delta_k_unc = (k_std**2 + k_ref_unc**2)**0.5 * 1e5

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
    model_desc = 'Single Pebble (k-inf)' if model_type == 'pin' else 'Full Core (k-eff)'
    print()
    print("=" * 60)
    print(f"HTR-10 Pebble Bed Reactor - {model_desc}")
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

    # Assessment
    if model_type == 'pin':
        # For k-inf, we expect larger variations
        if abs(delta_k) < 5000:
            verdict = "GOOD - within 5000 pcm of estimated k-inf"
        elif abs(delta_k) < 10000:
            verdict = "FAIR - within 10000 pcm of estimated k-inf"
        else:
            verdict = "CHECK - large deviation from estimated k-inf"
        print(f"  Assessment: {verdict}")
        print()
        print("  NOTE: The 'pin' model gives k-infinity of the infinite pebble")
        print("  lattice, not the critical k-eff of the HTR-10 core. Values")
        print("  of 1.4-1.7 are typical and depend on cross section library,")
        print("  temperature treatment, and whether TRISO is explicit or")
        print("  homogenized. The reference value is approximate.")
    else:
        if abs(delta_k) < 500:
            verdict = "EXCELLENT - within 500 pcm of reference"
        elif abs(delta_k) < 1000:
            verdict = "GOOD - within 1000 pcm of reference"
        elif abs(delta_k) < 2000:
            verdict = "FAIR - within 2000 pcm of reference"
        else:
            verdict = "POOR - more than 2000 pcm from reference"
        print(f"  Assessment: {verdict}")
        print()
        print("  NOTE: Differences from the reference may be due to:")
        print("    - Different cross section library")
        print("    - Fuel zone homogenization (loses TRISO self-shielding)")
        print("    - Simplified core geometry (no control rods, etc.)")
        print("    - Statistical uncertainty")
    print()

    # --- Save to JSON ---
    results = {
        'benchmark': 'HTR-10 Pebble Bed Reactor',
        'model_type': model_type,
        'description': model_desc,
        'reference': ref,
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
        description='Analyze HTR-10 OpenMC results'
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint file (auto-detect if omitted)')
    parser.add_argument('--model', type=str, default='pin',
                        choices=['pin', 'core'],
                        help='Model type for reference comparison')
    args = parser.parse_args()
    analyze(args.statepoint, args.model)


if __name__ == '__main__':
    main()
