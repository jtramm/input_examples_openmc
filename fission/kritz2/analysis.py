#!/usr/bin/env python3
"""
Analyze KRITZ-2 benchmark results.

Reads the OpenMC statepoint file, extracts k-eff, compares to the
benchmark reference value (average of 13 participant solutions from
NEA/NSC/DOC(2005)24), and writes a summary to results.json.

Usage:
    python analysis.py [--config 2:1-cold] [--statepoint statepoint.150.h5]
"""

import argparse
import json
import glob
import openmc

# Import the configuration data from the model script
from model import CONFIGS


def analyze(config_name, statepoint_path=None):
    """
    Compare OpenMC k-eff to the KRITZ-2 benchmark reference.

    Parameters
    ----------
    config_name : str
        Configuration key, e.g. '2:1-cold'
    statepoint_path : str or None
        Path to the statepoint HDF5 file.  If None, the most recent
        statepoint file in the current directory is used.

    Returns
    -------
    dict
        Results dictionary (also written to results.json)
    """
    cfg = CONFIGS[config_name]

    # Find the statepoint file if not specified
    if statepoint_path is None:
        candidates = sorted(glob.glob('statepoint.*.h5'))
        if not candidates:
            raise FileNotFoundError(
                "No statepoint files found.  Run OpenMC first."
            )
        statepoint_path = candidates[-1]
        print(f"Using statepoint: {statepoint_path}")

    # Open the statepoint and extract k-eff
    sp = openmc.StatePoint(statepoint_path)
    keff = sp.keff
    # keff is an openmc.stats.UFloat with nominal_value and std_dev
    k_calc = keff.nominal_value
    k_unc = keff.std_dev

    # Reference value from the benchmark (average of 13 participant solutions)
    k_ref = cfg['ref_keff']
    k_ref_unc = cfg['ref_keff_unc']

    # Compute the difference in pcm (per cent mille = 1e-5 dk)
    # A positive pcm means OpenMC predicts higher than the reference.
    delta_k = k_calc - k_ref
    delta_pcm = delta_k * 1e5  # Convert to pcm

    # Also compute difference from exact criticality (keff = 1.0)
    delta_from_critical_pcm = (k_calc - 1.0) * 1e5

    # Print summary
    print("\n" + "=" * 65)
    print(f"  KRITZ-2 Benchmark Results: {cfg['description']}")
    print("=" * 65)
    print(f"  Configuration:    KRITZ-{config_name}")
    print(f"  Fuel type:        {cfg['fuel_type']}")
    print(f"  Temperature:      {cfg['temperature_C']} C ({cfg['temperature_K']} K)")
    print(f"  Boron conc.:      {cfg['mod_boron_ppm']} ppm")
    print(f"  Lattice pitch:    {cfg['pitch']} cm")
    print("-" * 65)
    print(f"  OpenMC k-eff:     {k_calc:.5f} +/- {k_unc:.5f}")
    print(f"  Reference k-eff:  {k_ref:.5f}"
          + (f" +/- {k_ref_unc:.5f}" if k_ref_unc else " (unc. not reported)"))
    print("-" * 65)
    print(f"  Delta k-eff:      {delta_k:+.5f}")
    print(f"  Delta (pcm):      {delta_pcm:+.1f} pcm")
    print(f"  From critical:    {delta_from_critical_pcm:+.1f} pcm")
    print("-" * 65)

    # Interpret the result
    # The benchmark reference itself is not exactly 1.0 -- it is the
    # average of 13 calculations, most of which under-predict criticality.
    # So we compare to both the reference average and to exact criticality.
    if abs(delta_pcm) < 100:
        verdict = "EXCELLENT -- within 100 pcm of reference average"
    elif abs(delta_pcm) < 300:
        verdict = "GOOD -- within 300 pcm of reference average"
    elif abs(delta_pcm) < 500:
        verdict = "ACCEPTABLE -- within 500 pcm of reference average"
    else:
        verdict = "INVESTIGATE -- more than 500 pcm from reference average"
    print(f"  Assessment:       {verdict}")

    if k_ref_unc and abs(delta_k) < 2 * k_ref_unc:
        print(f"  Note: OpenMC result is within 2-sigma of the reference spread.")
    print("=" * 65)

    # Build results dictionary
    results = {
        'benchmark': 'KRITZ-2',
        'reference': 'NEA/NSC/DOC(2005)24',
        'configuration': f'KRITZ-{config_name}',
        'fuel_type': cfg['fuel_type'],
        'temperature_C': cfg['temperature_C'],
        'temperature_K': cfg['temperature_K'],
        'boron_ppm': cfg['mod_boron_ppm'],
        'pitch_cm': cfg['pitch'],
        'openmc_keff': k_calc,
        'openmc_keff_unc': k_unc,
        'reference_keff': k_ref,
        'reference_keff_unc': k_ref_unc,
        'delta_keff': delta_k,
        'delta_pcm': delta_pcm,
        'delta_from_critical_pcm': delta_from_critical_pcm,
        'statepoint': statepoint_path,
        'assessment': verdict,
    }

    # Write to JSON
    outfile = 'results.json'
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {outfile}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze KRITZ-2 benchmark results'
    )
    parser.add_argument(
        '--config', type=str, default='2:1-cold',
        choices=['2:1-cold', '2:1-hot', '2:13-cold', '2:13-hot',
                 '2:19-cold', '2:19-hot'],
        help='Configuration that was modelled (default: 2:1-cold)'
    )
    parser.add_argument(
        '--statepoint', type=str, default=None,
        help='Path to statepoint file (default: latest statepoint.*.h5)'
    )
    args = parser.parse_args()

    analyze(args.config, args.statepoint)


if __name__ == '__main__':
    main()
