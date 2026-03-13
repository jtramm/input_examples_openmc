"""
BEAVRS Full-Core Analysis Script
==================================

Post-processes OpenMC simulation results for the BEAVRS full-core benchmark.

This script:
1. Extracts the eigenvalue (k-eff) and compares to the measured critical value
2. Retrieves the assembly-wise fission rate distribution (proxy for power)
3. Normalizes the power map and prints it in a human-readable format
4. Computes radial peaking factors (Fxy)
5. Compares results to BEAVRS reference data

The BEAVRS benchmark measured critical condition at HZP, ARO, 975 ppm boron
gives k-eff ~ 0.99982. Monte Carlo codes with ENDF/B-VIII.0 data typically
predict k-eff within ~200 pcm of this value.

Usage:
    python analysis.py                           # Analyze latest statepoint
    python analysis.py --statepoint statepoint.100.h5  # Specific file
"""

import argparse
import glob
import os
import sys

import numpy as np

try:
    import openmc
except ImportError:
    print("ERROR: OpenMC Python API not found. Install with: pip install openmc")
    sys.exit(1)


# =============================================================================
# Reference values from BEAVRS Rev 2.0.1
# =============================================================================

# Measured critical eigenvalue at HZP, ARO, 975 ppm boron
# The reactor is essentially critical (k-eff = 1.0) but the benchmark
# reports a "measured" k-eff of 0.99982 after accounting for small
# reactivity biases from detector positioning and temperature gradients.
REFERENCE_KEFF = 0.99982

# 1 pcm = 1e-5 delta-k/k (percent mille)
# Typical Monte Carlo bias with ENDF/B-VIII.0: +100 to +300 pcm

# Core map for identifying which mesh cells contain fuel assemblies
# (reused from model.py)
CORE_MAP = [
    [0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0],
    [0, 0, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 0, 0],
    [0, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 0],
    [0, 3, 3, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 0],
    [3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3],
    [3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3],
    [3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
    [3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
    [3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3],
    [3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3],
    [3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3],
    [0, 3, 3, 2, 2, 2, 1, 1, 1, 2, 2, 2, 3, 3, 0],
    [0, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 0],
    [0, 0, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 0, 0],
    [0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0],
]


def find_latest_statepoint():
    """
    Find the most recent statepoint file in the current directory.

    OpenMC writes statepoint files as statepoint.<batch>.h5 where <batch>
    is the batch number. The final statepoint corresponds to the last batch.

    Returns
    -------
    str
        Path to the latest statepoint file

    Raises
    ------
    FileNotFoundError
        If no statepoint files are found
    """
    sp_files = sorted(glob.glob('statepoint.*.h5'))
    if not sp_files:
        raise FileNotFoundError(
            "No statepoint files found in the current directory.\n"
            "Run the simulation first: python model.py --run"
        )
    # Return the one with the highest batch number
    return sp_files[-1]


def analyze_keff(sp):
    """
    Extract and analyze the eigenvalue (k-effective) from the statepoint.

    The eigenvalue is the fundamental quantity in reactor criticality analysis.
    k-eff = 1.0 means the reactor is exactly critical (self-sustaining chain
    reaction). The BEAVRS benchmark should give k-eff very close to 1.0
    since it represents a measured critical configuration.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file

    Returns
    -------
    tuple
        (k_combined, k_std) -- combined k-eff estimate and its standard deviation
    """

    # OpenMC provides three k-eff estimators:
    #   - collision: based on collision rates in fissile material
    #   - track-length: based on track lengths through fissile material
    #   - absorption: based on absorption rates
    # The "combined" estimator optimally weights these for lowest variance.
    k_combined = sp.keff

    k_mean = k_combined.nominal_value
    k_std = k_combined.std_dev

    # Compute reactivity difference from reference in pcm
    # pcm = (k_calc - k_ref) / (k_calc * k_ref) * 1e5
    delta_k = k_mean - REFERENCE_KEFF
    reactivity_pcm = delta_k / (k_mean * REFERENCE_KEFF) * 1.0e5

    print("=" * 70)
    print("EIGENVALUE RESULTS")
    print("=" * 70)
    print(f"  k-eff (combined):     {k_mean:.5f} +/- {k_std:.5f}")
    print(f"  Reference k-eff:      {REFERENCE_KEFF:.5f}")
    print(f"  Delta k:              {delta_k:+.5f}")
    print(f"  Reactivity bias:      {reactivity_pcm:+.1f} pcm")
    print(f"  Statistical unc.:     {k_std * 1e5:.1f} pcm")
    print()

    # Interpret the result
    if abs(reactivity_pcm) < 200:
        print("  Assessment: EXCELLENT agreement with benchmark (< 200 pcm)")
    elif abs(reactivity_pcm) < 500:
        print("  Assessment: GOOD agreement with benchmark (< 500 pcm)")
    elif abs(reactivity_pcm) < 1000:
        print("  Assessment: FAIR agreement (500-1000 pcm) -- check model details")
    else:
        print("  Assessment: POOR agreement (> 1000 pcm) -- likely modeling issue")

    print()
    return k_mean, k_std


def analyze_assembly_power(sp):
    """
    Extract and analyze the assembly-wise power distribution.

    The fission rate tally on the 15x15 mesh gives a value proportional to
    power for each assembly position. We normalize so the average fuel-
    bearing assembly has a relative power of 1.0.

    Key metrics:
    - Radial peaking factor Fxy: ratio of maximum to average assembly power.
      For BEAVRS Cycle 1 at HZP, Fxy ~ 1.2-1.4 depending on the detailed
      loading pattern.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file

    Returns
    -------
    numpy.ndarray
        15x15 array of normalized assembly powers (0 for non-fuel positions)
    """

    # Find the assembly power tally
    tally = sp.get_tally(name='Assembly Power')

    # Extract the fission rate and reshape to the mesh dimensions
    # The mesh is 15x15, so we expect 225 values
    fission_rates = tally.get_values(scores=['fission']).reshape(15, 15)
    uncertainties = tally.get_values(scores=['fission'],
                                     value='std_dev').reshape(15, 15)

    # Create a boolean mask for fuel assembly positions
    core_mask = np.array(CORE_MAP) > 0

    # Normalize: set mean power of fuel-bearing assemblies to 1.0
    fuel_mean = fission_rates[core_mask].mean()
    if fuel_mean > 0:
        norm_power = fission_rates / fuel_mean
        norm_uncert = uncertainties / fuel_mean
    else:
        print("WARNING: No fission tallied. Check model and simulation.")
        return np.zeros((15, 15))

    # Zero out non-fuel positions for clarity
    norm_power[~core_mask] = 0.0
    norm_uncert[~core_mask] = 0.0

    # Compute peaking factor
    fxy = norm_power[core_mask].max()
    fxy_pos = np.unravel_index(np.argmax(norm_power), norm_power.shape)

    print("=" * 70)
    print("ASSEMBLY POWER DISTRIBUTION")
    print("=" * 70)
    print(f"  Radial peaking factor Fxy: {fxy:.4f}")
    print(f"  Peak assembly location:    row {fxy_pos[0]+1}, col {fxy_pos[1]+1}")
    print(f"  Min assembly power:        {norm_power[core_mask].min():.4f}")
    print()

    # Print the normalized power map
    print("  Normalized Assembly Power Map (rows x columns):")
    print("  (0.000 = water/no assembly, 1.000 = average assembly)")
    print()

    # Header row
    header = "      " + "".join(f"{j+1:>7d}" for j in range(15))
    print(header)
    print("      " + "-" * 105)

    for i in range(15):
        row_str = f"  {i+1:>2d} |"
        for j in range(15):
            if core_mask[i, j]:
                row_str += f" {norm_power[i, j]:5.3f} "
            else:
                row_str += "   -   "
        print(row_str)

    print()

    # Print enrichment zone averages
    for region in [1, 2, 3]:
        enrichments = {1: 1.6, 2: 2.4, 3: 3.1}
        region_mask = np.array(CORE_MAP) == region
        if np.any(region_mask):
            region_avg = norm_power[region_mask].mean()
            region_max = norm_power[region_mask].max()
            n_assy = np.sum(region_mask)
            print(f"  Region {region} ({enrichments[region]}% U-235, "
                  f"{n_assy} assemblies): "
                  f"avg power = {region_avg:.4f}, max = {region_max:.4f}")

    print()

    return norm_power


def analyze_symmetry(power_map):
    """
    Check the octant symmetry of the power distribution.

    A well-converged Monte Carlo simulation of a symmetric core should
    produce a nearly symmetric power distribution. Deviations indicate
    either insufficient statistics or modeling errors.

    Parameters
    ----------
    power_map : numpy.ndarray
        15x15 normalized assembly power map
    """

    print("=" * 70)
    print("SYMMETRY ANALYSIS")
    print("=" * 70)

    core_mask = np.array(CORE_MAP) > 0

    # Check left-right symmetry (flip about vertical axis, col 7)
    lr_flip = np.fliplr(power_map)
    lr_diff = np.abs(power_map - lr_flip)
    lr_max_diff = lr_diff[core_mask].max()
    lr_rms_diff = np.sqrt(np.mean(lr_diff[core_mask] ** 2))

    # Check top-bottom symmetry (flip about horizontal axis, row 7)
    tb_flip = np.flipud(power_map)
    tb_diff = np.abs(power_map - tb_flip)
    tb_max_diff = tb_diff[core_mask].max()
    tb_rms_diff = np.sqrt(np.mean(tb_diff[core_mask] ** 2))

    print(f"  Left-Right symmetry:")
    print(f"    Max difference: {lr_max_diff:.4f}")
    print(f"    RMS difference: {lr_rms_diff:.4f}")
    print(f"  Top-Bottom symmetry:")
    print(f"    Max difference: {tb_max_diff:.4f}")
    print(f"    RMS difference: {tb_rms_diff:.4f}")
    print()

    # Assess convergence based on symmetry
    if lr_rms_diff < 0.02 and tb_rms_diff < 0.02:
        print("  Assessment: Good statistical convergence (RMS < 2%)")
    elif lr_rms_diff < 0.05 and tb_rms_diff < 0.05:
        print("  Assessment: Moderate convergence -- consider more batches")
    else:
        print("  Assessment: Poor convergence -- significantly more histories needed")

    print()


def main():
    """Main analysis routine for BEAVRS benchmark results."""

    parser = argparse.ArgumentParser(
        description='Analyze BEAVRS full-core simulation results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script reads the statepoint file produced by an OpenMC eigenvalue
simulation of the BEAVRS benchmark and extracts key results:

  - k-effective (eigenvalue) compared to the measured critical value
  - Assembly power distribution and radial peaking factors
  - Symmetry check for convergence assessment

Examples:
    python analysis.py                                  # Latest statepoint
    python analysis.py --statepoint statepoint.200.h5   # Specific file
        """
    )
    parser.add_argument('--statepoint', type=str, default=None,
                        help='Path to statepoint HDF5 file (default: latest)')
    args = parser.parse_args()

    # Find and load the statepoint file
    if args.statepoint:
        sp_path = args.statepoint
    else:
        sp_path = find_latest_statepoint()

    print(f"\nLoading statepoint: {sp_path}")
    print()

    sp = openmc.StatePoint(sp_path)

    # Print simulation summary
    print("=" * 70)
    print("BEAVRS FULL-CORE BENCHMARK ANALYSIS")
    print("=" * 70)
    print(f"  Statepoint file:    {sp_path}")
    print(f"  Batches:            {sp.n_batches}")
    print(f"  Inactive batches:   {sp.n_inactive}")
    print(f"  Active batches:     {sp.n_batches - sp.n_inactive}")
    print(f"  Particles/batch:    {sp.n_particles}")
    total_active = (sp.n_batches - sp.n_inactive) * sp.n_particles
    print(f"  Total active hist.: {total_active:,}")
    print()

    # Analyze k-effective
    k_mean, k_std = analyze_keff(sp)

    # Analyze assembly power distribution
    power_map = analyze_assembly_power(sp)

    # Check symmetry
    if power_map.any():
        analyze_symmetry(power_map)

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    core_mask = np.array(CORE_MAP) > 0
    if power_map.any():
        fxy = power_map[core_mask].max()
        print(f"  k-eff = {k_mean:.5f} +/- {k_std:.5f}")
        print(f"  Fxy   = {fxy:.4f}")
        delta_pcm = (k_mean - REFERENCE_KEFF) / (k_mean * REFERENCE_KEFF) * 1e5
        print(f"  Bias  = {delta_pcm:+.0f} pcm vs measured critical")
    print("=" * 70)

    sp.close()


if __name__ == '__main__':
    main()
