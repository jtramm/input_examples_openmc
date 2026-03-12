#!/usr/bin/env python3
"""
TUD Iron Slab Benchmark -- Post-Processing & Analysis
=======================================================

This script reads the OpenMC statepoint file produced by the TUD iron slab
benchmark simulation and extracts:

  1. Transmitted neutron energy spectrum (surface current on rear face)
  2. Transmitted photon energy spectrum (surface current on rear face)
  3. Total transmitted neutron and photon currents (energy-integrated)
  4. Runtime statistics

Results are saved to results.json for comparison with experimental data
and MCNP reference calculations from SINBAD (FENDL-1 and EFF-2).

Usage:
    python analysis.py [--config A0|A1|A2]

Input:
    statepoint.*.h5  (latest statepoint file from OpenMC run)

Output:
    results.json     (structured results for benchmarking)
"""

import argparse
import glob
import json
import sys
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="TUD Iron Slab Benchmark -- analysis and post-processing"
)
parser.add_argument(
    "--config",
    type=str,
    choices=["A0", "A1", "A2"],
    default="A0",
    help="Configuration label for results metadata (default: A0)",
)
args = parser.parse_args()


def find_latest_statepoint():
    """Find the most recent statepoint file in the current directory.

    OpenMC writes statepoint files with the naming convention
    statepoint.<batch_number>.h5. This function finds the one with the
    highest batch number, which corresponds to the final state.

    Returns
    -------
    str
        Path to the latest statepoint file.
    """
    statepoint_files = sorted(glob.glob("statepoint.*.h5"))
    if not statepoint_files:
        print("ERROR: No statepoint files found. Run OpenMC first.")
        print("  Usage: python model.py --config A0 && openmc")
        sys.exit(1)
    # The last file (sorted lexicographically) has the highest batch number
    return statepoint_files[-1]


def extract_spectrum(tally):
    """Extract energy bin centres and tally values from a spectrum tally.

    Parameters
    ----------
    tally : openmc.Tally
        A tally object with an EnergyFilter.

    Returns
    -------
    energy_centres : numpy.ndarray
        Geometric mean of each energy bin (eV).
    mean : numpy.ndarray
        Mean tally values (current per source particle).
    std_dev : numpy.ndarray
        Standard deviation of tally values.
    """
    # Get the energy filter to find bin edges
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high] for each bin

    # Geometric mean of bin edges gives the representative energy for each bin
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])

    # Extract the tally data -- shape is (num_bins, 1, 1) for single surface
    mean = tally.mean.flatten()
    std_dev = tally.std_dev.flatten()

    return energy_centres, mean, std_dev


def main():
    # =========================================================================
    # Load the statepoint file
    # =========================================================================
    sp_path = find_latest_statepoint()
    print(f"Loading statepoint: {sp_path}")
    sp = openmc.StatePoint(sp_path)

    # Runtime information from the statepoint
    runtime_info = {
        "statepoint_file": sp_path,
        "n_batches": sp.n_batches,
        "particles_per_batch": sp.n_particles,
        "total_particles": sp.n_batches * sp.n_particles,
    }
    print(f"  Batches: {sp.n_batches}")
    print(f"  Particles/batch: {sp.n_particles:,}")
    print(f"  Total particles: {sp.n_batches * sp.n_particles:,}")

    # =========================================================================
    # Extract transmitted neutron spectrum
    # =========================================================================
    print("\n--- Transmitted Neutron Spectrum (Rear Face) ---")
    neutron_tally = sp.get_tally(name="transmitted_neutron_spectrum")
    n_energies, n_mean, n_std = extract_spectrum(neutron_tally)

    # Total transmitted neutron current (sum over all energy bins)
    total_neutron_current = n_mean.sum()
    print(f"  Total neutron current: {total_neutron_current:.4e} n/source-particle")

    # Peak energy bin (highest current contribution)
    peak_idx = n_mean.argmax()
    print(f"  Peak energy bin: {n_energies[peak_idx]/1e6:.3f} MeV "
          f"(current = {n_mean[peak_idx]:.4e})")

    # Store neutron spectrum results (convert energy to MeV for readability)
    neutron_spectrum = {
        "energy_MeV": (n_energies / 1.0e6).tolist(),
        "current_mean": n_mean.tolist(),
        "current_std_dev": n_std.tolist(),
        "total_current": float(total_neutron_current),
    }

    # =========================================================================
    # Extract transmitted photon spectrum
    # =========================================================================
    print("\n--- Transmitted Photon Spectrum (Rear Face) ---")
    photon_tally = sp.get_tally(name="transmitted_photon_spectrum")
    p_energies, p_mean, p_std = extract_spectrum(photon_tally)

    # Total transmitted photon current (sum over all energy bins)
    total_photon_current = p_mean.sum()
    print(f"  Total photon current: {total_photon_current:.4e} photons/source-particle")

    # Peak energy bin for photons
    if p_mean.max() > 0:
        p_peak_idx = p_mean.argmax()
        print(f"  Peak energy bin: {p_energies[p_peak_idx]/1e6:.3f} MeV "
              f"(current = {p_mean[p_peak_idx]:.4e})")
    else:
        print("  WARNING: No photons detected (check photon_transport setting)")

    # Store photon spectrum results
    photon_spectrum = {
        "energy_MeV": (p_energies / 1.0e6).tolist(),
        "current_mean": p_mean.tolist(),
        "current_std_dev": p_std.tolist(),
        "total_current": float(total_photon_current),
    }

    # =========================================================================
    # Extract total transmitted currents (energy-integrated tallies)
    # =========================================================================
    print("\n--- Total Transmitted Currents ---")

    # Total neutron current
    total_n_tally = sp.get_tally(name="total_transmitted_neutrons")
    total_n_mean = total_n_tally.mean.flatten()[0]
    total_n_std = total_n_tally.std_dev.flatten()[0]
    total_n_rel_unc = (total_n_std / total_n_mean * 100.0) if total_n_mean > 0 else 0.0
    print(f"  Neutrons: {total_n_mean:.4e} +/- {total_n_std:.4e} "
          f"({total_n_rel_unc:.2f}% rel. unc.)")

    # Total photon current
    total_p_tally = sp.get_tally(name="total_transmitted_photons")
    total_p_mean = total_p_tally.mean.flatten()[0]
    total_p_std = total_p_tally.std_dev.flatten()[0]
    total_p_rel_unc = (total_p_std / total_p_mean * 100.0) if total_p_mean > 0 else 0.0
    print(f"  Photons:  {total_p_mean:.4e} +/- {total_p_std:.4e} "
          f"({total_p_rel_unc:.2f}% rel. unc.)")

    # Neutron-to-photon ratio (informative for shielding analysis)
    if total_p_mean > 0:
        n_to_p_ratio = total_n_mean / total_p_mean
        print(f"  Neutron/photon ratio: {n_to_p_ratio:.2f}")
    else:
        n_to_p_ratio = None

    total_currents = {
        "neutron_current_mean": float(total_n_mean),
        "neutron_current_std_dev": float(total_n_std),
        "neutron_relative_uncertainty_pct": float(total_n_rel_unc),
        "photon_current_mean": float(total_p_mean),
        "photon_current_std_dev": float(total_p_std),
        "photon_relative_uncertainty_pct": float(total_p_rel_unc),
        "neutron_to_photon_ratio": n_to_p_ratio,
    }

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "TUD Iron Slab Experiment",
        "source": "SINBAD (Shielding Integral Benchmark Archive and Database)",
        "facility": "Technische Universitaet Dresden (TUD), Germany",
        "configuration": args.config,
        "url": "https://www.oecd-nea.org/science/wprs/shielding/sinbad/tud_fe/tufe-abs.htm",
        "runtime": runtime_info,
        "transmitted_neutron_spectrum": neutron_spectrum,
        "transmitted_photon_spectrum": photon_spectrum,
        "total_transmitted_currents": total_currents,
    }

    # Custom JSON encoder to handle numpy types (int32, float64, etc.)
    class NumpyEncoder(json.JSONEncoder):
        """JSON encoder that gracefully handles numpy numeric types."""
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    output_file = "results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)

    print(f"\nResults written to: {output_file}")
    print("Done.")


if __name__ == "__main__":
    main()
