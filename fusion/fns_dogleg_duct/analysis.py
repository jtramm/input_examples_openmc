#!/usr/bin/env python3
"""
FNS Dogleg Duct Streaming Experiment -- Post-Processing & Analysis
====================================================================

This script reads the OpenMC statepoint file produced by the FNS dogleg duct
benchmark simulation and extracts:

  1. Neutron flux spectra at the 4 detector positions within the duct
  2. Total (energy-integrated) flux at each detector position
  3. Streaming attenuation factors between detector positions
  4. Runtime statistics

Results are saved to results.json for comparison with experimental data.

The key physics observable is the attenuation of the neutron flux as particles
navigate around the two right-angle bends of the dogleg duct. The ratio of
flux at downstream positions to the first detector quantifies how effectively
the duct geometry attenuates the streaming radiation.

Usage:
    python analysis.py

Input:
    statepoint.*.h5  (latest statepoint file from OpenMC run)

Output:
    results.json     (structured results for benchmarking)
"""

import glob
import json
import sys
import numpy as np
import openmc


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
        sys.exit(1)
    return statepoint_files[-1]


def extract_spectrum(tally):
    """Extract energy bin centres and flux values from a flux tally.

    Parameters
    ----------
    tally : openmc.Tally
        A tally object with an EnergyFilter.

    Returns
    -------
    energy_centres : numpy.ndarray
        Geometric mean of each energy bin (eV).
    mean : numpy.ndarray
        Mean tally values (flux per source particle).
    std_dev : numpy.ndarray
        Standard deviation of tally values.
    """
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high]
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])

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
    # Detector descriptions (must match model.py)
    # =========================================================================
    detector_labels = [
        "Position 1: middle of 1st horizontal leg",
        "Position 2: bottom bend (1st leg -> vertical leg)",
        "Position 3: top bend (vertical leg -> 2nd leg)",
        "Position 4: middle of 2nd horizontal leg",
    ]
    n_detectors = len(detector_labels)

    # =========================================================================
    # Extract flux spectra at each detector position
    # =========================================================================
    print("\n--- Neutron Flux Spectra at Detector Positions ---")
    detector_spectra = {}

    for i in range(n_detectors):
        tally_name = f"flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        energy_centres, mean, std_dev = extract_spectrum(tally)

        peak_idx = mean.argmax()
        print(f"  {detector_labels[i]}:")
        print(f"    Peak flux = {mean[peak_idx]:.4e} +/- {std_dev[peak_idx]:.4e} "
              f"at {energy_centres[peak_idx]/1e6:.3f} MeV")

        detector_spectra[f"detector_{i+1}"] = {
            "description": detector_labels[i],
            "energy_MeV": (energy_centres / 1.0e6).tolist(),
            "flux_mean": mean.tolist(),
            "flux_std_dev": std_dev.tolist(),
        }

    # =========================================================================
    # Extract total (energy-integrated) flux at each detector
    # =========================================================================
    print("\n--- Total Flux at Detector Positions ---")
    total_fluxes = {}
    flux_values = []  # for attenuation calculation

    for i in range(n_detectors):
        tally_name = f"total_flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        mean_val = tally.mean.flatten()[0]
        std_val = tally.std_dev.flatten()[0]

        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  {detector_labels[i]}:")
        print(f"    Flux = {mean_val:.4e} +/- {std_val:.4e} "
              f"({rel_unc:.1f}% rel. unc.)")

        total_fluxes[f"detector_{i+1}"] = {
            "description": detector_labels[i],
            "flux_mean": float(mean_val),
            "flux_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }
        flux_values.append(mean_val)

    # =========================================================================
    # Compute streaming attenuation factors
    # =========================================================================
    # The attenuation factor is the ratio of flux at each downstream position
    # to the flux at the first detector. This quantifies how much the dogleg
    # geometry reduces the streaming radiation.
    print("\n--- Streaming Attenuation (relative to Detector 1) ---")
    attenuation = {}

    if flux_values[0] > 0:
        for i in range(1, n_detectors):
            ratio = flux_values[i] / flux_values[0]
            print(f"  Detector {i+1} / Detector 1: {ratio:.4e}")
            attenuation[f"detector_{i+1}_over_detector_1"] = float(ratio)

        # Overall attenuation: last detector / first detector
        overall = flux_values[-1] / flux_values[0]
        print(f"\n  Overall attenuation (Det 4 / Det 1): {overall:.4e}")
        attenuation["overall_det4_over_det1"] = float(overall)
    else:
        print("  WARNING: Detector 1 flux is zero; cannot compute attenuation.")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "FNS Dogleg Duct Streaming Experiment",
        "source": "SINBAD",
        "facility": "FNS/JAEA, Tokai-mura, Japan",
        "runtime": runtime_info,
        "detector_spectra": detector_spectra,
        "total_fluxes": total_fluxes,
        "streaming_attenuation": attenuation,
    }

    class NumpyEncoder(json.JSONEncoder):
        """JSON encoder that handles numpy types."""
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
