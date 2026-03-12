#!/usr/bin/env python3
"""
FNG Copper Benchmark -- Post-Processing & Analysis
====================================================

This script reads the OpenMC statepoint file produced by the FNG copper
benchmark simulation and extracts:

  1. Neutron flux spectra at the 7 detector positions
  2. Total (energy-integrated) flux at each detector position
  3. Attenuation profile through the copper assembly
  4. Runtime statistics

Results are saved to results.json for comparison with experimental data
and MCNP reference calculations.

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
    # The last file (sorted lexicographically) has the highest batch number
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
    # Get the energy filter to find bin edges
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high] for each bin

    # Geometric mean of bin edges gives the representative energy
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])

    # Extract the tally data -- shape is (num_bins, 1, 1) for single cell
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
    # Detector positions (must match model.py)
    # =========================================================================
    # 7 detector positions at plate boundaries: 0, 9.99, 19.98, ... 59.94 cm
    plate_thickness = 9.99  # cm
    n_detectors = 7
    detector_depths_cm = [i * plate_thickness for i in range(n_detectors)]

    # =========================================================================
    # Extract flux spectra at each detector position
    # =========================================================================
    print("\n--- Neutron Flux Spectra at Detector Positions ---")
    detector_spectra = {}

    for i in range(n_detectors):
        tally_name = f"flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        energy_centres, mean, std_dev = extract_spectrum(tally)

        depth_cm = detector_depths_cm[i]
        print(f"  Detector {i+1} (depth {depth_cm:5.2f} cm): "
              f"peak flux = {mean.max():.4e} +/- {std_dev[mean.argmax()]:.4e}")

        # Store results (convert energy to MeV for readability)
        detector_spectra[f"detector_{i+1}_depth_{depth_cm:.2f}cm"] = {
            "depth_cm": float(depth_cm),
            "energy_MeV": (energy_centres / 1.0e6).tolist(),
            "flux_mean": mean.tolist(),
            "flux_std_dev": std_dev.tolist(),
        }

    # =========================================================================
    # Extract total (energy-integrated) flux at each detector
    # =========================================================================
    print("\n--- Total Flux at Detector Positions ---")
    total_fluxes = {}

    for i in range(n_detectors):
        tally_name = f"total_flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        mean_val = tally.mean.flatten()[0]
        std_val = tally.std_dev.flatten()[0]
        depth_cm = detector_depths_cm[i]

        # Relative uncertainty (percent)
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  Detector {i+1} (depth {depth_cm:5.2f} cm): "
              f"flux = {mean_val:.4e} +/- {std_val:.4e} "
              f"({rel_unc:.1f}% rel. unc.)")

        total_fluxes[f"detector_{i+1}_depth_{depth_cm:.2f}cm"] = {
            "depth_cm": float(depth_cm),
            "flux_mean": float(mean_val),
            "flux_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

    # =========================================================================
    # Compute attenuation profile through assembly
    # =========================================================================
    print("\n--- Attenuation Profile ---")

    # Reference flux is at the front face (detector 1, depth 0)
    flux_front = total_fluxes["detector_1_depth_0.00cm"]["flux_mean"]
    attenuation_profile = {}

    if flux_front > 0:
        for i in range(n_detectors):
            depth_cm = detector_depths_cm[i]
            key = f"detector_{i+1}_depth_{depth_cm:.2f}cm"
            flux_val = total_fluxes[key]["flux_mean"]
            ratio = flux_val / flux_front
            print(f"  Depth {depth_cm:5.2f} cm: flux ratio = {ratio:.4e}")
            attenuation_profile[key] = {
                "depth_cm": float(depth_cm),
                "flux_ratio_to_front": float(ratio),
            }
    else:
        print("  WARNING: Front-face flux is zero; cannot compute attenuation.")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "FNG Copper Benchmark Experiment",
        "source": "SINBAD Database",
        "facility": "FNG, ENEA Frascati, Italy",
        "runtime": runtime_info,
        "detector_spectra": detector_spectra,
        "total_fluxes": total_fluxes,
        "attenuation_profile": attenuation_profile,
    }

    # Custom JSON encoder to handle numpy types (int32, float64, etc.)
    class NumpyEncoder(json.JSONEncoder):
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
