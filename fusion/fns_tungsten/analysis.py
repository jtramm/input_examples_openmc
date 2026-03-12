#!/usr/bin/env python3
"""
FNS Tungsten Benchmark -- Post-Processing & Analysis
=====================================================

This script reads the OpenMC statepoint file produced by the FNS tungsten
benchmark simulation and extracts:

  1. Neutron flux spectra at the 4 detector positions (0, 76, 228, 380 mm)
  2. Transmitted neutron spectrum through the rear face
  3. Total (energy-integrated) flux at each detector position
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
    """Extract energy bin centres and flux values from a flux/current tally.

    Parameters
    ----------
    tally : openmc.Tally
        A tally object with an EnergyFilter.

    Returns
    -------
    energy_centres : numpy.ndarray
        Geometric mean of each energy bin (eV).
    mean : numpy.ndarray
        Mean tally values (flux or current per source particle).
    std_dev : numpy.ndarray
        Standard deviation of tally values.
    """
    # Get the energy filter to find bin edges
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high] for each bin

    # Geometric mean of bin edges gives the representative energy
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])

    # Extract the tally data -- shape is (num_bins, 1, 1) for single cell/surface
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
    detector_depths_mm = [0, 76, 228, 380]  # mm from front face
    n_detectors = len(detector_depths_mm)

    # =========================================================================
    # Extract flux spectra at each detector position
    # =========================================================================
    print("\n--- Neutron Flux Spectra at Detector Positions ---")
    detector_spectra = {}

    for i in range(n_detectors):
        tally_name = f"flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        energy_centres, mean, std_dev = extract_spectrum(tally)

        depth_mm = detector_depths_mm[i]
        print(f"  Detector {i+1} (depth {depth_mm} mm): "
              f"peak flux = {mean.max():.4e} +/- {std_dev[mean.argmax()]:.4e}")

        # Store results (convert energy to MeV for readability)
        detector_spectra[f"detector_{i+1}_depth_{depth_mm}mm"] = {
            "energy_MeV": (energy_centres / 1.0e6).tolist(),
            "flux_mean": mean.tolist(),
            "flux_std_dev": std_dev.tolist(),
        }

    # =========================================================================
    # Extract transmitted spectrum through rear face
    # =========================================================================
    print("\n--- Transmitted Neutron Spectrum (Rear Face) ---")
    trans_tally = sp.get_tally(name="transmitted_spectrum")
    energy_centres, trans_mean, trans_std = extract_spectrum(trans_tally)

    # Total transmitted current (integral over all energies)
    total_transmitted = trans_mean.sum()
    print(f"  Total transmitted current: {total_transmitted:.4e} n/source-particle")
    print(f"  Peak energy bin: {energy_centres[trans_mean.argmax()]/1e6:.3f} MeV")

    transmitted_spectrum = {
        "energy_MeV": (energy_centres / 1.0e6).tolist(),
        "current_mean": trans_mean.tolist(),
        "current_std_dev": trans_std.tolist(),
        "total_transmitted": float(total_transmitted),
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
        depth_mm = detector_depths_mm[i]

        # Relative uncertainty (percent)
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  Detector {i+1} (depth {depth_mm:3d} mm): "
              f"flux = {mean_val:.4e} +/- {std_val:.4e} "
              f"({rel_unc:.1f}% rel. unc.)")

        total_fluxes[f"detector_{i+1}_depth_{depth_mm}mm"] = {
            "flux_mean": float(mean_val),
            "flux_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

    # =========================================================================
    # Compute attenuation through assembly
    # =========================================================================
    # Ratio of flux at deepest detector to flux at front face gives an
    # indication of overall attenuation through the tungsten assembly.
    flux_front = total_fluxes[f"detector_1_depth_0mm"]["flux_mean"]
    flux_deep = total_fluxes[f"detector_4_depth_380mm"]["flux_mean"]
    if flux_front > 0:
        attenuation_ratio = flux_deep / flux_front
        print(f"\n  Attenuation (380mm / 0mm): {attenuation_ratio:.4e}")
    else:
        attenuation_ratio = None
        print("\n  WARNING: Front-face flux is zero; cannot compute attenuation.")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "FNS Tungsten Clean Experiment",
        "source": "SINBAD NEA-1553/80",
        "facility": "FNS/JAEA, Tokai-mura, Japan",
        "runtime": runtime_info,
        "detector_spectra": detector_spectra,
        "transmitted_spectrum": transmitted_spectrum,
        "total_fluxes": total_fluxes,
        "attenuation_ratio_380mm_to_0mm": attenuation_ratio,
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
