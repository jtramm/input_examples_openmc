#!/usr/bin/env python3
"""
FNG-ITER Bulk Shield Mock-up -- Post-Processing & Analysis
============================================================

This script reads the OpenMC statepoint file produced by the FNG-ITER
bulk shield benchmark simulation and extracts:

  1. Neutron flux spectra at the 14 detector positions
  2. Total (energy-integrated) flux at each detector position
  3. Attenuation profile through the multi-layer shield
  4. Runtime statistics

Results are saved to results.json for comparison with experimental data
and MCNP reference calculations from the SINBAD database.

The key physics observable is the attenuation of neutron flux through
the ~94 cm mock-up, which should decrease by roughly 8-10 orders of
magnitude from front to rear. The different activation foil reactions
(Nb-93(n,2n), Al-27(n,alpha), Ni-58(n,p), Au-197(n,gamma)) probe
different parts of the neutron energy spectrum, providing detailed
validation of the transport calculation across energy groups.

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
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high] for each bin

    # Geometric mean of bin edges gives the representative energy
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
    # Detector positions (must match model.py)
    # =========================================================================
    # 14 detector positions at specific depths in the multi-layer shield
    detector_depths_cm = [
        3.43, 10.32, 17.15, 23.95, 30.80, 41.85, 46.85, 53.80,
        60.55, 67.40, 74.40, 81.10, 87.75, 92.15,
    ]
    n_detectors = len(detector_depths_cm)

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
        peak_idx = mean.argmax()
        print(f"  Detector {i+1:2d} (depth {depth_cm:6.2f} cm): "
              f"peak flux = {mean.max():.4e} +/- {std_dev[peak_idx]:.4e} "
              f"at {energy_centres[peak_idx]/1e6:.3f} MeV")

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
    print(f"  {'Det':>3s}  {'Depth (cm)':>10s}  {'Flux (n-cm/src)':>16s}  "
          f"{'Rel. Unc.':>10s}")
    print(f"  {'---':>3s}  {'----------':>10s}  {'---------------':>16s}  "
          f"{'----------':>10s}")

    total_fluxes = {}

    for i in range(n_detectors):
        tally_name = f"total_flux_detector_{i+1}"
        tally = sp.get_tally(name=tally_name)
        mean_val = tally.mean.flatten()[0]
        std_val = tally.std_dev.flatten()[0]
        depth_cm = detector_depths_cm[i]

        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  {i+1:3d}  {depth_cm:10.2f}  {mean_val:16.4e}  "
              f"{rel_unc:9.1f}%")

        total_fluxes[f"detector_{i+1}_depth_{depth_cm:.2f}cm"] = {
            "depth_cm": float(depth_cm),
            "flux_mean": float(mean_val),
            "flux_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

    # =========================================================================
    # Compute attenuation profile through the shield assembly
    # =========================================================================
    print("\n--- Attenuation Profile Through Shield ---")
    print("  (Flux ratio relative to first detector position)")

    flux_front_key = f"detector_1_depth_{detector_depths_cm[0]:.2f}cm"
    flux_front = total_fluxes[flux_front_key]["flux_mean"]
    attenuation_profile = {}

    if flux_front > 0:
        for i in range(n_detectors):
            depth_cm = detector_depths_cm[i]
            key = f"detector_{i+1}_depth_{depth_cm:.2f}cm"
            flux_val = total_fluxes[key]["flux_mean"]
            ratio = flux_val / flux_front

            # Compute approximate attenuation in decades
            decades = -np.log10(ratio) if ratio > 0 else float("inf")
            print(f"  Depth {depth_cm:6.2f} cm: ratio = {ratio:.4e}  "
                  f"({decades:.1f} decades attenuation)")

            attenuation_profile[key] = {
                "depth_cm": float(depth_cm),
                "flux_ratio_to_front": float(ratio),
                "attenuation_decades": float(decades),
            }
    else:
        print("  WARNING: Front detector flux is zero; cannot compute attenuation.")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "FNG-ITER Bulk Shield Mock-up",
        "source": "SINBAD Database",
        "facility": "FNG, ENEA Frascati, Italy",
        "description": (
            "Multi-layer mock-up (~94 cm) simulating the ITER inboard "
            "shielding: first wall (Cu), blanket (SS316/Perspex), vacuum "
            "vessel (SS316/Perspex), and TF coil (alternating Cu/SS316). "
            "Irradiated by 14.1 MeV D-T neutrons."
        ),
        "runtime": runtime_info,
        "detector_spectra": detector_spectra,
        "total_fluxes": total_fluxes,
        "attenuation_profile": attenuation_profile,
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
