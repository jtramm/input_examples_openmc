#!/usr/bin/env python3
"""
FNG-ITER Streaming Experiment -- Post-Processing & Analysis
=============================================================

This script reads the OpenMC statepoint file produced by the FNG-ITER
streaming benchmark and extracts:

  1. Total neutron flux at each detector position along the channel,
     inside the cavity, and behind the assembly
  2. Streaming ratio: flux inside the channel vs flux in surrounding
     bulk material at the same depth (quantifies streaming enhancement)
  3. Flux attenuation profile through the assembly
  4. Runtime statistics

Results are saved to results.json.

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


# =============================================================================
# Detector position labels -- must match model.py exactly
# =============================================================================
channel_det_z = [0.25, 12.95, 25.95, 38.65]
cavity_det_z = np.linspace(39.12, 43.82, 11).tolist()
behind_det_z = [46.35, 53.30, 60.05, 66.90, 73.90, 80.60, 87.25, 91.65]

all_det_z = channel_det_z + cavity_det_z + behind_det_z
all_det_labels = (
    [f"channel_{z:.2f}cm" for z in channel_det_z]
    + [f"cavity_{z:.2f}cm" for z in cavity_det_z]
    + [f"behind_{z:.2f}cm" for z in behind_det_z]
)


def find_latest_statepoint():
    """Find the most recent statepoint file in the current directory.

    OpenMC writes statepoint files named statepoint.<batch>.h5.
    This function returns the one with the highest batch number.

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


def main():
    # =========================================================================
    # Load statepoint
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
    # Extract total flux at each detector position
    # =========================================================================
    print("\n--- Total Flux at Detector Positions ---")
    print(f"{'Position':<25s} {'Depth (cm)':>10s} {'Flux':>12s} "
          f"{'Std Dev':>12s} {'Rel Unc (%)':>12s}")
    print("-" * 75)

    detector_fluxes = {}
    for i, (z_pos, label) in enumerate(zip(all_det_z, all_det_labels)):
        tally_name = f"total_flux_{label}"
        tally = sp.get_tally(name=tally_name)
        mean_val = float(tally.mean.flatten()[0])
        std_val = float(tally.std_dev.flatten()[0])
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        region = "channel" if z_pos <= 39.07 else ("cavity" if z_pos <= 44.53 else "behind")
        print(f"  {label:<23s} {z_pos:>10.2f} {mean_val:>12.4e} "
              f"{std_val:>12.4e} {rel_unc:>11.1f}%")

        detector_fluxes[label] = {
            "depth_cm": float(z_pos),
            "region": region,
            "flux_mean": mean_val,
            "flux_std_dev": std_val,
            "relative_uncertainty_pct": rel_unc,
        }

    # =========================================================================
    # Compute streaming ratio
    # =========================================================================
    # The streaming ratio compares the flux at positions along the channel
    # to the flux at the first behind-cavity position (which represents
    # flux through bulk material). A high ratio at channel positions
    # indicates strong streaming enhancement.
    print("\n--- Streaming Ratio ---")
    print("  (Ratio of flux at each position to first behind-cavity detector)")

    # Use first behind-cavity position as the bulk-material reference
    ref_label = all_det_labels[len(channel_det_z) + len(cavity_det_z)]  # behind_46.35cm
    ref_flux = detector_fluxes[ref_label]["flux_mean"]

    streaming_ratios = {}
    if ref_flux > 0:
        for label in all_det_labels:
            flux = detector_fluxes[label]["flux_mean"]
            ratio = flux / ref_flux
            depth = detector_fluxes[label]["depth_cm"]
            print(f"  {label:<23s} depth={depth:6.2f} cm  "
                  f"ratio = {ratio:.4e}")
            streaming_ratios[label] = {
                "depth_cm": depth,
                "streaming_ratio": ratio,
            }
    else:
        print("  WARNING: Reference flux is zero; cannot compute streaming ratio.")
        print("  This is expected if statistics are very poor (streaming problem).")

    # =========================================================================
    # Flux attenuation profile
    # =========================================================================
    print("\n--- Attenuation Profile (normalised to front detector) ---")
    front_label = all_det_labels[0]
    front_flux = detector_fluxes[front_label]["flux_mean"]

    attenuation = {}
    if front_flux > 0:
        for label in all_det_labels:
            flux = detector_fluxes[label]["flux_mean"]
            ratio = flux / front_flux
            depth = detector_fluxes[label]["depth_cm"]
            print(f"  {label:<23s} depth={depth:6.2f} cm  "
                  f"flux/flux_front = {ratio:.4e}")
            attenuation[label] = {
                "depth_cm": depth,
                "flux_ratio_to_front": ratio,
            }
    else:
        print("  WARNING: Front detector flux is zero.")

    # =========================================================================
    # Write results.json
    # =========================================================================
    results = {
        "benchmark": "FNG-ITER Streaming Experiment",
        "source": "SINBAD Database",
        "facility": "FNG, ENEA Frascati, Italy",
        "description": "Neutron streaming through 28 mm duct in ITER-like "
                       "shielding assembly with detector cavity",
        "runtime": runtime_info,
        "detector_fluxes": detector_fluxes,
        "streaming_ratios": streaming_ratios,
        "attenuation_profile": attenuation,
    }

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
