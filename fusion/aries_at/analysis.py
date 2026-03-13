#!/usr/bin/env python3
"""
ARIES-AT Advanced Tokamak -- Post-Processing & Analysis
==========================================================

This script reads the OpenMC statepoint file produced by the ARIES-AT
benchmark simulation and extracts:

  1. Tritium Breeding Ratio (TBR) in blanket and back wall
  2. Nuclear heating in each component (first wall, blanket, shield, VV, TF)
  3. Energy multiplication factor
  4. Neutron flux spectrum at the TF coil
  5. Pb(n,2n) neutron multiplication rate
  6. Runtime statistics

Results are saved to results.json for comparison with published ARIES-AT
neutronics analyses.

Reference results (from Sawan & Abdou, Fusion Eng. Des., 2006):
  - TBR (with 90% Li-6): ~1.1
  - Energy multiplication: ~1.15
  - Blanket captures ~70% of fusion energy
  - Peak NWL: ~3.2 MW/m2 (outboard midplane)

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
    # TBR (Tritium Breeding Ratio)
    # =========================================================================
    print("\n--- Tritium Breeding Ratio ---")
    print("  (n,Xt) score: total triton production per source neutron")
    print("  Contributions from Li-6(n,t) + Li-7(n,n't) in Pb-17Li")
    print()

    tbr_tally = sp.get_tally(name="TBR")
    tbr_mean = tbr_tally.mean.flatten().sum()
    tbr_std = np.sqrt((tbr_tally.std_dev.flatten() ** 2).sum())
    tbr_rel = (tbr_std / tbr_mean * 100.0) if tbr_mean > 0 else 0.0

    tbr_blanket_tally = sp.get_tally(name="TBR_blanket")
    tbr_blanket_mean = tbr_blanket_tally.mean.flatten().sum()
    tbr_blanket_std = np.sqrt((tbr_blanket_tally.std_dev.flatten() ** 2).sum())

    tbr_backwall_tally = sp.get_tally(name="TBR_backwall")
    tbr_backwall_mean = tbr_backwall_tally.mean.flatten().sum()
    tbr_backwall_std = np.sqrt((tbr_backwall_tally.std_dev.flatten() ** 2).sum())

    print(f"  Total TBR:    {tbr_mean:.4f} +/- {tbr_std:.4f} ({tbr_rel:.1f}%)")
    print(f"  Blanket:      {tbr_blanket_mean:.4f} +/- {tbr_blanket_std:.4f}")
    print(f"  Back wall:    {tbr_backwall_mean:.4f} +/- {tbr_backwall_std:.4f}")
    print(f"  Reference:    ~1.1 (Sawan & Abdou 2006)")

    # =========================================================================
    # Nuclear heating
    # =========================================================================
    print("\n--- Nuclear Heating (eV per source neutron) ---")
    print("  'heating' score includes neutron + photon energy deposition")
    print()

    component_names = [
        "First wall", "Blanket", "Back wall",
        "Shield", "Vacuum vessel", "TF coil",
    ]

    heating_results = {}
    total_heating_sum = 0.0
    total_heating_var = 0.0

    for name in component_names:
        tally = sp.get_tally(name=f"heating_{name}")
        mean_val = tally.mean.flatten()[0]
        std_val = tally.std_dev.flatten()[0]
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        # Convert from eV to MeV for readability
        mean_mev = mean_val / 1.0e6
        std_mev = std_val / 1.0e6

        print(f"  {name:18s}: {mean_mev:.6e} MeV/src +/- {std_mev:.6e} "
              f"({rel_unc:.1f}%)")

        heating_results[name] = {
            "heating_eV_per_source": float(mean_val),
            "heating_std_eV": float(std_val),
            "heating_MeV_per_source": float(mean_mev),
            "relative_uncertainty_pct": float(rel_unc),
        }

        total_heating_sum += mean_val
        total_heating_var += std_val ** 2

    # Total heating from dedicated tally
    total_tally = sp.get_tally(name="heating_total")
    total_heating_mean = total_tally.mean.flatten()[0]
    total_heating_std = total_tally.std_dev.flatten()[0]

    print(f"\n  {'Total':18s}: {total_heating_mean/1e6:.6e} MeV/src "
          f"+/- {total_heating_std/1e6:.6e}")

    # Energy multiplication factor
    # E_mult = total deposited energy / source neutron energy (14.1 MeV)
    source_energy_ev = 14.1e6  # eV
    e_mult = total_heating_mean / source_energy_ev if total_heating_mean > 0 else 0.0
    e_mult_std = total_heating_std / source_energy_ev if total_heating_std > 0 else 0.0

    print(f"\n  Energy multiplication: {e_mult:.4f} +/- {e_mult_std:.4f}")
    print(f"  Reference: ~1.15 (Sawan & Abdou 2006)")

    # Blanket fraction of total heating
    blanket_heating = heating_results.get("Blanket", {}).get("heating_eV_per_source", 0)
    if total_heating_mean > 0:
        blanket_frac = blanket_heating / total_heating_mean * 100.0
        print(f"  Blanket fraction: {blanket_frac:.1f}% (reference: ~70%)")

    # =========================================================================
    # Pb(n,2n) neutron multiplication
    # =========================================================================
    print("\n--- Pb(n,2n) Neutron Multiplication ---")
    n2n_tally = sp.get_tally(name="Pb_n2n_blanket")
    n2n_mean = n2n_tally.mean.flatten()[0]
    n2n_std = n2n_tally.std_dev.flatten()[0]
    n2n_rel = (n2n_std / n2n_mean * 100.0) if n2n_mean > 0 else 0.0

    print(f"  (n,2n) rate: {n2n_mean:.6e} +/- {n2n_std:.6e} ({n2n_rel:.1f}%)")
    print("  This is the primary neutron multiplication mechanism in SiC/LiPb blankets.")

    # =========================================================================
    # TF coil neutron flux
    # =========================================================================
    print("\n--- Neutron Flux at TF Coil ---")

    # Total flux
    tf_total = sp.get_tally(name="flux_TF_coil_total")
    tf_flux_mean = tf_total.mean.flatten()[0]
    tf_flux_std = tf_total.std_dev.flatten()[0]
    tf_rel = (tf_flux_std / tf_flux_mean * 100.0) if tf_flux_mean > 0 else 0.0

    print(f"  Total flux: {tf_flux_mean:.6e} +/- {tf_flux_std:.6e} "
          f"({tf_rel:.1f}%) [n-cm/source]")

    # Spectrum
    tf_spec = sp.get_tally(name="flux_TF_coil")
    energy_centres, spec_mean, spec_std = extract_spectrum(tf_spec)

    if spec_mean.max() > 0:
        peak_idx = spec_mean.argmax()
        peak_energy_mev = energy_centres[peak_idx] / 1.0e6
        print(f"  Peak flux energy: {peak_energy_mev:.4f} MeV")
    else:
        print("  No flux recorded at TF coil (may need more particles or weight windows)")

    # =========================================================================
    # Blanket flux
    # =========================================================================
    print("\n--- Blanket Neutron Flux ---")
    blanket_flux_tally = sp.get_tally(name="flux_blanket")
    bl_flux_mean = blanket_flux_tally.mean.flatten()[0]
    bl_flux_std = blanket_flux_tally.std_dev.flatten()[0]
    bl_rel = (bl_flux_std / bl_flux_mean * 100.0) if bl_flux_mean > 0 else 0.0
    print(f"  Total flux: {bl_flux_mean:.6e} +/- {bl_flux_std:.6e} "
          f"({bl_rel:.1f}%) [n-cm/source]")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "ARIES-AT Advanced Tokamak Fusion Power Plant",
        "reference": "Najmabadi et al., Fusion Eng. Des. 80 (2006) 3-23",
        "description": (
            "Simplified CSG toroidal sector model of the ARIES-AT advanced "
            "tokamak with SiC/LiPb blanket. 22.5-degree sector with reflective "
            "boundaries, representing 1/16 of the full torus."
        ),
        "device_parameters": {
            "major_radius_cm": 520.0,
            "minor_radius_cm": 130.0,
            "fusion_power_MW": 1755,
            "NWL_outboard_MW_m2": 3.2,
            "n_TF_coils": 16,
        },
        "runtime": runtime_info,
        "TBR": {
            "total_mean": float(tbr_mean),
            "total_std_dev": float(tbr_std),
            "total_relative_uncertainty_pct": float(tbr_rel),
            "blanket_mean": float(tbr_blanket_mean),
            "blanket_std_dev": float(tbr_blanket_std),
            "backwall_mean": float(tbr_backwall_mean),
            "backwall_std_dev": float(tbr_backwall_std),
            "reference_value": 1.1,
        },
        "nuclear_heating": heating_results,
        "energy_multiplication": {
            "mean": float(e_mult),
            "std_dev": float(e_mult_std),
            "reference_value": 1.15,
        },
        "Pb_n2n_multiplication": {
            "rate_mean": float(n2n_mean),
            "rate_std_dev": float(n2n_std),
        },
        "TF_coil_flux": {
            "total_flux_mean": float(tf_flux_mean),
            "total_flux_std_dev": float(tf_flux_std),
            "spectrum_energy_MeV": (energy_centres / 1.0e6).tolist(),
            "spectrum_flux_mean": spec_mean.tolist(),
            "spectrum_flux_std_dev": spec_std.tolist(),
        },
        "blanket_flux": {
            "total_flux_mean": float(bl_flux_mean),
            "total_flux_std_dev": float(bl_flux_std),
        },
    }

    # Custom JSON encoder to handle numpy types
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
