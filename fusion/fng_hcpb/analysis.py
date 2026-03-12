#!/usr/bin/env python3
"""
FNG HCPB Tritium Breeder Module Mock-up -- Post-Processing & Analysis
=======================================================================

This script reads the OpenMC statepoint file produced by the FNG HCPB
benchmark simulation and extracts:

  1. Tritium production rates in each Li2CO3 breeder layer
  2. Total tritium production across all breeder layers (TBR indicator)
  3. Neutron flux spectra in breeder layers (to assess spectrum softening)
  4. Be-9(n,2n) reaction rates in beryllium zones (multiplication efficiency)
  5. Neutron flux in beryllium zones
  6. Runtime statistics

Results are saved to results.json for comparison with experimental data
and MCNP-4C reference calculations from the SINBAD database.

The key experimental observables are:
  - Tritium production rate from Li-6(n,t) at 16 pellet positions
  - Activation foil reaction rates at 4 depths (y=4.2, 10.5, 16.8, 23.1 cm)
  - Total tritium production in each breeder section

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
    # Breeder layer names (must match model.py)
    # =========================================================================
    breeder_names = [
        "breeder_1a", "breeder_1b",
        "breeder_2a", "breeder_2b",
        "rear_cassette",
    ]
    breeder_descriptions = [
        "First double-layer, front (y ~ 3.65-4.85 cm)",
        "First double-layer, rear (y ~ 4.95-6.15 cm)",
        "Second double-layer, front (y ~ 16.05-17.25 cm)",
        "Second double-layer, rear (y ~ 17.35-18.55 cm)",
        "Rear cassette (y ~ 29.5-43.3 cm)",
    ]

    # =========================================================================
    # Extract tritium production rates in each breeder layer
    # =========================================================================
    print("\n--- Tritium Production Rates (per source neutron) ---")
    print("  (n,Xt) score: total triton production from Li-6(n,t) + Li-7(n,n't)")
    print()

    tritium_results = {}
    total_tritium_mean = 0.0
    total_tritium_var = 0.0

    for bname, bdesc in zip(breeder_names, breeder_descriptions):
        tally_name = f"tritium_production_{bname}"
        tally = sp.get_tally(name=tally_name)
        mean_val = tally.mean.flatten()[0]
        std_val = tally.std_dev.flatten()[0]
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  {bname:20s}: {mean_val:.6e} +/- {std_val:.6e} "
              f"({rel_unc:.1f}% rel. unc.)  [{bdesc}]")

        tritium_results[bname] = {
            "description": bdesc,
            "tritium_production_mean": float(mean_val),
            "tritium_production_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

        total_tritium_mean += mean_val
        total_tritium_var += std_val ** 2

    # Total tritium from the dedicated tally
    total_tally = sp.get_tally(name="total_tritium_production")
    total_mean = total_tally.mean.flatten().sum()
    total_std = np.sqrt((total_tally.std_dev.flatten() ** 2).sum())
    total_rel_unc = (total_std / total_mean * 100.0) if total_mean > 0 else 0.0

    print(f"\n  {'TOTAL':20s}: {total_mean:.6e} +/- {total_std:.6e} "
          f"({total_rel_unc:.1f}% rel. unc.)")
    print(f"  (Sum of individual layers: {total_tritium_mean:.6e})")

    # =========================================================================
    # Extract neutron flux spectra in breeder layers
    # =========================================================================
    print("\n--- Neutron Flux Spectra in Breeder Layers ---")
    spectra_results = {}

    for bname in breeder_names:
        tally_name = f"flux_spectrum_{bname}"
        tally = sp.get_tally(name=tally_name)
        energy_centres, mean, std_dev = extract_spectrum(tally)

        # Report peak flux energy and integrated flux
        if mean.max() > 0:
            peak_idx = mean.argmax()
            peak_energy_mev = energy_centres[peak_idx] / 1.0e6
            print(f"  {bname:20s}: peak flux at {peak_energy_mev:.4f} MeV, "
                  f"integrated flux = {mean.sum():.4e}")
        else:
            print(f"  {bname:20s}: no flux recorded")

        spectra_results[bname] = {
            "energy_MeV": (energy_centres / 1.0e6).tolist(),
            "flux_mean": mean.tolist(),
            "flux_std_dev": std_dev.tolist(),
        }

    # =========================================================================
    # Extract Be-9(n,2n) reaction rates in beryllium zones
    # =========================================================================
    print("\n--- Be-9(n,2n) Neutron Multiplication Rates ---")
    print("  (n,2n) score: neutron multiplication in beryllium zones")
    print()

    be_names = ["be_zone1", "be_zone2", "be_zone3"]
    be_descriptions = [
        "Front zone (before first breeder)",
        "Central zone (between breeder sections, main multiplier)",
        "Rear zone (after second breeder)",
    ]
    be_results = {}

    for bname, bdesc in zip(be_names, be_descriptions):
        # (n,2n) rate
        n2n_tally = sp.get_tally(name=f"n2n_rate_{bname}")
        n2n_mean = n2n_tally.mean.flatten()[0]
        n2n_std = n2n_tally.std_dev.flatten()[0]

        # Total flux
        flux_tally = sp.get_tally(name=f"flux_{bname}")
        flux_mean = flux_tally.mean.flatten()[0]
        flux_std = flux_tally.std_dev.flatten()[0]

        print(f"  {bname:10s}: (n,2n) = {n2n_mean:.6e} +/- {n2n_std:.6e}, "
              f"flux = {flux_mean:.6e}  [{bdesc}]")

        be_results[bname] = {
            "description": bdesc,
            "n2n_rate_mean": float(n2n_mean),
            "n2n_rate_std_dev": float(n2n_std),
            "flux_mean": float(flux_mean),
            "flux_std_dev": float(flux_std),
        }

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "FNG HCPB Tritium Breeder Module Mock-up",
        "source": "SINBAD Database (NEA-1553)",
        "facility": "FNG, ENEA Frascati, Italy",
        "description": (
            "Mock-up of the HCPB breeding blanket concept with beryllium "
            "neutron multiplier and Li2CO3 breeder layers. Key benchmark "
            "for validating tritium breeding calculations."
        ),
        "runtime": runtime_info,
        "tritium_production": tritium_results,
        "total_tritium_production": {
            "mean": float(total_mean),
            "std_dev": float(total_std),
            "relative_uncertainty_pct": float(total_rel_unc),
        },
        "breeder_flux_spectra": spectra_results,
        "beryllium_zones": be_results,
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
