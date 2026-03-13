#!/usr/bin/env python3
"""
ST-FNSF Benchmark -- Post-Processing & Analysis
=================================================

This script reads the OpenMC statepoint file produced by the ST-FNSF
benchmark simulation and extracts the key neutronics quantities:

  1. Tritium Breeding Ratio (TBR) -- the most critical figure of merit
  2. Nuclear heating in each component (center post, shield, FW, blanket, VV)
  3. Neutron flux and spectrum at the center column
  4. Fast neutron flux at the vacuum vessel (proxy for TF coil exposure)
  5. Inboard shield heating

For a spherical tokamak like ST-FNSF, the TBR is expected to be marginal
(~1.0) because there is no breeding blanket on the inboard side. Center
column heating and radiation damage are critical design constraints.

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
    energy_bins = energy_filter.bins
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
    # 1. Tritium Breeding Ratio (TBR)
    # =========================================================================
    print("\n--- Tritium Breeding Ratio (TBR) ---")
    tbr_tally = sp.get_tally(name="TBR")
    tbr_mean = tbr_tally.mean.flatten()[0]
    tbr_std = tbr_tally.std_dev.flatten()[0]
    tbr_rel_unc = (tbr_std / tbr_mean * 100.0) if tbr_mean > 0 else 0.0

    # The TBR from a 30-degree sector with reflective boundaries represents
    # the full-device TBR (per source neutron). For ST-FNSF, the expected
    # value is around 1.0 due to limited blanket coverage.
    print(f"  TBR = {tbr_mean:.4f} +/- {tbr_std:.4f} ({tbr_rel_unc:.1f}% rel. unc.)")

    if tbr_mean < 1.0:
        print("  WARNING: TBR < 1.0 -- tritium self-sufficiency NOT achieved")
        print("  This is expected for ST-FNSF due to absence of inboard blanket")
    else:
        print("  TBR >= 1.0 -- tritium self-sufficiency achieved")

    # =========================================================================
    # 2. Nuclear heating in all components
    # =========================================================================
    print("\n--- Nuclear Heating by Component ---")
    print("  (Units: eV per source neutron)")

    heating_tally = sp.get_tally(name="nuclear_heating")
    cell_filter = heating_tally.find_filter(openmc.CellFilter)
    cell_ids = cell_filter.bins

    # Component names corresponding to cell IDs in the filter
    component_names = [
        "Center post (Cu)",
        "Inboard shield (WC+water)",
        "First wall (F82H)",
        "DCLL blanket",
        "Outboard shield (WC+water)",
        "Vacuum vessel (SS316)",
    ]

    heating_results = {}
    total_heating = 0.0

    for i, name in enumerate(component_names):
        mean_val = heating_tally.mean.flatten()[i]
        std_val = heating_tally.std_dev.flatten()[i]
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0
        total_heating += mean_val

        print(f"  {name:35s}: {mean_val:.4e} +/- {std_val:.4e} "
              f"({rel_unc:.1f}%)")

        heating_results[name] = {
            "heating_eV_per_source": float(mean_val),
            "std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

    print(f"  {'Total':35s}: {total_heating:.4e} eV/source")

    # =========================================================================
    # 3. Center column neutron flux
    # =========================================================================
    print("\n--- Center Column Neutron Flux ---")
    cc_flux_tally = sp.get_tally(name="center_column_flux")
    cc_flux_mean = cc_flux_tally.mean.flatten()[0]
    cc_flux_std = cc_flux_tally.std_dev.flatten()[0]
    cc_flux_rel = (cc_flux_std / cc_flux_mean * 100.0) if cc_flux_mean > 0 else 0.0

    print(f"  Total flux = {cc_flux_mean:.4e} +/- {cc_flux_std:.4e} "
          f"({cc_flux_rel:.1f}% rel. unc.)")
    print(f"  [neutron-cm / source particle]")

    # =========================================================================
    # 4. Center column neutron spectrum
    # =========================================================================
    print("\n--- Center Column Neutron Spectrum ---")
    cc_spectrum_tally = sp.get_tally(name="center_column_spectrum")
    energies, spec_mean, spec_std = extract_spectrum(cc_spectrum_tally)

    # Find the peak energy
    peak_idx = np.argmax(spec_mean)
    peak_energy_MeV = energies[peak_idx] / 1.0e6
    print(f"  Peak flux at E = {peak_energy_MeV:.3f} MeV")
    print(f"  Peak value = {spec_mean[peak_idx]:.4e}")

    # Fraction of flux above 0.1 MeV (fast neutrons)
    fast_mask = energies > 0.1e6
    fast_fraction = spec_mean[fast_mask].sum() / spec_mean.sum() if spec_mean.sum() > 0 else 0
    print(f"  Fast fraction (E > 0.1 MeV) = {fast_fraction:.1%}")

    cc_spectrum_results = {
        "energy_MeV": (energies / 1.0e6).tolist(),
        "flux_mean": spec_mean.tolist(),
        "flux_std_dev": spec_std.tolist(),
        "peak_energy_MeV": float(peak_energy_MeV),
        "fast_fraction": float(fast_fraction),
    }

    # =========================================================================
    # 5. Fast flux at vacuum vessel
    # =========================================================================
    print("\n--- Fast Flux at Vacuum Vessel (TF Coil Proxy) ---")
    fast_flux_tally = sp.get_tally(name="fast_flux_at_vv")
    ff_mean = fast_flux_tally.mean.flatten()[0]
    ff_std = fast_flux_tally.std_dev.flatten()[0]
    ff_rel = (ff_std / ff_mean * 100.0) if ff_mean > 0 else 0.0

    print(f"  Fast flux (E > 0.1 MeV) = {ff_mean:.4e} +/- {ff_std:.4e} "
          f"({ff_rel:.1f}% rel. unc.)")

    # =========================================================================
    # 6. Inboard shield heating
    # =========================================================================
    print("\n--- Inboard Shield Heating ---")
    ib_heat_tally = sp.get_tally(name="inboard_shield_heating")
    ib_heat_mean = ib_heat_tally.mean.flatten()[0]
    ib_heat_std = ib_heat_tally.std_dev.flatten()[0]
    ib_heat_rel = (ib_heat_std / ib_heat_mean * 100.0) if ib_heat_mean > 0 else 0.0

    print(f"  Heating = {ib_heat_mean:.4e} +/- {ib_heat_std:.4e} "
          f"({ib_heat_rel:.1f}% rel. unc.) eV/source")

    # =========================================================================
    # 7. Blanket neutron spectrum
    # =========================================================================
    print("\n--- Blanket Neutron Spectrum ---")
    bl_spectrum_tally = sp.get_tally(name="blanket_spectrum")
    bl_energies, bl_mean, bl_std = extract_spectrum(bl_spectrum_tally)

    bl_peak_idx = np.argmax(bl_mean)
    bl_peak_MeV = bl_energies[bl_peak_idx] / 1.0e6
    print(f"  Peak flux at E = {bl_peak_MeV:.3f} MeV")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "ST-FNSF Spherical Tokamak (Menard et al. 2016)",
        "device_parameters": {
            "major_radius_cm": 170.0,
            "minor_radius_cm": 100.0,
            "aspect_ratio": 1.7,
            "sector_angle_deg": 30.0,
            "n_tf_coils": 12,
        },
        "runtime": runtime_info,
        "tbr": {
            "mean": float(tbr_mean),
            "std_dev": float(tbr_std),
            "relative_uncertainty_pct": float(tbr_rel_unc),
            "expected_range": "~1.0 (marginal, no inboard blanket)",
        },
        "heating_by_component": heating_results,
        "total_heating_eV_per_source": float(total_heating),
        "center_column": {
            "total_flux": float(cc_flux_mean),
            "total_flux_std": float(cc_flux_std),
            "spectrum": cc_spectrum_results,
        },
        "fast_flux_at_vv": {
            "mean": float(ff_mean),
            "std_dev": float(ff_std),
            "relative_uncertainty_pct": float(ff_rel),
        },
        "inboard_shield_heating": {
            "mean_eV_per_source": float(ib_heat_mean),
            "std_dev": float(ib_heat_std),
        },
        "blanket_spectrum": {
            "energy_MeV": (bl_energies / 1.0e6).tolist(),
            "flux_mean": bl_mean.tolist(),
            "flux_std_dev": bl_std.tolist(),
            "peak_energy_MeV": float(bl_peak_MeV),
        },
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

    # =========================================================================
    # Summary assessment
    # =========================================================================
    print("\n" + "=" * 70)
    print("ST-FNSF Neutronics Assessment Summary")
    print("=" * 70)
    print(f"  TBR:                  {tbr_mean:.4f} +/- {tbr_std:.4f}")
    if tbr_mean < 1.0:
        print(f"    -> Below unity (expected for ST with no inboard blanket)")
    print(f"  Center post flux:     {cc_flux_mean:.4e} n-cm/source")
    print(f"  Center post heating:  "
          f"{heating_results.get('Center post (Cu)', {}).get('heating_eV_per_source', 0):.4e} eV/source")
    print(f"  Inboard shield heat:  {ib_heat_mean:.4e} eV/source")
    print(f"  Fast flux at VV:      {ff_mean:.4e} n-cm/source")
    print(f"  Total heating:        {total_heating:.4e} eV/source")
    print("=" * 70)
    print("Done.")


if __name__ == "__main__":
    main()
