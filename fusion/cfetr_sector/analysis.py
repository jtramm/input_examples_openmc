#!/usr/bin/env python3
"""
CFETR 22.5-Degree Toroidal Sector -- Post-Processing & Analysis
==================================================================

This script reads the OpenMC statepoint file produced by the CFETR sector
benchmark and extracts:

  1. Tritium Breeding Ratio (TBR): total triton production per source neutron
  2. Nuclear heating distribution across all major components
  3. Be(n,2n) neutron multiplication rate in the blanket
  4. Fast neutron flux (E > 0.1 MeV) at the vacuum vessel
  5. Neutron energy spectrum in the blanket

Results are saved to results.json for comparison with published CFETR
neutronics analyses.

The key figure of merit is the TBR, which must exceed ~1.05 for tritium
self-sufficiency (accounting for processing losses). The CFETR design
target is TBR > 1.2 to provide adequate margin.

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
    print("\n" + "=" * 60)
    print("TRITIUM BREEDING RATIO (TBR)")
    print("=" * 60)
    print("  Score: (n,Xt) = total triton production per source neutron")
    print("  Breeding reactions:")
    print("    Li-6(n,t)He-4: dominant (sigma_th ~ 940 b, 1/v law)")
    print("    Li-7(n,n't)He-4: secondary (threshold 2.47 MeV)")
    print()

    tbr_tally = sp.get_tally(name="TBR")
    tbr_mean = tbr_tally.mean.flatten().sum()
    tbr_std = np.sqrt((tbr_tally.std_dev.flatten() ** 2).sum())
    tbr_rel_unc = (tbr_std / tbr_mean * 100.0) if tbr_mean > 0 else 0.0

    print(f"  TBR = {tbr_mean:.4f} +/- {tbr_std:.4f} "
          f"({tbr_rel_unc:.2f}% rel. unc.)")

    if tbr_mean > 1.2:
        print("  Status: EXCEEDS design target (TBR > 1.2)")
    elif tbr_mean > 1.05:
        print("  Status: Meets minimum requirement (TBR > 1.05) but below target")
    else:
        print("  Status: BELOW minimum requirement (TBR < 1.05)")

    # =========================================================================
    # 2. Nuclear heating in each component
    # =========================================================================
    print("\n" + "=" * 60)
    print("NUCLEAR HEATING (per source neutron)")
    print("=" * 60)
    print("  Score: heating = total energy deposition (eV/source)")
    print()

    component_names = [
        "First wall", "Blanket", "Back support",
        "VV inner", "VV fill", "VV outer",
        "Thermal shield", "TF coil",
    ]

    heating_tally = sp.get_tally(name="nuclear_heating")
    heating_mean = heating_tally.mean.flatten()
    heating_std = heating_tally.std_dev.flatten()

    heating_results = {}
    total_heating = 0.0

    for i, name in enumerate(component_names):
        mean_val = heating_mean[i]
        std_val = heating_std[i]
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        # Convert from eV to MeV for readability
        mean_mev = mean_val / 1.0e6
        std_mev = std_val / 1.0e6

        print(f"  {name:20s}: {mean_mev:.4e} +/- {std_mev:.4e} MeV/source "
              f"({rel_unc:.1f}%)")

        heating_results[name] = {
            "heating_eV_per_source": float(mean_val),
            "heating_std_dev_eV": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }
        total_heating += mean_val

    print(f"\n  {'TOTAL':20s}: {total_heating/1e6:.4e} MeV/source")

    # =========================================================================
    # 3. Be(n,2n) neutron multiplication rate
    # =========================================================================
    print("\n" + "=" * 60)
    print("NEUTRON MULTIPLICATION -- Be-9(n,2n)")
    print("=" * 60)
    print("  Score: (n,2n) = neutron multiplication in blanket")
    print("  Reaction: Be-9 + n -> 2 He-4 + 2n (threshold 1.85 MeV)")
    print()

    n2n_tally = sp.get_tally(name="Be_n2n_rate")
    n2n_mean = n2n_tally.mean.flatten()[0]
    n2n_std = n2n_tally.std_dev.flatten()[0]
    n2n_rel_unc = (n2n_std / n2n_mean * 100.0) if n2n_mean > 0 else 0.0

    print(f"  Be(n,2n) rate = {n2n_mean:.4e} +/- {n2n_std:.4e} "
          f"per source neutron ({n2n_rel_unc:.1f}%)")
    print(f"  Each (n,2n) event produces one EXTRA neutron for breeding")

    # =========================================================================
    # 4. Fast neutron flux at vacuum vessel
    # =========================================================================
    print("\n" + "=" * 60)
    print("FAST NEUTRON FLUX AT VACUUM VESSEL (E > 0.1 MeV)")
    print("=" * 60)

    fast_flux_tally = sp.get_tally(name="fast_flux_VV")
    fast_flux_mean = fast_flux_tally.mean.flatten()
    fast_flux_std = fast_flux_tally.std_dev.flatten()

    vv_names = ["VV inner wall", "VV fill", "VV outer wall"]
    fast_flux_results = {}

    for i, name in enumerate(vv_names):
        mean_val = fast_flux_mean[i]
        std_val = fast_flux_std[i]
        rel_unc = (std_val / mean_val * 100.0) if mean_val > 0 else 0.0

        print(f"  {name:20s}: {mean_val:.4e} +/- {std_val:.4e} "
              f"cm/source ({rel_unc:.1f}%)")

        fast_flux_results[name] = {
            "fast_flux_mean": float(mean_val),
            "fast_flux_std_dev": float(std_val),
            "relative_uncertainty_pct": float(rel_unc),
        }

    # =========================================================================
    # 5. Blanket neutron spectrum
    # =========================================================================
    print("\n" + "=" * 60)
    print("BLANKET NEUTRON SPECTRUM")
    print("=" * 60)

    spectrum_tally = sp.get_tally(name="blanket_spectrum")
    energy_filter = spectrum_tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])
    spec_mean = spectrum_tally.mean.flatten()
    spec_std = spectrum_tally.std_dev.flatten()

    if spec_mean.max() > 0:
        peak_idx = spec_mean.argmax()
        peak_energy_mev = energy_centres[peak_idx] / 1.0e6
        print(f"  Peak flux energy: {peak_energy_mev:.4f} MeV")
        print(f"  Integrated flux:  {spec_mean.sum():.4e} cm/source")
    else:
        print("  No flux recorded in blanket")

    # =========================================================================
    # Assemble and write results.json
    # =========================================================================
    results = {
        "benchmark": "CFETR 22.5-degree toroidal sector (WCCB blanket)",
        "reference": "Zhu et al., Plasma Sci. Technol. 18 (2016) 13",
        "doi": "10.1088/1009-0630/18/7/13",
        "description": (
            "Simplified 22.5-degree sector of the CFETR tokamak with "
            "homogenized Water-Cooled Ceramic Breeder (WCCB) blanket. "
            "Key metrics: TBR, nuclear heating, neutron multiplication."
        ),
        "runtime": runtime_info,
        "tritium_breeding_ratio": {
            "TBR_mean": float(tbr_mean),
            "TBR_std_dev": float(tbr_std),
            "relative_uncertainty_pct": float(tbr_rel_unc),
        },
        "nuclear_heating": heating_results,
        "total_heating_eV_per_source": float(total_heating),
        "Be_n2n_multiplication": {
            "rate_mean": float(n2n_mean),
            "rate_std_dev": float(n2n_std),
            "relative_uncertainty_pct": float(n2n_rel_unc),
        },
        "fast_flux_VV": fast_flux_results,
        "blanket_spectrum": {
            "energy_MeV": (energy_centres / 1.0e6).tolist(),
            "flux_mean": spec_mean.tolist(),
            "flux_std_dev": spec_std.tolist(),
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
