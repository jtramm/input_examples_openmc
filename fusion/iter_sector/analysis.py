#!/usr/bin/env python3
"""
Simplified ITER 40-Degree Sector -- Post-Processing & Analysis
=================================================================

This script reads the OpenMC statepoint file produced by the ITER sector
simulation and extracts:

  1. Nuclear heating in each shielding component (blanket, VV, TF coils)
  2. Fast neutron flux (E > 0.1 MeV) at the TF coil winding pack
  3. Neutron flux spectra at multiple radial depths (if available)
  4. First wall flux (proxy for neutron wall loading)
  5. Comparison with ITER reference values from the literature

Results are written to results.json for further analysis and benchmarking.

ITER Reference Values (from literature)
----------------------------------------
  - Peak neutron wall loading (outboard midplane): ~1.0 MW/m^2
  - Average neutron wall loading: ~0.56 MW/m^2
  - Total nuclear heating in blanket: ~85% of 400 MW neutron power
  - Fast neutron fluence at TF coil: target < 10^22 n/m^2 over lifetime
  - Nuclear heating in TF coils: < 10-20 kW total (to maintain 4.5 K)

Note that this simplified model uses circular toroidal shells rather than
the actual elongated, D-shaped cross-section, so absolute values will
differ from detailed ITER neutronics analyses. The relative attenuation
through the shielding layers should be qualitatively correct.

Usage:
    python analysis.py [--statepoint FILE]

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
# Command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="ITER sector benchmark -- post-processing and analysis"
)
parser.add_argument(
    "--statepoint",
    type=str,
    default=None,
    help="Path to statepoint file (default: latest statepoint.*.h5)",
)
cli_args = parser.parse_args()


def find_latest_statepoint():
    """Find the most recent statepoint file in the current directory.

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
    """Extract energy bin centres and flux values from a spectral tally.

    Parameters
    ----------
    tally : openmc.Tally
        A tally with an EnergyFilter.

    Returns
    -------
    energy_centres : numpy.ndarray
        Geometric mean of each energy bin [eV].
    mean : numpy.ndarray
        Mean flux values [neutrons-cm / source-particle].
    std_dev : numpy.ndarray
        Standard deviation of flux values.
    """
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2)
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])
    mean = tally.mean.flatten()
    std_dev = tally.std_dev.flatten()
    return energy_centres, mean, std_dev


def safe_rel_unc(mean_val, std_val):
    """Compute relative uncertainty in percent, handling zero mean."""
    if mean_val > 0:
        return std_val / mean_val * 100.0
    return 0.0


def main():
    # =========================================================================
    # Load the statepoint file
    # =========================================================================
    sp_path = cli_args.statepoint if cli_args.statepoint else find_latest_statepoint()
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

    results = {
        "benchmark": "Simplified ITER 40-degree sector",
        "description": (
            "Simplified CSG model of an ITER tokamak 40-degree sector with "
            "concentric toroidal shells: Be armor, Cu heat sink, SS316LN "
            "first wall, shielding blanket (SS316LN+H2O), vacuum vessel, "
            "thermal shield, and TF coils. D-T neutron source at 14.1 MeV."
        ),
        "runtime": runtime_info,
    }

    # =========================================================================
    # Try to extract heating tallies (full tally set)
    # =========================================================================
    print("\n--- Nuclear Heating by Component ---")
    heating_names = [
        ("heating_be_armor",         "Be armor (first wall)"),
        ("heating_cu_heatsink",      "Cu heat sink"),
        ("heating_fw_structure",     "First wall SS316LN"),
        ("heating_blanket",          "Shielding blanket"),
        ("heating_vv_inner",         "VV inner shell"),
        ("heating_vv_fill",          "VV fill"),
        ("heating_vv_outer",         "VV outer shell"),
        ("heating_thermal_shield",   "Thermal shield"),
        ("heating_tf_case_inner",    "TF inner case"),
        ("heating_tf_winding",       "TF winding pack"),
        ("heating_tf_case_outer",    "TF outer case"),
    ]

    heating_results = {}
    total_heating = 0.0
    has_component_heating = False

    for tally_name, component_name in heating_names:
        try:
            tally = sp.get_tally(name=tally_name)
            mean_val = float(tally.mean.flatten()[0])
            std_val = float(tally.std_dev.flatten()[0])
            rel_unc = safe_rel_unc(mean_val, std_val)

            print(f"  {component_name:30s}: {mean_val:12.4e} +/- {std_val:.4e} "
                  f"eV/src  ({rel_unc:.1f}%)")

            heating_results[tally_name] = {
                "component": component_name,
                "heating_eV_per_source": mean_val,
                "std_dev_eV_per_source": std_val,
                "relative_uncertainty_pct": rel_unc,
            }
            total_heating += mean_val
            has_component_heating = True
        except Exception:
            pass

    # Try the combined heating tally (small tally set)
    if not has_component_heating:
        try:
            tally = sp.get_tally(name="total_heating")
            mean_val = float(tally.mean.flatten().sum())
            std_val = float(np.sqrt(np.sum(tally.std_dev.flatten()**2)))
            rel_unc = safe_rel_unc(mean_val, std_val)
            print(f"  {'Total (all components)':30s}: {mean_val:12.4e} +/- "
                  f"{std_val:.4e} eV/src  ({rel_unc:.1f}%)")
            heating_results["total_heating"] = {
                "component": "All components (combined)",
                "heating_eV_per_source": mean_val,
                "std_dev_eV_per_source": std_val,
                "relative_uncertainty_pct": rel_unc,
            }
            total_heating = mean_val
        except Exception:
            print("  No heating tallies found.")

    if has_component_heating:
        print(f"\n  {'Total heating':30s}: {total_heating:12.4e} eV/src")

    results["nuclear_heating"] = heating_results

    # =========================================================================
    # Fast neutron flux at TF winding pack
    # =========================================================================
    print("\n--- Fast Neutron Flux at TF Winding Pack (E > 0.1 MeV) ---")
    try:
        tally_tf = sp.get_tally(name="tf_coil_fast_flux")
        mean_val = float(tally_tf.mean.flatten()[0])
        std_val = float(tally_tf.std_dev.flatten()[0])
        rel_unc = safe_rel_unc(mean_val, std_val)

        print(f"  Fast flux: {mean_val:.4e} +/- {std_val:.4e} "
              f"n-cm/src  ({rel_unc:.1f}%)")

        results["tf_coil_fast_flux"] = {
            "flux_n_cm_per_source": mean_val,
            "std_dev": std_val,
            "relative_uncertainty_pct": rel_unc,
            "note": (
                "ITER design limit: cumulative fast fluence < 10^22 n/m^2 "
                "at TF coil winding pack over device lifetime"
            ),
        }

        if rel_unc > 50.0:
            print("  WARNING: High statistical uncertainty. Use weight windows "
                  "(--generate-ww) for reliable TF coil flux estimates.")
    except Exception:
        print("  TF coil fast flux tally not found.")

    # =========================================================================
    # Neutron flux spectra at multiple depths
    # =========================================================================
    spectrum_names = [
        ("spectrum_first_wall", "First wall"),
        ("spectrum_blanket",    "Blanket"),
        ("spectrum_vv",         "Vacuum vessel"),
        ("spectrum_tf_coil",    "TF coil"),
    ]

    spectra_results = {}
    has_spectra = False

    print("\n--- Neutron Flux Spectra ---")
    for tally_name, location in spectrum_names:
        try:
            tally = sp.get_tally(name=tally_name)
            energy_centres, mean, std_dev = extract_spectrum(tally)

            total_flux = mean.sum()
            peak_idx = mean.argmax()
            peak_energy = energy_centres[peak_idx]

            print(f"  {location:20s}: total flux = {total_flux:.4e} n-cm/src, "
                  f"peak at {peak_energy/1e6:.3f} MeV")

            spectra_results[tally_name] = {
                "location": location,
                "energy_MeV": (energy_centres / 1.0e6).tolist(),
                "flux_mean": mean.tolist(),
                "flux_std_dev": std_dev.tolist(),
                "total_flux": float(total_flux),
                "peak_energy_MeV": float(peak_energy / 1.0e6),
            }
            has_spectra = True
        except Exception:
            pass

    if not has_spectra:
        print("  No spectral tallies found (run without --small-tallies).")

    results["neutron_spectra"] = spectra_results

    # =========================================================================
    # First wall flux (proxy for neutron wall loading)
    # =========================================================================
    print("\n--- First Wall Neutron Flux ---")
    try:
        tally_fw = sp.get_tally(name="first_wall_flux")
        df = tally_fw.get_pandas_dataframe()

        for score in ["flux", "current"]:
            score_data = df[df["score"] == score]
            if not score_data.empty:
                mean_val = float(score_data["mean"].values[0])
                std_val = float(score_data["std. dev."].values[0])
                rel_unc = safe_rel_unc(mean_val, std_val)
                print(f"  {score:10s}: {mean_val:.4e} +/- {std_val:.4e}  "
                      f"({rel_unc:.1f}%)")

        results["first_wall_flux"] = {
            "note": "Proxy for neutron wall loading (NWL). "
                    "ITER reference: peak NWL ~1.0 MW/m^2 at outboard midplane, "
                    "average ~0.56 MW/m^2.",
        }
    except Exception:
        print("  First wall flux tally not found.")

    # =========================================================================
    # Comparison with ITER reference values
    # =========================================================================
    print("\n--- ITER Reference Comparison ---")
    print("  Note: This simplified model uses circular toroidal cross-sections")
    print("  rather than the actual D-shaped elongated plasma. Absolute values")
    print("  will differ from detailed ITER neutronics analyses.")
    print()
    print("  ITER design targets:")
    print("    Peak neutron wall loading:        ~1.0 MW/m^2")
    print("    Total nuclear heating (blanket):   ~340 MW (85% of 400 MW)")
    print("    TF coil nuclear heating:           < 10-20 kW total")
    print("    TF coil fast fluence (lifetime):   < 10^22 n/m^2")
    print("    Shield attenuation factor:         ~10^10 (plasma to TF coil)")

    results["reference_values"] = {
        "peak_nwl_MW_per_m2": 1.0,
        "average_nwl_MW_per_m2": 0.56,
        "blanket_heating_fraction": 0.85,
        "tf_coil_heating_limit_kW": 20.0,
        "tf_coil_fast_fluence_limit_n_per_m2": 1.0e22,
        "source": "ITER Organization, ITR-18-003 (2018)",
    }

    # =========================================================================
    # Write results to JSON
    # =========================================================================
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
