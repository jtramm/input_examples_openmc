#!/usr/bin/env python3
"""
ARC Reactor Fusion Neutronics Benchmark -- Post-Processing & Analysis
=======================================================================

This script reads the OpenMC statepoint file produced by the ARC reactor
simulation and extracts key neutronics figures of merit:

  1. Tritium Breeding Ratio (TBR) -- the most critical metric for a
     self-sustaining D-T fusion reactor. TBR > 1 is required.
  2. Nuclear heating in the FLiBe blanket -- determines thermal power
     extraction and temperature profiles.
  3. Neutron flux at the TF coil -- determines radiation damage rate
     and lifetime of the HTS magnets.
  4. Energy-dependent neutron spectrum at the coil -- characterises
     the radiation environment for magnet qualification.

Results are compared against published reference values from:
  - Sorbom et al. (2015): ARC design paper
  - Kuang et al. (2018): ARC heat exhaust study
  - Bae, Peterson, Shimwell (2022): TBR benchmark with 3 codes

Usage:
    python analysis.py [--statepoint PATH]

Input:
    statepoint.*.h5  (latest statepoint file from OpenMC run)

Output:
    results.json     (structured results for comparison / archiving)
"""

import argparse
import glob
import json
import sys

import numpy as np
import openmc


# =============================================================================
# Reference values from literature
# =============================================================================
# TBR reference values for the ARC reactor with natural lithium FLiBe blanket.
# These come from Bae, Peterson, and Shimwell (2022), who benchmarked three
# Monte Carlo codes (MCNP, Serpent, OpenMC) on the ARC geometry.
TBR_REFERENCE = {
    "natural_li_low": 1.05,
    "natural_li_high": 1.10,
    "enriched_li6_90pct": 1.30,
    "code_agreement_pct": 0.6,  # three codes agreed within 0.6%
}

# Fusion power for normalisation (525 MW total D-T fusion power)
FUSION_POWER_MW = 525.0
# Energy per D-T fusion reaction: 17.6 MeV (14.1 MeV neutron + 3.5 MeV alpha)
ENERGY_PER_FUSION_EV = 17.6e6  # eV
# Source rate: fusions per second = P_fusion / E_per_fusion
# (Used to convert per-source-particle tallies to absolute quantities)
EV_TO_JOULE = 1.602176634e-19
FUSIONS_PER_SECOND = (FUSION_POWER_MW * 1.0e6) / (ENERGY_PER_FUSION_EV * EV_TO_JOULE)


def find_latest_statepoint():
    """Find the most recent statepoint file in the current directory.

    OpenMC writes statepoint files named statepoint.<batch>.h5. This
    function returns the one with the highest batch number.

    Returns
    -------
    str
        Path to the latest statepoint file.
    """
    statepoint_files = sorted(glob.glob("statepoint.*.h5"))
    if not statepoint_files:
        print("ERROR: No statepoint files found in the current directory.")
        print("       Run 'python model.py' first to generate simulation output.")
        sys.exit(1)
    return statepoint_files[-1]


def extract_tbr(sp):
    """Extract the Tritium Breeding Ratio from the statepoint.

    The TBR is the mean value of the (n,Xt) tally in the FLiBe blanket.
    Since the tally is per source neutron and one source neutron
    corresponds to one fusion reaction, the TBR is directly the tally
    mean value (no normalisation needed).

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file.

    Returns
    -------
    dict
        Dictionary with TBR mean, standard deviation, and comparison
        to reference values.
    """
    tally = sp.get_tally(name="TBR")
    mean = float(tally.mean.flatten()[0])
    std_dev = float(tally.std_dev.flatten()[0])
    rel_unc_pct = (std_dev / mean * 100.0) if mean > 0 else 0.0

    # Compare to reference range
    in_range = TBR_REFERENCE["natural_li_low"] <= mean <= TBR_REFERENCE["natural_li_high"]

    result = {
        "tbr_mean": mean,
        "tbr_std_dev": std_dev,
        "tbr_relative_uncertainty_pct": rel_unc_pct,
        "reference_range": [TBR_REFERENCE["natural_li_low"],
                            TBR_REFERENCE["natural_li_high"]],
        "within_reference_range": in_range,
    }

    print("\n--- Tritium Breeding Ratio (TBR) ---")
    print(f"  TBR = {mean:.4f} +/- {std_dev:.4f} ({rel_unc_pct:.2f}% rel. unc.)")
    print(f"  Reference range (natural Li): "
          f"{TBR_REFERENCE['natural_li_low']:.2f} - "
          f"{TBR_REFERENCE['natural_li_high']:.2f}")
    if in_range:
        print("  STATUS: Within reference range")
    else:
        print("  STATUS: Outside reference range (may indicate geometry or "
              "material differences)")

    return result


def extract_blanket_heating(sp):
    """Extract nuclear heating in the FLiBe blanket.

    The 'heating' score gives energy deposited in eV per source particle.
    We convert to absolute power using the fusion reaction rate.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file.

    Returns
    -------
    dict
        Dictionary with heating values in eV/source and MW.
    """
    tally = sp.get_tally(name="Blanket Heating")
    mean_ev = float(tally.mean.flatten()[0])
    std_dev_ev = float(tally.std_dev.flatten()[0])

    # Convert from eV/source-particle to watts, then to MW
    # Power [W] = (eV/source) * (sources/second) * (J/eV)
    # Note: This is for the 20-degree sector; multiply by 18 for full torus.
    sector_power_w = mean_ev * FUSIONS_PER_SECOND * EV_TO_JOULE
    full_torus_power_mw = sector_power_w * 18.0 / 1.0e6

    rel_unc_pct = (std_dev_ev / mean_ev * 100.0) if mean_ev > 0 else 0.0

    result = {
        "heating_ev_per_source": mean_ev,
        "heating_std_dev_ev": std_dev_ev,
        "heating_relative_uncertainty_pct": rel_unc_pct,
        "sector_power_mw": sector_power_w / 1.0e6,
        "full_torus_blanket_power_mw": full_torus_power_mw,
        "fusion_power_mw": FUSION_POWER_MW,
        "energy_multiplication": full_torus_power_mw / FUSION_POWER_MW if FUSION_POWER_MW > 0 else 0.0,
    }

    print("\n--- Nuclear Heating in FLiBe Blanket ---")
    print(f"  Heating per source neutron: {mean_ev:.4e} eV "
          f"({rel_unc_pct:.2f}% rel. unc.)")
    print(f"  Sector blanket power:       {sector_power_w / 1.0e6:.2f} MW "
          f"(20-deg sector)")
    print(f"  Full torus blanket power:   {full_torus_power_mw:.1f} MW")
    print(f"  Energy multiplication:      "
          f"{full_torus_power_mw / FUSION_POWER_MW:.3f}")

    return result


def extract_coil_flux(sp):
    """Extract neutron flux at the TF coil region.

    This is a critical design quantity: the neutron flux determines
    the radiation damage rate in the HTS magnets. Typical design
    limits are ~1e18 n/cm^2 total fluence over the magnet lifetime.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file.

    Returns
    -------
    dict
        Dictionary with flux values per source particle and absolute.
    """
    tally = sp.get_tally(name="Coil Neutron Flux")
    mean = float(tally.mean.flatten()[0])
    std_dev = float(tally.std_dev.flatten()[0])
    rel_unc_pct = (std_dev / mean * 100.0) if mean > 0 else 0.0

    # Convert to absolute flux: flux [n/cm^2/s] = (tally/source) * (source/s)
    # Note: tally is track-length in cm per source, need to divide by volume
    # for volumetric flux. But track-length/volume = scalar flux.
    # For comparison purposes, we report the raw tally value.
    abs_flux = mean * FUSIONS_PER_SECOND

    result = {
        "flux_per_source": mean,
        "flux_std_dev": std_dev,
        "flux_relative_uncertainty_pct": rel_unc_pct,
        "absolute_flux_n_cm2_s": abs_flux,
    }

    print("\n--- Neutron Flux at TF Coil ---")
    print(f"  Flux (per source): {mean:.4e} +/- {std_dev:.4e} cm "
          f"({rel_unc_pct:.2f}% rel. unc.)")
    print(f"  Absolute flux:     {abs_flux:.4e} n*cm/s (track-length)")
    print("  Note: Divide by coil cell volume [cm^3] for scalar flux [n/cm^2/s]")

    return result


def extract_spectrum(sp):
    """Extract the energy-dependent neutron spectrum at the TF coil.

    This shows the spectral character of the radiation field reaching
    the magnets: how much is thermalised by the blanket vs. how much
    penetrates at high energy.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file.

    Returns
    -------
    dict
        Dictionary with energy bin centres and flux values.
    """
    tally = sp.get_tally(name="Shield Leakage Spectrum")
    energy_filter = tally.find_filter(openmc.EnergyFilter)
    energy_bins = energy_filter.bins  # shape (N, 2): [low, high] per bin

    # Geometric mean of bin edges gives representative energy
    energy_centres = np.sqrt(energy_bins[:, 0] * energy_bins[:, 1])

    mean = tally.mean.flatten()
    std_dev = tally.std_dev.flatten()

    # Find the peak energy bin
    peak_idx = np.argmax(mean)
    peak_energy_mev = energy_centres[peak_idx] / 1.0e6

    result = {
        "energy_centres_eV": energy_centres.tolist(),
        "energy_centres_MeV": (energy_centres / 1.0e6).tolist(),
        "flux_mean": mean.tolist(),
        "flux_std_dev": std_dev.tolist(),
        "peak_energy_MeV": float(peak_energy_mev),
        "peak_flux": float(mean[peak_idx]),
        "n_energy_bins": len(energy_centres),
    }

    print("\n--- Neutron Spectrum at TF Coil ---")
    print(f"  Number of energy bins: {len(energy_centres)}")
    print(f"  Peak flux energy:      {peak_energy_mev:.4f} MeV")
    print(f"  Peak flux value:       {mean[peak_idx]:.4e} per source")

    return result


def extract_shield_flux(sp):
    """Extract neutron flux in the shield region.

    Useful for evaluating the shielding effectiveness of the borated
    steel layer.

    Parameters
    ----------
    sp : openmc.StatePoint
        Loaded statepoint file.

    Returns
    -------
    dict
        Dictionary with shield flux values.
    """
    tally = sp.get_tally(name="Shield Neutron Flux")
    mean = float(tally.mean.flatten()[0])
    std_dev = float(tally.std_dev.flatten()[0])
    rel_unc_pct = (std_dev / mean * 100.0) if mean > 0 else 0.0

    result = {
        "flux_per_source": mean,
        "flux_std_dev": std_dev,
        "flux_relative_uncertainty_pct": rel_unc_pct,
    }

    print("\n--- Neutron Flux in Shield Region ---")
    print(f"  Flux (per source): {mean:.4e} +/- {std_dev:.4e} cm "
          f"({rel_unc_pct:.2f}% rel. unc.)")

    return result


def main():
    """Main analysis routine."""
    # =========================================================================
    # Parse arguments
    # =========================================================================
    ap = argparse.ArgumentParser(
        description="ARC Reactor Benchmark -- Post-Processing"
    )
    ap.add_argument(
        "--statepoint",
        type=str,
        default=None,
        help="Path to statepoint file (default: latest statepoint.*.h5)",
    )
    cli_args = ap.parse_args()

    # =========================================================================
    # Load the statepoint file
    # =========================================================================
    sp_path = cli_args.statepoint or find_latest_statepoint()
    print(f"Loading statepoint: {sp_path}")
    sp = openmc.StatePoint(sp_path)

    runtime_info = {
        "statepoint_file": sp_path,
        "n_batches": sp.n_batches,
        "particles_per_batch": sp.n_particles,
        "total_particles": sp.n_batches * sp.n_particles,
    }
    print(f"  Batches:           {sp.n_batches}")
    print(f"  Particles/batch:   {sp.n_particles:,}")
    print(f"  Total particles:   {sp.n_batches * sp.n_particles:,}")

    # =========================================================================
    # Extract all tallies
    # =========================================================================
    tbr_result = extract_tbr(sp)
    heating_result = extract_blanket_heating(sp)
    coil_flux_result = extract_coil_flux(sp)
    spectrum_result = extract_spectrum(sp)
    shield_flux_result = extract_shield_flux(sp)

    # =========================================================================
    # Summary comparison with reference
    # =========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY: ARC Reactor Neutronics Results")
    print("=" * 70)
    print(f"  TBR:                    {tbr_result['tbr_mean']:.4f} "
          f"(ref: {TBR_REFERENCE['natural_li_low']:.2f}-"
          f"{TBR_REFERENCE['natural_li_high']:.2f})")
    print(f"  Blanket power (full):   "
          f"{heating_result['full_torus_blanket_power_mw']:.1f} MW "
          f"(fusion: {FUSION_POWER_MW:.0f} MW)")
    print(f"  Energy multiplication:  "
          f"{heating_result['energy_multiplication']:.3f}")
    print(f"  Coil flux (per source): {coil_flux_result['flux_per_source']:.4e}")
    print("=" * 70)

    # =========================================================================
    # Write results to JSON
    # =========================================================================
    class NumpyEncoder(json.JSONEncoder):
        """Custom JSON encoder that handles numpy types."""
        def default(self, obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)

    results = {
        "benchmark": "ARC Reactor Fusion Neutronics",
        "design": "MIT ARC compact high-field tokamak",
        "fusion_power_mw": FUSION_POWER_MW,
        "runtime": runtime_info,
        "tbr": tbr_result,
        "blanket_heating": heating_result,
        "coil_neutron_flux": coil_flux_result,
        "shield_neutron_flux": shield_flux_result,
        "coil_spectrum": spectrum_result,
        "references": {
            "sorbom_2015": "Sorbom et al., Fusion Eng. Des. 100, 378-405 (2015)",
            "kuang_2018": "Kuang et al., Fusion Eng. Des. 137, 221-242 (2018)",
            "bae_2022": "Bae, Peterson, Shimwell, Fusion Sci. Technol. (2022)",
        },
    }

    output_file = "results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)

    print(f"\nResults written to: {output_file}")
    print("Done.")


if __name__ == "__main__":
    main()
