#!/usr/bin/env python3
"""
================================================================================
IPPE Iron Spherical Shell Transmission Benchmark - OpenMC Model
================================================================================

Source: SINBAD Database (Shielding Integral Benchmark Archive and Database)
        OECD Nuclear Energy Agency
        https://www.oecd-nea.org/science/wprs/shielding/sinbad/ippe-fe/

Experiment Description:
    Neutron transmission measurements through spherical iron shells of five
    different sizes were performed at the Institute of Physics and Power
    Engineering (IPPE), Obninsk, Russia, between 1989 and 1995. A 14.1 MeV
    D-T neutron source was placed at the center of each spherical shell, and
    the leakage neutron spectrum was measured using a fast scintillator
    detector (paraterphenyl crystal, 5 cm dia x 5 cm height) positioned at
    a flight path of 6.8 m. The detector covered an energy range from
    50 keV to 15 MeV.

Neutron Source:
    A Cockroft-Walton accelerator (KG-0.3 generator) produced 14.1 MeV
    neutrons via the D-T fusion reaction. The deuteron kinetic energy was
    280 keV maximum, incident on a titanium-tritide (Ti-T) target mounted
    on a 0.8 mm copper backing. The target diameter was 11 mm with a 5 mm
    beam spot. For simulation purposes, the source is modeled as an
    isotropic point source of 14.1 MeV neutrons at the origin.

Shell Configurations (5 spherical iron shells):
    Shell 1: Inner radius =  4.5 cm, thickness =  2.5 cm, outer radius =  7.0 cm
    Shell 2: Inner radius = 12.0 cm, thickness =  7.5 cm, outer radius = 19.5 cm
    Shell 3: Inner radius = 12.0 cm, thickness = 10.0 cm, outer radius = 22.0 cm
    Shell 4: Inner radius = 20.0 cm, thickness = 18.1 cm, outer radius = 38.1 cm
    Shell 5: Inner radius = 30.0 cm, thickness = 28.0 cm, outer radius = 58.0 cm

Material:
    Pure iron (Fe), density approximately 7.874 g/cm^3. The actual atom
    densities were calculated from measured sphere weights for each
    configuration. For this simulation, natural iron at standard density
    is used.

Reference:
    SINBAD database, "Neutron Transmission Through Iron Spherical Shells"
    IPPE, Obninsk, Russia
================================================================================
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Shell geometry data: (inner_radius_cm, wall_thickness_cm)
# Outer radius = inner_radius + wall_thickness
# =============================================================================
SHELL_DATA = {
    1: {"inner_radius": 4.5,  "thickness": 2.5,  "description": "Small thin shell"},
    2: {"inner_radius": 12.0, "thickness": 7.5,  "description": "Medium shell, 7.5 cm wall"},
    3: {"inner_radius": 12.0, "thickness": 10.0, "description": "Medium shell, 10 cm wall"},
    4: {"inner_radius": 20.0, "thickness": 18.1, "description": "Large thick shell"},
    5: {"inner_radius": 30.0, "thickness": 28.0, "description": "Largest shell, 28 cm wall"},
}


def build_model(shell_number: int, particles: int) -> openmc.Model:
    """
    Build the OpenMC model for the specified IPPE iron shell configuration.

    Parameters
    ----------
    shell_number : int
        Shell configuration number (1-5).
    particles : int
        Number of source particles to simulate.

    Returns
    -------
    openmc.Model
        Complete OpenMC model ready to run.
    """

    # =========================================================================
    # Validate shell selection
    # =========================================================================
    if shell_number not in SHELL_DATA:
        raise ValueError(f"Shell number must be 1-5, got {shell_number}")

    shell = SHELL_DATA[shell_number]
    r_inner = shell["inner_radius"]     # Inner radius in cm
    r_outer = r_inner + shell["thickness"]  # Outer radius in cm

    print(f"Building IPPE Iron Shell #{shell_number}: {shell['description']}")
    print(f"  Inner radius:   {r_inner:.1f} cm")
    print(f"  Wall thickness: {shell['thickness']:.1f} cm")
    print(f"  Outer radius:   {r_outer:.1f} cm")
    print(f"  Particles:      {particles:,}")

    # =========================================================================
    # Materials
    # =========================================================================
    # Pure iron - the shells were made of "pure iron" per the SINBAD spec.
    # Density is approximately 7.874 g/cm^3 (standard density of iron).
    # The actual experiment derived atom densities from measured sphere weights;
    # here we use the standard value as a representative approximation.
    # =========================================================================
    iron = openmc.Material(name="Iron (Fe)")
    iron.add_element("Fe", 1.0)            # 100% natural iron
    iron.set_density("g/cm3", 7.874)       # Standard iron density

    materials = openmc.Materials([iron])

    # =========================================================================
    # Geometry
    # =========================================================================
    # The geometry consists of three regions:
    #   1. Inner void (vacuum) - inside the spherical shell
    #   2. Iron shell wall     - the spherical shell itself
    #   3. Outer void (vacuum) - outside the shell, extends to boundary
    #
    # A vacuum boundary sphere is placed well beyond the outer shell radius
    # to terminate the geometry. Particles crossing this boundary are killed.
    # =========================================================================

    # --- Surfaces ---

    # Inner surface of the iron shell
    inner_sphere = openmc.Sphere(
        r=r_inner,
        name=f"Inner sphere (r={r_inner:.1f} cm)"
    )

    # Outer surface of the iron shell
    outer_sphere = openmc.Sphere(
        r=r_outer,
        name=f"Outer sphere (r={r_outer:.1f} cm)"
    )

    # Boundary sphere - vacuum boundary condition
    # Place it 50 cm beyond the outer shell to give sufficient room,
    # but this is just a geometric termination for the Monte Carlo transport.
    boundary_radius = r_outer + 50.0
    boundary_sphere = openmc.Sphere(
        r=boundary_radius,
        boundary_type="vacuum",
        name=f"Vacuum boundary (r={boundary_radius:.1f} cm)"
    )

    # --- Cells ---

    # Cell 1: Inner void cavity (vacuum inside the shell)
    # This is where the point source is located.
    inner_void = openmc.Cell(name="Inner void (source region)")
    inner_void.region = -inner_sphere     # Region inside the inner sphere
    # No material assigned -> vacuum

    # Cell 2: Iron spherical shell
    # The shell occupies the region between inner and outer sphere surfaces.
    shell_cell = openmc.Cell(name=f"Iron shell (t={shell['thickness']:.1f} cm)")
    shell_cell.region = +inner_sphere & -outer_sphere  # Between the two spheres
    shell_cell.fill = iron                              # Filled with iron

    # Cell 3: Outer void (vacuum outside the shell, up to boundary)
    # Neutrons that leak through the shell traverse this region.
    outer_void = openmc.Cell(name="Outer void (leakage region)")
    outer_void.region = +outer_sphere & -boundary_sphere  # Between shell and boundary
    # No material assigned -> vacuum

    # --- Universe and Geometry ---
    root_universe = openmc.Universe(
        cells=[inner_void, shell_cell, outer_void]
    )
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # Source Definition
    # =========================================================================
    # The D-T neutron source is modeled as an isotropic point source at the
    # origin emitting monoenergetic 14.1 MeV neutrons. In the actual
    # experiment, the source was a Ti-T target with a 5 mm beam spot, but
    # for this benchmark the point source approximation is standard practice
    # (the shell radii are much larger than the source size).
    # =========================================================================

    # Point source at the origin (center of the spherical shell)
    point_source = openmc.stats.Point((0.0, 0.0, 0.0))

    # Monoenergetic 14.1 MeV neutrons (D-T fusion energy)
    mono_energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV in eV

    # Isotropic angular distribution (default for IndependentSource)
    source = openmc.IndependentSource(
        space=point_source,
        energy=mono_energy,
        strength=1.0,
    )
    source.particle = "neutron"

    # =========================================================================
    # Settings
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = "fixed source"       # Fixed source (not eigenvalue)
    settings.particles = particles           # Number of histories per batch
    settings.batches = 10                    # 10 batches for uncertainty estimation
    settings.source = source                 # Assign the D-T point source
    settings.photon_transport = False        # Neutron-only transport

    # =========================================================================
    # Tallies
    # =========================================================================
    # The primary observable is the neutron leakage spectrum through the
    # outer surface of the iron shell. We use a surface current tally on
    # the outer sphere surface.
    #
    # Energy bins: logarithmic spacing from 0.05 MeV (50 keV) to 15 MeV.
    # This matches the experimental detector range (50 keV to 15 MeV).
    # We use 100 energy bins for good spectral resolution.
    # =========================================================================

    # Energy filter: 100 logarithmically-spaced bins from 50 keV to 15 MeV
    energy_bins = np.logspace(
        np.log10(0.05e6),   # 50 keV = 0.05 MeV, in eV
        np.log10(15.0e6),   # 15 MeV, in eV
        101                  # 101 edges -> 100 bins
    )
    energy_filter = openmc.EnergyFilter(energy_bins)

    # Surface filter: tally on the outer sphere of the iron shell
    surface_filter = openmc.SurfaceFilter(outer_sphere)

    # Current tally: counts neutrons crossing the outer surface
    # This gives us the neutron leakage spectrum per source particle.
    leakage_tally = openmc.Tally(name="Leakage spectrum")
    leakage_tally.filters = [surface_filter, energy_filter]
    leakage_tally.scores = ["current"]  # Surface current (particles/src)

    # Total leakage tally (without energy binning) for quick integral check
    total_leakage_tally = openmc.Tally(name="Total leakage")
    total_leakage_tally.filters = [surface_filter]
    total_leakage_tally.scores = ["current"]

    tallies = openmc.Tallies([leakage_tally, total_leakage_tally])

    # =========================================================================
    # Assemble and return the model
    # =========================================================================
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )

    return model


def main():
    """Parse arguments and run the IPPE iron shell benchmark."""

    # =========================================================================
    # Command-line argument parsing
    # =========================================================================
    parser = argparse.ArgumentParser(
        description=(
            "IPPE Iron Spherical Shell Transmission Benchmark.\n"
            "Simulates neutron leakage through iron shells with a "
            "14.1 MeV D-T source at the center."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Shell configurations:\n"
            "  1: r_in=4.5 cm,  t=2.5 cm   (thin, small)\n"
            "  2: r_in=12.0 cm, t=7.5 cm   (medium)\n"
            "  3: r_in=12.0 cm, t=10.0 cm  (medium, thicker)\n"
            "  4: r_in=20.0 cm, t=18.1 cm  (large, thick)\n"
            "  5: r_in=30.0 cm, t=28.0 cm  (largest, thickest)\n"
        ),
    )

    parser.add_argument(
        "--shell", type=int, default=1, choices=[1, 2, 3, 4, 5],
        help="Shell configuration number (1-5). Default: 1"
    )
    parser.add_argument(
        "--particles", type=int, default=1_000_000,
        help="Number of source particles. Default: 1,000,000"
    )

    args = parser.parse_args()

    # =========================================================================
    # Build and run the model
    # =========================================================================
    model = build_model(args.shell, args.particles)

    # Export XML input files (for inspection / debugging)
    model.export_to_xml()
    print("\nXML files exported. Running OpenMC...")

    # Run the simulation
    statepoint_path = model.run()
    print(f"\nSimulation complete. Statepoint: {statepoint_path}")


if __name__ == "__main__":
    main()
