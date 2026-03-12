#!/usr/bin/env python3
"""
OKTAVIAN Iron Sphere Benchmark -- OpenMC Input Model
=====================================================

Background
----------
The OKTAVIAN (Osaka-university Toolkit for Verification and Analysis of
Neutronics) facility is located at Osaka University, Japan.  It was designed
specifically for benchmarking nuclear data and neutron transport codes used in
fusion reactor design.

The experiment is conceptually simple:

  1. A deuterium-tritium (D-T) neutron source is placed at the exact center
     of a hollow sphere made of the test material (in this case, iron /
     carbon steel).

  2. The D-T reaction produces nearly monoenergetic 14.1 MeV neutrons that
     are emitted isotropically.

  3. As these neutrons travel outward through the sphere material they
     undergo elastic scattering, inelastic scattering, and absorption
     reactions.  Each interaction removes energy from the neutron or removes
     the neutron entirely.  The result is that neutrons "leak" out of the
     sphere surface with a broad energy spectrum -- from the original
     14.1 MeV all the way down to thermal energies.

  4. A time-of-flight (TOF) detector located approximately 9.5 m from the
     sphere measures the energy spectrum of the leaking neutrons.

The measured leakage spectrum is a sensitive test of the accuracy of iron
cross sections, particularly:
  - Elastic scattering angular distributions
  - Inelastic scattering to discrete levels (e.g., the first excited state
    of Fe-56 at 0.847 MeV)
  - The (n,2n) threshold reaction
  - Resonance structure in the keV range

This benchmark is part of the SINBAD (Shielding Integral Benchmark Archive
Database) collection maintained by the OECD Nuclear Energy Agency.

This script builds the complete OpenMC model for the iron sphere case:
  - Sphere radius: 50.32 cm
  - Material: carbon steel (approximated as natural iron, 7.874 g/cc)
  - Source: 14.1 MeV isotropic point source at the origin
  - Tallies: surface current leakage spectrum and cell flux spectrum

Reference
---------
SINBAD: OKTAVIAN Iron Sphere
https://www.oecd-nea.org/science/wprs/shielding/sinbad/oktav_fe/okfe-abs.htm
"""

import argparse
import numpy as np

import openmc


def build_model(particles: int = 1_000_000) -> openmc.Model:
    """
    Construct and return a complete openmc.Model for the OKTAVIAN iron sphere.

    Parameters
    ----------
    particles : int
        Number of source particles to simulate (default 1 000 000).

    Returns
    -------
    openmc.Model
        Fully configured model ready to export or run.
    """

    # =========================================================================
    # Materials
    # =========================================================================
    # The actual OKTAVIAN sphere was made of carbon steel (JIS SS41).
    # The composition by weight is approximately:
    #   Fe  98.69 %
    #   Mn   0.50 %
    #   Si   0.25 %
    #   C    0.20 %
    # with trace amounts of P and S (ignored here).
    #
    # The density of carbon steel is about 7.874 g/cm^3.
    # We include the four major constituents and use natural isotopic
    # abundances for each element.

    iron_steel = openmc.Material(name="Carbon Steel (JIS SS41)")
    iron_steel.set_density("g/cm3", 7.874)

    # Weight fractions of the major constituents
    iron_steel.add_element("Fe", 0.9869, percent_type="wo")  # 98.69 wt%
    iron_steel.add_element("Mn", 0.0050, percent_type="wo")  #  0.50 wt%
    iron_steel.add_element("Si", 0.0025, percent_type="wo")  #  0.25 wt%
    iron_steel.add_element("C",  0.0020, percent_type="wo")  #  0.20 wt%
    # NOTE: the remaining 0.36 wt% (P, S, etc.) is lumped into Fe since
    # those trace elements have negligible neutronic impact at this level.

    # Collect all materials into an openmc.Materials object
    materials = openmc.Materials([iron_steel])

    # =========================================================================
    # Geometry
    # =========================================================================
    # The geometry consists of three regions:
    #
    #   1. The iron sphere itself (radius 50.32 cm, centered at origin).
    #   2. A void region outside the sphere (vacuum -- no material).
    #   3. An outer bounding sphere that serves as the "graveyard" -- any
    #      particle that reaches this surface is terminated.  We place it
    #      at 200 cm radius, well beyond the iron sphere.
    #
    # In OpenMC, surfaces divide space into half-spaces.  For a sphere
    # centered at the origin:
    #   - The "negative" side (operator -) is the interior (r < R).
    #   - The "positive" side (operator +) is the exterior (r > R).

    # --- Surfaces ---

    # Inner sphere surface (the physical iron sphere boundary)
    sphere_surface = openmc.Sphere(
        r=50.32,
        name="Iron sphere outer surface (R = 50.32 cm)",
    )

    # Outer bounding sphere with vacuum boundary condition.
    # Any particle crossing this surface is killed (removed from tracking).
    outer_surface = openmc.Sphere(
        r=200.0,
        boundary_type="vacuum",
        name="Outer vacuum boundary (R = 200 cm)",
    )

    # --- Cells ---

    # Cell 1: The iron sphere.
    # Region is the interior of the sphere surface (negative half-space).
    iron_cell = openmc.Cell(name="Iron sphere")
    iron_cell.region = -sphere_surface          # inside the sphere
    iron_cell.fill = iron_steel                  # filled with carbon steel

    # Cell 2: Void region between the iron sphere and the outer boundary.
    # Region is outside the iron sphere AND inside the outer boundary.
    # No material fill means vacuum.
    void_cell = openmc.Cell(name="Void (outside sphere)")
    void_cell.region = +sphere_surface & -outer_surface  # between the two spheres

    # Build the universe and geometry
    root_universe = openmc.Universe(cells=[iron_cell, void_cell])
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # Source Definition
    # =========================================================================
    # The D-T neutron source produces 14.1 MeV neutrons isotropically from
    # a point at the center of the sphere (the origin).
    #
    # In a real D-T source the neutron energy has a slight spread around
    # 14.06-14.08 MeV depending on kinematics, but the standard benchmark
    # specification uses a monoenergetic 14.1 MeV source.

    # Monoenergetic energy distribution at 14.1 MeV
    energy_dist = openmc.stats.Discrete([14.1e6], [1.0])  # energy in eV

    # Point spatial distribution at the origin
    spatial_dist = openmc.stats.Point((0.0, 0.0, 0.0))

    # Isotropic angular distribution (default, but stated explicitly for clarity)
    angle_dist = openmc.stats.Isotropic()

    # Assemble the independent source
    source = openmc.IndependentSource(
        space=spatial_dist,
        angle=angle_dist,
        energy=energy_dist,
        strength=1.0,          # relative source strength (only one source)
    )

    # =========================================================================
    # Settings
    # =========================================================================
    # This is a FIXED SOURCE problem, not an eigenvalue (k-effective) problem.
    # In fixed-source mode OpenMC tracks a fixed number of source particles
    # per batch (there is no fission source iteration).

    settings = openmc.Settings()
    settings.run_mode = "fixed source"     # NOT 'eigenvalue'
    settings.batches = 10                  # number of batches
    settings.particles = particles // 10   # particles per batch
    # Total particles = batches * particles_per_batch = particles
    # Splitting into batches allows us to estimate statistical uncertainty
    # via batch-to-batch variation.

    settings.source = source               # attach the D-T source

    # Photon transport is off by default; we keep it off since the benchmark
    # focuses on neutron leakage spectra.
    settings.photon_transport = False

    # =========================================================================
    # Tallies
    # =========================================================================
    # We define two tallies:
    #
    #   1. Surface current tally on the iron sphere surface.
    #      This counts every neutron that crosses the sphere boundary
    #      (i.e., leaks out of the sphere) binned by energy.  This is the
    #      primary observable of the OKTAVIAN experiment.
    #
    #   2. Cell flux tally in the iron sphere.
    #      This gives the scalar neutron flux (track-length estimator) inside
    #      the iron, binned by energy.  It is useful for understanding the
    #      in-sphere spectrum and for comparison with other codes.

    # --- Energy bin structure ---
    # Logarithmic bins from 0.01 MeV (10 keV) to 15 MeV.
    # We use 100 bins for reasonable energy resolution.
    # Note: OpenMC energies are in eV, so 0.01 MeV = 1e4 eV, 15 MeV = 1.5e7 eV.
    energy_bins = np.logspace(
        np.log10(1.0e4),    # 10 keV  = 1e4 eV  (lower bound)
        np.log10(1.5e7),    # 15 MeV  = 1.5e7 eV (upper bound)
        101,                 # 101 edges -> 100 bins
    )

    energy_filter = openmc.EnergyFilter(energy_bins)

    # --- Tally 1: Neutron leakage spectrum (surface current) ---
    # A surface current tally counts particles crossing a specified surface.
    # By adding an energy filter we get the leakage spectrum.
    surface_filter = openmc.SurfaceFilter([sphere_surface])

    leakage_tally = openmc.Tally(name="Neutron leakage spectrum")
    leakage_tally.filters = [surface_filter, energy_filter]
    leakage_tally.scores = ["current"]
    # "current" counts particles crossing the surface in either direction.
    # For this geometry (source inside, vacuum outside) essentially all
    # crossings are outward, so current ~ leakage.

    # --- Tally 2: Cell flux spectrum in the iron sphere ---
    cell_filter = openmc.CellFilter([iron_cell])

    flux_tally = openmc.Tally(name="Iron sphere flux spectrum")
    flux_tally.filters = [cell_filter, energy_filter]
    flux_tally.scores = ["flux"]
    # "flux" uses the track-length estimator: sum of track lengths in the
    # cell divided by the cell volume, per source particle.

    tallies = openmc.Tallies([leakage_tally, flux_tally])

    # =========================================================================
    # Assemble the Model
    # =========================================================================
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )

    return model


# =============================================================================
# Main entry point
# =============================================================================
if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # Command-line arguments
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="OKTAVIAN Iron Sphere benchmark -- OpenMC fixed-source model",
    )
    parser.add_argument(
        "--particles",
        type=int,
        default=1_000_000,
        help="Total number of source particles to simulate (default: 1000000)",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        default=False,
        help="Run the simulation immediately after exporting XML files",
    )
    args = parser.parse_args()

    # Build the model
    print("=" * 70)
    print("OKTAVIAN Iron Sphere Benchmark")
    print("=" * 70)
    print(f"  Sphere radius : 50.32 cm")
    print(f"  Material      : Carbon steel (7.874 g/cc)")
    print(f"  Source energy  : 14.1 MeV (D-T)")
    print(f"  Particles     : {args.particles:,}")
    print("=" * 70)

    model = build_model(particles=args.particles)

    # Export XML input files (materials.xml, geometry.xml, settings.xml, tallies.xml)
    model.export_to_xml()
    print("\nXML files exported successfully.")

    # Optionally run the simulation
    if args.run:
        print("\nStarting OpenMC simulation...")
        statepoint_path = model.run()
        print(f"\nSimulation complete. Statepoint: {statepoint_path}")
