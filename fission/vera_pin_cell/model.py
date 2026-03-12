#!/usr/bin/env python3
"""
VERA Core Physics Benchmark - Problem 1A: Single Pin Cell at Hot Zero Power
============================================================================

This script builds an OpenMC model for VERA (Virtual Environment for Reactor
Applications) Benchmark Problem 1A, as defined in:

    CASL-U-2012-0131-004
    "VERA Core Physics Benchmark Progression Problem Specifications"
    A. Godfrey, Oak Ridge National Laboratory (2014)

VERA was developed by the Consortium for Advanced Simulation of Light Water
Reactors (CASL), a U.S. Department of Energy Innovation Hub. The progression
problems provide a graded set of benchmarks for validating reactor physics
codes, starting from simple pin cells and building up to full-core models.

Problem 1A: Single Fuel Pin Cell
---------------------------------
This is the simplest problem in the VERA suite. It models a single
Westinghouse 17x17 PWR fuel pin surrounded by borated water, with
reflective boundary conditions on all six faces to simulate an infinite
lattice of identical pins. The conditions are Hot Zero Power (HZP):

  - Fuel temperature:      565 K
  - Moderator temperature: 565 K
  - Moderator density:     0.743 g/cc
  - Soluble boron:         1300 ppm

Physical geometry (radial, from center outward):
  1. UO2 fuel pellet       r = 0.4096 cm
  2. Helium gap            r = 0.418  cm  (between pellet and clad)
  3. Zircaloy-4 cladding   r = 0.475  cm  (outer radius)
  4. Borated water         fills the rest of the pin cell

The pin pitch is 1.26 cm (center-to-center distance in the 17x17 lattice).

Reference k-effective (CE KENO-VI, ENDF/B-VII.0): 1.187038 +/- 0.000054

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import openmc


def build_model(particles=50000, batches=200, inactive=50):
    """
    Construct the OpenMC model for VERA Problem 1A.

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch.
    batches : int
        Total number of batches (active + inactive).
    inactive : int
        Number of inactive (discarded) batches for source convergence.

    Returns
    -------
    openmc.Model
        The complete OpenMC model ready for export or execution.
    """

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # All number densities are in atoms/barn-cm, taken directly from the
    # VERA benchmark specification (Table 3).

    # ---- UO2 Fuel (3.1 wt% enriched) ----
    # Density: 10.257 g/cc (theoretical density of UO2 pellet)
    # Enrichment: 3.1 w/o U-235
    # The benchmark specifies explicit isotopic number densities including
    # trace amounts of U-234 and U-236 from the enrichment process.
    fuel = openmc.Material(name='UO2 Fuel 3.1%')
    fuel.set_density('sum')  # density determined by sum of number densities
    fuel.add_nuclide('U234', 6.11864e-06)   # trace from enrichment
    fuel.add_nuclide('U235', 7.18132e-04)   # fissile isotope (3.1 w/o)
    fuel.add_nuclide('U236', 3.29861e-06)   # trace from enrichment
    fuel.add_nuclide('U238', 2.21546e-02)   # fertile isotope (bulk)
    fuel.add_nuclide('O16',  4.57642e-02)   # oxygen in UO2
    # Temperature set on cell, not material

    # ---- Helium Gap ----
    # A thin gas-filled gap exists between the fuel pellet and the inner
    # surface of the cladding. It accommodates fuel swelling and fission gas
    # release. At beginning of life it's filled with helium at low pressure.
    # We model it as void (near-zero density helium) since it has negligible
    # neutronic effect.
    gap = openmc.Material(name='Helium Gap')
    gap.set_density('g/cc', 0.001598)  # helium at ~1 atm, 565K
    gap.add_nuclide('He4', 1.0)

    # ---- Zircaloy-4 Cladding ----
    # Zircaloy-4 is the standard PWR cladding material. It is a zirconium
    # alloy chosen for its low neutron absorption cross section, good
    # corrosion resistance, and mechanical strength at operating temperatures.
    # Composition (typical weight percents):
    #   Zr: balance (~98.2%), Sn: 1.5%, Fe: 0.2%, Cr: 0.1%
    clad = openmc.Material(name='Zircaloy-4')
    clad.set_density('g/cc', 6.56)
    clad.add_element('Zr', 0.982, 'wo')  # zirconium (balance)
    clad.add_element('Sn', 0.015, 'wo')  # tin - improves strength
    clad.add_element('Fe', 0.002, 'wo')  # iron - improves corrosion resistance
    clad.add_element('Cr', 0.001, 'wo')  # chromium - improves corrosion resistance

    # ---- Borated Water Moderator ----
    # The moderator is light water with dissolved boric acid. Soluble boron
    # (specifically B-10, which has a large thermal absorption cross section)
    # is used for reactivity control. At 1300 ppm boron and 0.743 g/cc:
    moderator = openmc.Material(name='Borated Water')
    moderator.set_density('sum')  # density from sum of number densities
    moderator.add_nuclide('H1',  4.96224e-02)  # hydrogen in water
    moderator.add_nuclide('O16', 2.48112e-02)  # oxygen in water
    moderator.add_nuclide('B10', 1.07070e-05)  # soluble boron (absorber)
    moderator.add_nuclide('B11', 4.30971e-05)  # soluble boron (inert)
    # Add S(alpha,beta) thermal scattering data for hydrogen in water.
    # This is critical for correctly modeling neutron thermalization.
    # Without it, the free-gas scattering kernel would be used for hydrogen,
    # which poorly represents the bound-state effects in water molecules
    # and would give an incorrect thermal spectrum and biased k-eff.
    moderator.add_s_alpha_beta('c_H_in_H2O')

    # Collect all materials
    materials = openmc.Materials([fuel, gap, clad, moderator])

    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # The geometry consists of concentric cylinders (fuel, gap, clad) centered
    # in a rectangular prism (the pin cell). All surfaces are along the z-axis.

    # --- Radial surfaces (cylinders along z-axis) ---

    # Fuel pellet outer surface at r = 0.4096 cm
    # This is the outer radius of the UO2 fuel pellet.
    fuel_or = openmc.ZCylinder(r=0.4096, name='Fuel pellet outer radius')

    # Cladding inner surface at r = 0.418 cm
    # The gap between fuel_or (0.4096) and clad_ir (0.418) is only 84 microns
    # wide - just enough for a thin helium gas gap.
    clad_ir = openmc.ZCylinder(r=0.418, name='Cladding inner radius')

    # Cladding outer surface at r = 0.475 cm
    # The cladding wall thickness is 0.475 - 0.418 = 0.057 cm = 570 microns.
    clad_or = openmc.ZCylinder(r=0.475, name='Cladding outer radius')

    # --- Pin cell boundary (rectangular prism) ---
    # The pin pitch is 1.26 cm, so the cell extends from -0.63 to +0.63 cm
    # in both x and y. Reflective boundary conditions simulate an infinite
    # lattice of identical pins - neutrons that hit the boundary are reflected
    # back, as if surrounded by identical pin cells on all sides.
    pitch = 1.26  # cm, center-to-center pin spacing in 17x17 lattice
    half_pitch = pitch / 2.0

    # Left boundary at x = -0.63 cm (reflective)
    left = openmc.XPlane(x0=-half_pitch, boundary_type='reflective',
                         name='Left boundary')
    # Right boundary at x = +0.63 cm (reflective)
    right = openmc.XPlane(x0=+half_pitch, boundary_type='reflective',
                          name='Right boundary')
    # Bottom boundary at y = -0.63 cm (reflective)
    bottom = openmc.YPlane(y0=-half_pitch, boundary_type='reflective',
                           name='Bottom boundary')
    # Top boundary at y = +0.63 cm (reflective)
    top = openmc.YPlane(y0=+half_pitch, boundary_type='reflective',
                        name='Top boundary')

    # --- Axial surfaces ---
    # We model a thin axial slice. The actual VERA problem is 2D (infinite
    # in z), which we approximate with reflective top/bottom planes.
    z_min = openmc.ZPlane(z0=-0.5, boundary_type='reflective',
                          name='Axial bottom (reflective)')
    z_max = openmc.ZPlane(z0=+0.5, boundary_type='reflective',
                          name='Axial top (reflective)')

    # --- Cell definitions ---
    # Each cell is a region of space filled with a material. Regions are
    # defined by the intersection (&) and complement (~) of half-spaces.

    # Cell 1: Fuel pellet
    # Region: inside the fuel outer radius cylinder, between axial planes
    # Material: UO2 fuel
    fuel_cell = openmc.Cell(name='UO2 Fuel Pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or & +z_min & -z_max
    fuel_cell.temperature = 565.0  # Hot Zero Power fuel temperature [K]

    # Cell 2: Helium gap
    # Region: between fuel outer radius and cladding inner radius
    # Material: low-density helium (essentially void)
    gap_cell = openmc.Cell(name='Helium Gap')
    gap_cell.fill = gap
    gap_cell.region = +fuel_or & -clad_ir & +z_min & -z_max
    gap_cell.temperature = 565.0

    # Cell 3: Zircaloy-4 cladding
    # Region: between cladding inner and outer radii
    # Material: Zircaloy-4
    clad_cell = openmc.Cell(name='Zircaloy-4 Cladding')
    clad_cell.fill = clad
    clad_cell.region = +clad_ir & -clad_or & +z_min & -z_max
    clad_cell.temperature = 565.0

    # Cell 4: Moderator (borated water)
    # Region: outside the cladding, inside the pin cell box
    # Material: borated light water at 0.743 g/cc with 1300 ppm boron
    mod_cell = openmc.Cell(name='Borated Water Moderator')
    mod_cell.fill = moderator
    mod_cell.region = +clad_or & +left & -right & +bottom & -top & +z_min & -z_max
    mod_cell.temperature = 565.0  # Hot Zero Power moderator temperature [K]

    # Create the root universe containing all four cells
    root_universe = openmc.Universe(cells=[fuel_cell, gap_cell,
                                           clad_cell, mod_cell])

    # Build the geometry from the root universe
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    # Configure the eigenvalue (k-effective) calculation.

    settings = openmc.Settings()

    # Eigenvalue mode: solve the neutron transport eigenvalue problem to
    # find k-effective, the ratio of neutron production to loss per generation.
    settings.run_mode = 'eigenvalue'

    # Number of neutron histories (particles) simulated per batch.
    # More particles = lower statistical uncertainty but longer runtime.
    # 50,000 is a reasonable balance for this simple pin cell.
    settings.particles = particles

    # Total number of active batches. k-eff is tallied over these batches
    # and averaged to get the final result. More batches = lower uncertainty.
    settings.batches = batches

    # Inactive batches: these are run first to converge the fission source
    # distribution before tallying begins. For a simple pin cell with
    # reflective boundaries, the source converges quickly, but 50 inactive
    # batches provides a comfortable margin.
    settings.inactive = inactive

    # Temperature settings: use nearest available temperature data.
    # The cross section library may not have data at exactly 565K,
    # so we tell OpenMC to use the nearest available temperature.
    settings.temperature = {'method': 'nearest', 'tolerance': 1000}

    # Initial source distribution: uniform in the fuel region.
    # We define a box source that covers the fuel pellet. OpenMC will
    # sample initial neutron positions uniformly within this box and
    # reject any that fall outside the fissile region.
    settings.source = [openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-0.4096, -0.4096, -0.5),
            upper_right=(0.4096, 0.4096, 0.5)
        )
    )]

    # =========================================================================
    # ASSEMBLE AND RETURN MODEL
    # =========================================================================
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings)
    return model


def main():
    """Parse command-line arguments, build the model, and export/run."""
    parser = argparse.ArgumentParser(
        description='VERA Problem 1A: PWR Pin Cell (OpenMC)'
    )
    parser.add_argument('--particles', type=int, default=50000,
                        help='Neutron histories per batch (default: 50000)')
    parser.add_argument('--batches', type=int, default=200,
                        help='Total active batches (default: 200)')
    parser.add_argument('--inactive', type=int, default=50,
                        help='Inactive batches for source convergence (default: 50)')
    parser.add_argument('--run', action='store_true',
                        help='Run OpenMC after exporting XML files')
    args = parser.parse_args()

    # Build the model
    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive
    )

    # Export XML input files (geometry.xml, materials.xml, settings.xml)
    model.export_to_xml()
    print("Exported OpenMC XML input files.")

    # Optionally run the simulation
    if args.run:
        print(f"Running OpenMC: {args.particles} particles/batch, "
              f"{args.batches} active + {args.inactive} inactive batches")
        model.run()


if __name__ == '__main__':
    main()
