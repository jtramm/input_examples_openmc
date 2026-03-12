#!/usr/bin/env python3
"""
VERA Core Physics Benchmark - Problem 2A: Single Fuel Assembly at Hot Zero Power
=================================================================================

This script builds an OpenMC model for VERA (Virtual Environment for Reactor
Applications) Benchmark Problem 2A, as defined in:

    CASL-U-2012-0131-004
    "VERA Core Physics Benchmark Progression Problem Specifications"
    A. Godfrey, Oak Ridge National Laboratory (2014)

VERA was developed by the Consortium for Advanced Simulation of Light Water
Reactors (CASL), a U.S. Department of Energy Innovation Hub. The progression
problems provide a graded set of benchmarks for validating reactor physics
codes, starting from simple pin cells and building up to full-core models.

Problem 2A: Single Fuel Assembly (2D, Hot Zero Power)
------------------------------------------------------
This problem models a full 17x17 Westinghouse-type PWR fuel assembly with
reflective boundary conditions on all faces, simulating an infinite lattice
of identical assemblies. It is representative of the fuel assemblies used
in the Watts Bar Nuclear Unit 1 reactor (a Westinghouse 4-loop PWR).

The assembly contains three types of pin locations:
  - 264 fuel pins: UO2 fuel pellet + He gap + Zircaloy-4 clad + moderator
  -  24 guide tubes: empty (water-filled) tubes where control rods can be
    inserted. In this problem, control rods are withdrawn, so the guide
    tubes are filled with moderator water.
  -   1 instrument tube: located at the center of the assembly (position
    8,8 in 0-indexed coordinates). Houses in-core neutron detectors.
    Similar to a guide tube but with slightly different dimensions.

Together: 264 + 24 + 1 = 289 = 17 x 17 positions.

Hot Zero Power (HZP) conditions:
  - Fuel temperature:      565 K
  - Moderator temperature: 565 K
  - Moderator density:     0.743 g/cc
  - Soluble boron:         1300 ppm
  - No control rods inserted (all rods out, or "ARO")

The pin pitch is 1.26 cm, giving an active lattice width of 17 * 1.26 =
21.42 cm. The assembly pitch is 21.50 cm, which includes a 0.04 cm
inter-assembly half-gap on each side (total gap of 0.08 cm between
adjacent assemblies). This gap is filled with moderator water.

Reference k-effective (CE KENO-VI, ENDF/B-VII.0): 1.18703 +/- 0.00009
(from Table 11 of CASL-U-2012-0131-004)

Usage:
    python model.py [--particles N] [--batches N] [--inactive N]
                    [--small-tallies] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# GUIDE TUBE AND INSTRUMENT TUBE POSITIONS
# =============================================================================
# The standard Westinghouse 17x17 fuel assembly layout has 24 guide tube
# positions and 1 instrument tube position. These positions are given in
# (row, column) format with 0-based indexing, where (0,0) is the top-left
# corner of the assembly when viewed from above.
#
# The guide tubes are symmetrically placed in the assembly. During power
# operation, 24 of these tubes house control rod fingers (in rodded
# assemblies) or are simply water-filled (in unrodded assemblies). The
# center position (8,8) holds the instrument thimble for in-core flux
# measurement.
#
# In this HZP benchmark with all rods out, all guide tubes are water-filled.

GUIDE_TUBE_POSITIONS = [
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5),          (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
]

# The instrument tube is at the very center of the assembly.
INSTRUMENT_TUBE_POSITION = (8, 8)


def build_model(particles=50000, batches=200, inactive=50, small_tallies=False):
    """
    Construct the OpenMC model for VERA Problem 2A.

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch.
    batches : int
        Total number of batches (active + inactive).
    inactive : int
        Number of inactive (discarded) batches for source convergence.
    small_tallies : bool
        If True, skip the pin power mesh tally (faster runtime).

    Returns
    -------
    openmc.Model
        The complete OpenMC model ready for export or execution.
    """

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # All number densities are in atoms/barn-cm, taken directly from the
    # VERA benchmark specification (Table 3). These are identical to
    # Problem 1A since the fuel enrichment and moderator conditions are
    # the same.

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

    # ---- Helium Gap ----
    # A thin gas-filled gap exists between the fuel pellet and the inner
    # surface of the cladding. It accommodates fuel swelling and fission gas
    # release. At beginning of life it's filled with helium at low pressure.
    gap = openmc.Material(name='Helium Gap')
    gap.set_density('g/cc', 0.001598)  # helium at ~1 atm, 565K
    gap.add_nuclide('He4', 1.0)

    # ---- Zircaloy-4 Cladding ----
    # Zircaloy-4 is the standard PWR cladding material. It is a zirconium
    # alloy chosen for its low neutron absorption cross section, good
    # corrosion resistance, and mechanical strength at operating temperatures.
    # Used for both fuel rod cladding and guide/instrument tube walls.
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

    # Collect all materials into a Materials object
    materials = openmc.Materials([fuel, gap, clad, moderator])

    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # The geometry uses OpenMC's lattice feature to arrange 17x17 pin
    # universes into the fuel assembly. We define three types of pin
    # universes (fuel pin, guide tube, instrument tube) and then place
    # them into a RectLattice.

    # --- Pin pitch ---
    # Center-to-center distance between adjacent pins in the lattice.
    pin_pitch = 1.26  # cm

    # --- Assembly dimensions ---
    # The assembly pitch is 21.50 cm. The active lattice region is
    # 17 * 1.26 = 21.42 cm, leaving (21.50 - 21.42) / 2 = 0.04 cm
    # of inter-assembly water gap on each side.
    assembly_pitch = 21.50  # cm

    # =====================================================================
    # FUEL PIN UNIVERSE
    # =====================================================================
    # The fuel pin consists of concentric cylinders: fuel pellet, helium
    # gap, Zircaloy-4 cladding, and surrounding moderator water.
    #
    # Cross-section (not to scale):
    #
    #     Moderator (water)
    #     +---------------------------+
    #     |                           |
    #     |     Clad (Zirc-4)         |
    #     |     +-----------------+   |
    #     |     |   Gap (He)      |   |
    #     |     |   +---------+   |   |
    #     |     |   | Fuel    |   |   |
    #     |     |   | (UO2)   |   |   |
    #     |     |   +---------+   |   |
    #     |     +-----------------+   |
    #     +---------------------------+

    # Fuel pellet outer radius: 0.4096 cm
    fuel_or = openmc.ZCylinder(r=0.4096, name='Fuel pellet OR')
    # Cladding inner radius: 0.418 cm (gap is 84 microns thick)
    clad_ir = openmc.ZCylinder(r=0.418, name='Fuel clad IR')
    # Cladding outer radius: 0.475 cm (clad wall is 570 microns thick)
    clad_or = openmc.ZCylinder(r=0.475, name='Fuel clad OR')

    # Fuel pellet cell: inside the fuel outer radius
    fuel_cell = openmc.Cell(name='Fuel pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or
    fuel_cell.temperature = 565.0  # HZP fuel temperature [K]

    # Helium gap cell: annular region between fuel and clad
    gap_cell = openmc.Cell(name='Helium gap')
    gap_cell.fill = gap
    gap_cell.region = +fuel_or & -clad_ir
    gap_cell.temperature = 565.0

    # Cladding cell: annular region of Zircaloy-4
    clad_cell = openmc.Cell(name='Fuel cladding')
    clad_cell.fill = clad
    clad_cell.region = +clad_ir & -clad_or
    clad_cell.temperature = 565.0

    # Moderator cell: everything outside the cladding within this universe.
    # When placed in a lattice, OpenMC automatically clips each universe
    # to its lattice cell boundaries, so we don't need to define the
    # outer rectangular boundary explicitly.
    mod_cell = openmc.Cell(name='Fuel pin moderator')
    mod_cell.fill = moderator
    mod_cell.region = +clad_or
    mod_cell.temperature = 565.0

    # Assemble the fuel pin universe
    fuel_pin_universe = openmc.Universe(
        name='Fuel Pin',
        cells=[fuel_cell, gap_cell, clad_cell, mod_cell]
    )

    # =====================================================================
    # GUIDE TUBE UNIVERSE
    # =====================================================================
    # Guide tubes are larger-diameter Zircaloy-4 tubes that are normally
    # empty (filled with moderator water) when control rods are withdrawn.
    # In this problem (all rods out), the guide tubes contain only water
    # inside and outside the tube wall.
    #
    # Guide tube dimensions:
    #   Inner radius: 0.561 cm
    #   Outer radius: 0.602 cm
    #   Wall thickness: 0.041 cm = 410 microns
    #
    # These are larger than fuel pin cladding (OR = 0.475 cm) to allow
    # insertion of control rod absorber pins.

    # Guide tube inner radius
    gt_ir = openmc.ZCylinder(r=0.561, name='Guide tube IR')
    # Guide tube outer radius
    gt_or = openmc.ZCylinder(r=0.602, name='Guide tube OR')

    # Water inside the guide tube (where control rods would go)
    gt_inner_cell = openmc.Cell(name='GT inner water')
    gt_inner_cell.fill = moderator
    gt_inner_cell.region = -gt_ir
    gt_inner_cell.temperature = 565.0

    # Guide tube wall (Zircaloy-4)
    gt_wall_cell = openmc.Cell(name='GT wall')
    gt_wall_cell.fill = clad
    gt_wall_cell.region = +gt_ir & -gt_or
    gt_wall_cell.temperature = 565.0

    # Water outside the guide tube
    gt_outer_cell = openmc.Cell(name='GT outer water')
    gt_outer_cell.fill = moderator
    gt_outer_cell.region = +gt_or
    gt_outer_cell.temperature = 565.0

    # Assemble the guide tube universe
    guide_tube_universe = openmc.Universe(
        name='Guide Tube',
        cells=[gt_inner_cell, gt_wall_cell, gt_outer_cell]
    )

    # =====================================================================
    # INSTRUMENT TUBE UNIVERSE
    # =====================================================================
    # The instrument tube at the center of the assembly (position 8,8)
    # houses in-core neutron flux detectors. It has slightly different
    # dimensions than the guide tubes:
    #   Inner radius: 0.559 cm
    #   Outer radius: 0.605 cm
    #   Wall thickness: 0.046 cm = 460 microns
    #
    # In this benchmark, the instrument tube interior is modeled as
    # containing moderator water (the detector is not explicitly modeled
    # as it has negligible neutronic effect in this context).

    # Instrument tube inner radius
    it_ir = openmc.ZCylinder(r=0.559, name='Instrument tube IR')
    # Instrument tube outer radius
    it_or = openmc.ZCylinder(r=0.605, name='Instrument tube OR')

    # Water inside the instrument tube
    it_inner_cell = openmc.Cell(name='IT inner water')
    it_inner_cell.fill = moderator
    it_inner_cell.region = -it_ir
    it_inner_cell.temperature = 565.0

    # Instrument tube wall (Zircaloy-4)
    it_wall_cell = openmc.Cell(name='IT wall')
    it_wall_cell.fill = clad
    it_wall_cell.region = +it_ir & -it_or
    it_wall_cell.temperature = 565.0

    # Water outside the instrument tube
    it_outer_cell = openmc.Cell(name='IT outer water')
    it_outer_cell.fill = moderator
    it_outer_cell.region = +it_or
    it_outer_cell.temperature = 565.0

    # Assemble the instrument tube universe
    instrument_tube_universe = openmc.Universe(
        name='Instrument Tube',
        cells=[it_inner_cell, it_wall_cell, it_outer_cell]
    )

    # =====================================================================
    # 17x17 ASSEMBLY LATTICE
    # =====================================================================
    # Build the 17x17 lattice by populating a 2D array of universes.
    # Most positions get the fuel pin universe; the guide tube and
    # instrument tube positions are overridden.
    #
    # The lattice is centered at the origin. Each element spans one
    # pin pitch (1.26 cm x 1.26 cm).

    lattice = openmc.RectLattice(name='17x17 Fuel Assembly Lattice')

    # The lattice pitch is the pin-to-pin spacing in both x and y.
    lattice.pitch = (pin_pitch, pin_pitch)

    # Center the lattice at the origin. The lattice extends from
    # -8.5*pitch to +8.5*pitch in both x and y.
    lattice.lower_left = (-17 * pin_pitch / 2.0, -17 * pin_pitch / 2.0)

    # Build the 17x17 universe array.
    # OpenMC RectLattice universes are specified as a 2D list where
    # universes[i][j] corresponds to row i (from top), column j (from left).
    # Row 0 is the top row (highest y), row 16 is the bottom row (lowest y).
    universes = []
    for row in range(17):
        row_list = []
        for col in range(17):
            if (row, col) == INSTRUMENT_TUBE_POSITION:
                # Center position: instrument tube
                row_list.append(instrument_tube_universe)
            elif (row, col) in GUIDE_TUBE_POSITIONS:
                # Guide tube position (control rod location, rods withdrawn)
                row_list.append(guide_tube_universe)
            else:
                # Standard fuel pin position
                row_list.append(fuel_pin_universe)
        universes.append(row_list)

    lattice.universes = universes

    # Set the outer universe for the lattice. Any point outside the
    # lattice boundaries will be assigned to this universe. We use
    # moderator water since the inter-assembly gap is water-filled.
    outer_universe = openmc.Universe(name='Outer (water)')
    outer_cell = openmc.Cell(name='Outer water cell', fill=moderator)
    outer_cell.temperature = 565.0
    outer_universe.add_cell(outer_cell)
    lattice.outer = outer_universe

    # =====================================================================
    # ASSEMBLY CELL AND ROOT UNIVERSE
    # =====================================================================
    # The assembly is bounded by a rectangular prism with reflective
    # boundary conditions. The prism dimensions match the assembly pitch
    # (21.50 cm) to properly represent the inter-assembly gap.
    #
    # The 0.04 cm water gap between the lattice edge (at 17*1.26/2 =
    # 10.71 cm) and the assembly boundary (at 21.50/2 = 10.75 cm) is
    # automatically filled by the lattice's outer universe (moderator).

    half_assembly = assembly_pitch / 2.0  # 10.75 cm

    # Assembly boundary planes with reflective BCs.
    # Reflective boundaries simulate an infinite 2D array of identical
    # assemblies - appropriate for Problem 2A which has no inter-assembly
    # heterogeneity.
    left = openmc.XPlane(x0=-half_assembly, boundary_type='reflective',
                         name='Assembly left boundary')
    right = openmc.XPlane(x0=+half_assembly, boundary_type='reflective',
                          name='Assembly right boundary')
    bottom = openmc.YPlane(y0=-half_assembly, boundary_type='reflective',
                           name='Assembly bottom boundary')
    top = openmc.YPlane(y0=+half_assembly, boundary_type='reflective',
                        name='Assembly top boundary')

    # Axial boundaries - 2D problem approximated with a thin slice.
    # Reflective BCs in z simulate infinite axial extent.
    z_min = openmc.ZPlane(z0=-0.5, boundary_type='reflective',
                          name='Axial bottom (reflective)')
    z_max = openmc.ZPlane(z0=+0.5, boundary_type='reflective',
                          name='Axial top (reflective)')

    # Create the main assembly cell containing the lattice
    assembly_cell = openmc.Cell(name='Fuel Assembly')
    assembly_cell.fill = lattice
    assembly_cell.region = +left & -right & +bottom & -top & +z_min & -z_max

    # Create root universe and geometry
    root_universe = openmc.Universe(name='Root', cells=[assembly_cell])
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    # Configure the eigenvalue (k-effective) calculation.

    settings = openmc.Settings()

    # Eigenvalue mode: solve the neutron transport eigenvalue problem.
    settings.run_mode = 'eigenvalue'

    # Number of neutron histories per batch.
    # 50,000 is a reasonable default for a fuel assembly.
    settings.particles = particles

    # Total number of batches (active + inactive).
    settings.batches = batches

    # Inactive batches for fission source convergence.
    # Assembly problems need more inactive batches than pin cells because
    # the source must equilibrate across all 264 fuel pins.
    settings.inactive = inactive

    # Temperature settings: use nearest available temperature data.
    # The cross section library may not have data at exactly 565K,
    # so we tell OpenMC to use the nearest available temperature.
    settings.temperature = {'method': 'nearest', 'tolerance': 1000}

    # Initial source distribution: uniform box spanning the fuel region.
    # The box covers the entire lattice extent so that all fuel pins
    # receive initial source particles. OpenMC will reject any samples
    # that don't land in fissile material.
    lattice_half = 17 * pin_pitch / 2.0  # 10.71 cm
    settings.source = [openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-lattice_half, -lattice_half, -0.5),
            upper_right=(+lattice_half, +lattice_half, +0.5)
        )
    )]

    # =========================================================================
    # TALLIES (optional pin power distribution)
    # =========================================================================
    tallies = openmc.Tallies()

    if not small_tallies:
        # Create a regular mesh tally over the assembly to capture the
        # pin-by-pin fission rate (power) distribution. Each mesh element
        # covers one pin cell, giving a 17x17 map of relative pin powers.
        #
        # Pin power distributions are a key validation quantity for assembly
        # calculations. They test the spatial fidelity of the transport solution
        # and are important for ensuring fuel pins don't exceed thermal limits.

        pin_mesh = openmc.RegularMesh(name='Pin power mesh')
        pin_mesh.dimension = (17, 17)  # one mesh element per pin position
        pin_mesh.lower_left = (-lattice_half, -lattice_half)
        pin_mesh.upper_right = (+lattice_half, +lattice_half)

        # Tally fission rate (proportional to power) on the pin mesh
        pin_power_tally = openmc.Tally(name='Pin Powers')
        pin_power_tally.filters = [openmc.MeshFilter(pin_mesh)]
        pin_power_tally.scores = ['fission']
        tallies.append(pin_power_tally)

    # =========================================================================
    # ASSEMBLE AND RETURN MODEL
    # =========================================================================
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings, tallies=tallies)
    return model


def main():
    """Parse command-line arguments, build the model, and export/run."""
    parser = argparse.ArgumentParser(
        description='VERA Problem 2A: 17x17 PWR Fuel Assembly (OpenMC)'
    )
    parser.add_argument('--particles', type=int, default=50000,
                        help='Neutron histories per batch (default: 50000)')
    parser.add_argument('--batches', type=int, default=200,
                        help='Total active batches (default: 200)')
    parser.add_argument('--inactive', type=int, default=50,
                        help='Inactive batches for source convergence (default: 50)')
    parser.add_argument('--small-tallies', action='store_true',
                        help='Skip pin power mesh tally (faster runtime)')
    parser.add_argument('--run', action='store_true',
                        help='Run OpenMC after exporting XML files')
    args = parser.parse_args()

    # Build the model
    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
        small_tallies=args.small_tallies
    )

    # Export XML input files (geometry.xml, materials.xml, settings.xml, tallies.xml)
    model.export_to_xml()
    print("Exported OpenMC XML input files.")

    # Print a summary of the assembly layout
    n_fuel = 264
    n_gt = len(GUIDE_TUBE_POSITIONS)
    n_it = 1
    print(f"Assembly layout: {n_fuel} fuel pins, {n_gt} guide tubes, "
          f"{n_it} instrument tube = {n_fuel + n_gt + n_it} positions")

    # Optionally run the simulation
    if args.run:
        print(f"Running OpenMC: {args.particles} particles/batch, "
              f"{args.batches} active + {args.inactive} inactive batches")
        model.run()


if __name__ == '__main__':
    main()
