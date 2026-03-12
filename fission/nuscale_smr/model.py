#!/usr/bin/env python3
"""
McSAFER NuScale-like SMR Benchmark - Single Fuel Assembly (2D Infinite Lattice)
================================================================================

This script builds an OpenMC model for a NuScale-like Small Modular Reactor
(SMR) fuel assembly, based on the McSAFER benchmark project specifications:

    "Neutronics Benchmark of the NuScale-like SMR in the McSAFER Project"
    J. Bousquet et al., J. Nucl. Eng. 2025, 6(4), 44
    https://www.mdpi.com/2673-4362/6/4/44

NuScale Power Module Overview
------------------------------
The NuScale SMR is a light water reactor designed as a self-contained power
module producing approximately 77 MWe (250 MWth per module). The entire
primary system -- reactor core, steam generators, and pressurizer -- is
housed inside a single integrated vessel, which is itself enclosed in a
steel containment submerged in a large below-grade pool. Up to 12 modules
can operate in a single NuScale VOYGR power plant.

Key design features that distinguish the NuScale from conventional PWRs:
  - Compact core: 37 fuel assemblies (vs. ~157-193 in a standard PWR)
  - Short active fuel length: 200 cm (vs. ~366 cm in a standard PWR)
  - Natural circulation cooling: no reactor coolant pumps needed
  - BORON-FREE operation: unlike conventional PWRs which use dissolved
    boric acid in the coolant for reactivity control, the NuScale design
    relies entirely on burnable absorbers (Gd2O3) and control rods. This
    simplifies the chemical volume control system and eliminates boron
    dilution transient concerns.

The core uses standard 17x17 Westinghouse-type fuel assemblies with
264 fuel rods, 24 guide tubes, and 1 instrument tube. Seven different
assembly types are loaded, varying in enrichment and gadolinium content.

This Model
----------
We implement the simplest assembly type: 3.1 wt% U-235 enriched UO2
without gadolinium burnable absorber pins. Two model options are available:

  'pin-cell'  -- A single fuel pin cell with reflective BCs, representing
                 an infinite lattice of identical pins. Useful for quick
                 verification of material definitions and cross sections.

  'assembly'  -- The full 17x17 fuel assembly with reflective BCs on all
                 faces, representing an infinite 2D lattice of identical
                 assemblies. This is the primary benchmark configuration.

Operating Conditions (Hot Full Power)
--------------------------------------
  - Fuel temperature:      900 K (higher than HZP due to power generation)
  - Moderator temperature: 557 K (~284 degC, NuScale operating temperature)
  - Moderator pressure:    12.76 MPa (NuScale operating pressure)
  - Moderator density:     0.7527 g/cc (subcooled water at 557K, 12.76 MPa)
  - Soluble boron:         0 ppm (BORON-FREE design!)

Geometry Dimensions
--------------------
  - Pin pitch:              1.2598 cm
  - Assembly pitch:         21.5036 cm (includes inter-assembly water gap)
  - Active fuel length:     200 cm (modeled as 2D with reflective axial BCs)
  - Fuel pellet OR:         0.4096 cm
  - Fuel clad IR:           0.418 cm (He gap = 84 microns)
  - Fuel clad OR:           0.475 cm
  - Guide tube IR:          0.561 cm
  - Guide tube OR:          0.602 cm
  - Instrument tube IR:     0.559 cm
  - Instrument tube OR:     0.605 cm

Reference Results
------------------
  Serpent 2 reference k-eff for this assembly type:
    k-eff ~ 1.18 (exact value depends on cross section library and
    depletion state; fresh fuel at beginning of life)
  OpenMC vs Serpent difference: ~355 pcm reported in the benchmark paper

Usage:
    python model.py [--particles N] [--batches N] [--inactive N]
                    [--model {pin-cell,assembly}] [--run]
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
# The NuScale SMR uses the same 17x17 layout as conventional Westinghouse
# assemblies. The 24 guide tube locations can house control rods or remain
# water-filled. Since this is a boron-free design, the control rods and
# burnable absorbers carry the full reactivity control burden that would
# otherwise be shared with soluble boron.
#
# In this benchmark, we model the assembly with all control rods withdrawn
# (all rods out, ARO), so all guide tubes are water-filled.

GUIDE_TUBE_POSITIONS = [
    # Row 2: three guide tubes
    (2, 5), (2, 8), (2, 11),
    # Row 3: two guide tubes
    (3, 3), (3, 13),
    # Row 5: five guide tubes
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    # Row 8: four guide tubes (center position is instrument tube)
    (8, 2), (8, 5),          (8, 11), (8, 14),
    # Row 11: five guide tubes (mirror of row 5)
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    # Row 13: two guide tubes (mirror of row 3)
    (13, 3), (13, 13),
    # Row 14: three guide tubes (mirror of row 2)
    (14, 5), (14, 8), (14, 11),
]

# The instrument tube is at the very center of the assembly (position 8,8
# in 0-indexed coordinates). It houses in-core neutron flux detectors.
INSTRUMENT_TUBE_POSITION = (8, 8)


def build_pin_cell(particles=50000, batches=200, inactive=50):
    """
    Construct a single NuScale-like SMR fuel pin cell model.

    This is the simplest possible model: one fuel pin with reflective
    boundary conditions on all faces, simulating an infinite array of
    identical fuel pins. Useful for verifying material definitions,
    cross section processing, and basic transport before moving to the
    full assembly model.

    The pin cell consists of:
      - UO2 fuel pellet (3.1% enriched, 0.4096 cm radius)
      - Helium gap (0.418 cm outer radius)
      - Zircaloy-4 cladding (0.475 cm outer radius)
      - Unborated light water moderator (fills the square cell)

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
        The complete OpenMC pin cell model.
    """

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # Material definitions are shared between pin-cell and assembly models.
    # See build_assembly() for detailed material descriptions.

    fuel, gap, clad, moderator = _create_materials()
    materials = openmc.Materials([fuel, gap, clad, moderator])

    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # Single fuel pin cell: concentric cylinders in a rectangular box.

    # --- Pin pitch ---
    # NuScale-like SMR pin pitch is 1.2598 cm, slightly smaller than the
    # standard Westinghouse 17x17 pitch of 1.26 cm.
    pin_pitch = 1.2598  # cm

    # --- Radial surfaces (concentric cylinders) ---
    # Fuel pellet outer radius
    fuel_or = openmc.ZCylinder(r=0.4096, name='Fuel pellet OR')
    # Cladding inner radius (gap between pellet and clad)
    clad_ir = openmc.ZCylinder(r=0.418, name='Fuel clad IR')
    # Cladding outer radius
    clad_or = openmc.ZCylinder(r=0.475, name='Fuel clad OR')

    # --- Cells ---
    # Fuel pellet: innermost region, filled with enriched UO2
    fuel_cell = openmc.Cell(name='Fuel pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or
    fuel_cell.temperature = 900.0  # Hot full power fuel temperature [K]

    # Helium gap: thin annular region between pellet and clad
    gap_cell = openmc.Cell(name='Helium gap')
    gap_cell.fill = gap
    gap_cell.region = +fuel_or & -clad_ir
    gap_cell.temperature = 900.0

    # Zircaloy-4 cladding: structural tube around the fuel
    clad_cell = openmc.Cell(name='Fuel cladding')
    clad_cell.fill = clad
    clad_cell.region = +clad_ir & -clad_or
    clad_cell.temperature = 600.0  # Clad temperature between fuel and moderator

    # Moderator: unborated light water surrounding the pin
    # Note: NO BORON -- this is a key NuScale design feature!
    mod_cell = openmc.Cell(name='Moderator')
    mod_cell.fill = moderator
    mod_cell.region = +clad_or
    mod_cell.temperature = 557.0  # NuScale moderator temperature [K]

    # Create pin cell universe
    pin_universe = openmc.Universe(
        name='Pin Cell',
        cells=[fuel_cell, gap_cell, clad_cell, mod_cell]
    )

    # --- Pin cell boundary ---
    # Rectangular prism with reflective BCs, one pin pitch on each side.
    half_pitch = pin_pitch / 2.0

    left = openmc.XPlane(x0=-half_pitch, boundary_type='reflective',
                         name='Left boundary')
    right = openmc.XPlane(x0=+half_pitch, boundary_type='reflective',
                          name='Right boundary')
    bottom = openmc.YPlane(y0=-half_pitch, boundary_type='reflective',
                           name='Bottom boundary')
    top = openmc.YPlane(y0=+half_pitch, boundary_type='reflective',
                        name='Top boundary')
    z_min = openmc.ZPlane(z0=-0.5, boundary_type='reflective',
                          name='Axial bottom')
    z_max = openmc.ZPlane(z0=+0.5, boundary_type='reflective',
                          name='Axial top')

    # Main cell containing the pin universe
    main_cell = openmc.Cell(name='Pin Cell')
    main_cell.fill = pin_universe
    main_cell.region = +left & -right & +bottom & -top & +z_min & -z_max

    # Root universe and geometry
    root_universe = openmc.Universe(name='Root', cells=[main_cell])
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = _create_settings(particles, batches, inactive,
                                half_extent=half_pitch)

    # =========================================================================
    # ASSEMBLE MODEL
    # =========================================================================
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings)
    return model


def build_assembly(particles=50000, batches=200, inactive=50):
    """
    Construct the full 17x17 NuScale-like SMR fuel assembly model.

    This is the primary benchmark configuration: a 17x17 fuel assembly with
    reflective boundary conditions on all faces, representing an infinite
    2D lattice of identical assemblies. The assembly contains:
      - 264 UO2 fuel pins (3.1 wt% U-235, no gadolinium)
      -  24 guide tubes (water-filled, control rods withdrawn)
      -   1 instrument tube (water-filled, center position)

    The NuScale SMR assembly layout follows the standard Westinghouse 17x17
    pattern but with slightly different pin pitch (1.2598 cm vs 1.26 cm)
    and a boron-free moderator.

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
        The complete OpenMC assembly model.
    """

    # =========================================================================
    # MATERIALS
    # =========================================================================
    fuel, gap, clad, moderator = _create_materials()
    materials = openmc.Materials([fuel, gap, clad, moderator])

    # =========================================================================
    # GEOMETRY
    # =========================================================================

    # --- Key dimensions ---
    # Pin pitch: center-to-center spacing of fuel rods in the 17x17 lattice.
    # The NuScale-like SMR uses 1.2598 cm, marginally different from the
    # standard Westinghouse 1.26 cm. This small difference affects the
    # moderator-to-fuel ratio and thus the neutron spectrum.
    pin_pitch = 1.2598  # cm

    # Assembly pitch: center-to-center distance between adjacent assemblies.
    # This is larger than 17 * pin_pitch = 21.4166 cm because it includes
    # the inter-assembly water gap:
    #   gap = (21.5036 - 21.4166) / 2 = 0.0435 cm per side
    assembly_pitch = 21.5036  # cm

    # =====================================================================
    # FUEL PIN UNIVERSE
    # =====================================================================
    # The fuel pin is the most common cell type in the assembly (264 of 289
    # positions). It consists of concentric cylinders:
    #
    #   1. UO2 fuel pellet (r = 0.4096 cm)
    #   2. Helium gap      (r = 0.418 cm)  -- accommodates fuel swelling
    #   3. Zircaloy-4 clad (r = 0.475 cm)  -- structural containment
    #   4. Moderator water  -- fills remaining space in lattice cell
    #
    # Cross-section (not to scale):
    #
    #     Moderator (unborated water!)
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
    # Cladding inner radius: 0.418 cm (84 micron He gap)
    clad_ir = openmc.ZCylinder(r=0.418, name='Fuel clad IR')
    # Cladding outer radius: 0.475 cm (570 micron clad wall)
    clad_or = openmc.ZCylinder(r=0.475, name='Fuel clad OR')

    # Fuel pellet cell
    fuel_cell = openmc.Cell(name='Fuel pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or
    fuel_cell.temperature = 900.0  # HFP fuel temperature [K]

    # Helium gap cell
    gap_cell = openmc.Cell(name='Helium gap')
    gap_cell.fill = gap
    gap_cell.region = +fuel_or & -clad_ir
    gap_cell.temperature = 900.0

    # Zircaloy-4 cladding cell
    clad_cell = openmc.Cell(name='Fuel cladding')
    clad_cell.fill = clad
    clad_cell.region = +clad_ir & -clad_or
    clad_cell.temperature = 600.0  # Clad temperature

    # Moderator water cell (NO BORON!)
    # When placed in the lattice, OpenMC clips this universe to the
    # lattice cell boundaries automatically.
    mod_cell = openmc.Cell(name='Fuel pin moderator')
    mod_cell.fill = moderator
    mod_cell.region = +clad_or
    mod_cell.temperature = 557.0  # NuScale moderator temperature [K]

    # Assemble the fuel pin universe
    fuel_pin_universe = openmc.Universe(
        name='Fuel Pin',
        cells=[fuel_cell, gap_cell, clad_cell, mod_cell]
    )

    # =====================================================================
    # GUIDE TUBE UNIVERSE
    # =====================================================================
    # Guide tubes are larger-diameter Zircaloy-4 tubes normally used for
    # control rod insertion. In the NuScale boron-free design, control rods
    # bear a larger share of reactivity control than in conventional PWRs.
    #
    # With all rods withdrawn (ARO), the guide tubes are simply water-filled
    # tubes. Their larger inner diameter (0.561 cm vs. fuel clad OR of
    # 0.475 cm) provides more moderation locally, which creates a slight
    # increase in thermal flux (and hence fission rate) in adjacent fuel pins.
    #
    # Guide tube dimensions:
    #   Inner radius: 0.561 cm
    #   Outer radius: 0.602 cm
    #   Wall thickness: 0.041 cm = 410 microns

    gt_ir = openmc.ZCylinder(r=0.561, name='Guide tube IR')
    gt_or = openmc.ZCylinder(r=0.602, name='Guide tube OR')

    # Water inside the guide tube (where control rods would go when inserted)
    gt_inner_cell = openmc.Cell(name='GT inner water')
    gt_inner_cell.fill = moderator
    gt_inner_cell.region = -gt_ir
    gt_inner_cell.temperature = 557.0

    # Guide tube wall (Zircaloy-4)
    gt_wall_cell = openmc.Cell(name='GT wall')
    gt_wall_cell.fill = clad
    gt_wall_cell.region = +gt_ir & -gt_or
    gt_wall_cell.temperature = 557.0

    # Water outside the guide tube (fills rest of lattice cell)
    gt_outer_cell = openmc.Cell(name='GT outer water')
    gt_outer_cell.fill = moderator
    gt_outer_cell.region = +gt_or
    gt_outer_cell.temperature = 557.0

    guide_tube_universe = openmc.Universe(
        name='Guide Tube',
        cells=[gt_inner_cell, gt_wall_cell, gt_outer_cell]
    )

    # =====================================================================
    # INSTRUMENT TUBE UNIVERSE
    # =====================================================================
    # The instrument tube at the center of the assembly (position 8,8)
    # houses in-core neutron flux detectors for power monitoring. It has
    # slightly different dimensions than the guide tubes:
    #   Inner radius: 0.559 cm
    #   Outer radius: 0.605 cm
    #   Wall thickness: 0.046 cm = 460 microns
    #
    # The detector itself is not modeled; the tube interior is water-filled.
    # This has negligible neutronic impact since the detector material
    # volume is tiny and has a small cross section compared to the
    # surrounding moderator.

    it_ir = openmc.ZCylinder(r=0.559, name='Instrument tube IR')
    it_or = openmc.ZCylinder(r=0.605, name='Instrument tube OR')

    # Water inside the instrument tube
    it_inner_cell = openmc.Cell(name='IT inner water')
    it_inner_cell.fill = moderator
    it_inner_cell.region = -it_ir
    it_inner_cell.temperature = 557.0

    # Instrument tube wall (Zircaloy-4)
    it_wall_cell = openmc.Cell(name='IT wall')
    it_wall_cell.fill = clad
    it_wall_cell.region = +it_ir & -it_or
    it_wall_cell.temperature = 557.0

    # Water outside the instrument tube
    it_outer_cell = openmc.Cell(name='IT outer water')
    it_outer_cell.fill = moderator
    it_outer_cell.region = +it_or
    it_outer_cell.temperature = 557.0

    instrument_tube_universe = openmc.Universe(
        name='Instrument Tube',
        cells=[it_inner_cell, it_wall_cell, it_outer_cell]
    )

    # =====================================================================
    # 17x17 ASSEMBLY LATTICE
    # =====================================================================
    # Populate the 17x17 lattice with the appropriate universe at each
    # position. The standard Westinghouse layout places:
    #   - 264 fuel pin universes
    #   -  24 guide tube universes
    #   -   1 instrument tube universe (at center)
    #
    # The lattice is centered at the origin. Each cell spans one pin pitch
    # (1.2598 cm x 1.2598 cm).
    #
    # OpenMC RectLattice convention: universes[i][j] where
    #   i = row index (0 = top/highest y, 16 = bottom/lowest y)
    #   j = column index (0 = leftmost x, 16 = rightmost x)

    lattice = openmc.RectLattice(name='17x17 NuScale SMR Assembly')

    # Pin-to-pin spacing
    lattice.pitch = (pin_pitch, pin_pitch)

    # Lower-left corner of the lattice (centers lattice at origin)
    lattice.lower_left = (-17 * pin_pitch / 2.0, -17 * pin_pitch / 2.0)

    # Build the 17x17 universe array
    universes = []
    for row in range(17):
        row_list = []
        for col in range(17):
            if (row, col) == INSTRUMENT_TUBE_POSITION:
                # Center position: instrument tube for flux monitoring
                row_list.append(instrument_tube_universe)
            elif (row, col) in GUIDE_TUBE_POSITIONS:
                # Guide tube position (control rod location, rods withdrawn)
                row_list.append(guide_tube_universe)
            else:
                # Standard fuel pin position (264 of these)
                row_list.append(fuel_pin_universe)
        universes.append(row_list)

    lattice.universes = universes

    # Outer universe: any point outside the 17x17 lattice boundary but
    # inside the assembly boundary is filled with moderator water. This
    # represents the inter-assembly water gap.
    outer_universe = openmc.Universe(name='Outer (water gap)')
    outer_cell = openmc.Cell(name='Inter-assembly water gap', fill=moderator)
    outer_cell.temperature = 557.0
    outer_universe.add_cell(outer_cell)
    lattice.outer = outer_universe

    # =====================================================================
    # ASSEMBLY CELL AND ROOT UNIVERSE
    # =====================================================================
    # The assembly boundary is a rectangular prism at the assembly pitch
    # (21.5036 cm). Reflective BCs on all faces simulate an infinite 2D
    # lattice of identical assemblies.
    #
    # The water gap between the lattice edge (17 * 1.2598 / 2 = 10.7083 cm)
    # and the assembly boundary (21.5036 / 2 = 10.7518 cm) is 0.0435 cm
    # on each side. This thin gap is filled by the lattice's outer universe.

    half_assembly = assembly_pitch / 2.0  # 10.7518 cm

    # Assembly boundary planes with reflective boundary conditions.
    # Reflective BCs are appropriate because this benchmark models a single
    # assembly type in an infinite lattice (no inter-assembly heterogeneity).
    left = openmc.XPlane(x0=-half_assembly, boundary_type='reflective',
                         name='Assembly left boundary')
    right = openmc.XPlane(x0=+half_assembly, boundary_type='reflective',
                          name='Assembly right boundary')
    bottom = openmc.YPlane(y0=-half_assembly, boundary_type='reflective',
                           name='Assembly bottom boundary')
    top = openmc.YPlane(y0=+half_assembly, boundary_type='reflective',
                        name='Assembly top boundary')

    # Axial boundaries: reflective BCs simulate infinite axial extent.
    # This is a 2D benchmark; the axial dimension is a thin slab.
    z_min = openmc.ZPlane(z0=-0.5, boundary_type='reflective',
                          name='Axial bottom (reflective)')
    z_max = openmc.ZPlane(z0=+0.5, boundary_type='reflective',
                          name='Axial top (reflective)')

    # Create the main assembly cell containing the lattice
    assembly_cell = openmc.Cell(name='NuScale SMR Fuel Assembly')
    assembly_cell.fill = lattice
    assembly_cell.region = +left & -right & +bottom & -top & +z_min & -z_max

    # Root universe and geometry
    root_universe = openmc.Universe(name='Root', cells=[assembly_cell])
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    lattice_half = 17 * pin_pitch / 2.0  # 10.7083 cm
    settings = _create_settings(particles, batches, inactive,
                                half_extent=lattice_half)

    # =========================================================================
    # TALLIES
    # =========================================================================
    # Pin-by-pin fission rate tally on a 17x17 mesh. Each mesh element
    # covers one pin position, giving us relative pin powers for validation.
    tallies = openmc.Tallies()

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


def _create_materials():
    """
    Create all materials for the NuScale-like SMR benchmark.

    Returns a tuple of (fuel, gap, clad, moderator) materials.

    The critical difference from a conventional PWR model is the moderator:
    NuScale operates with ZERO soluble boron. Reactivity control is achieved
    entirely through burnable absorbers (Gd2O3 in some assembly types) and
    mechanical control rods. This significantly simplifies the plant's
    chemical and volume control system.
    """

    # ---- UO2 Fuel (3.1 wt% U-235 enriched) ----
    # Density: 10.257 g/cc (standard UO2 pellet density ~95% theoretical)
    # Enrichment: 3.1 w/o U-235 (the highest enrichment zone in the
    #   NuScale-like SMR benchmark; other zones use 2.0% and 2.6%)
    #
    # Number densities are specified directly in atoms/barn-cm for
    # consistency with the benchmark specification. These include trace
    # amounts of U-234 and U-236 that arise from the isotopic enrichment
    # process (centrifuge cascade tails).
    fuel = openmc.Material(name='UO2 Fuel 3.1%')
    fuel.set_density('sum')  # density determined by sum of number densities
    fuel.add_nuclide('U234', 6.11864e-06)   # trace from enrichment process
    fuel.add_nuclide('U235', 7.18132e-04)   # fissile isotope (3.1 w/o)
    fuel.add_nuclide('U236', 3.29861e-06)   # trace from enrichment process
    fuel.add_nuclide('U238', 2.21546e-02)   # fertile isotope (bulk of uranium)
    fuel.add_nuclide('O16',  4.57642e-02)   # oxygen in UO2 (two per molecule)

    # ---- Helium Gap ----
    # A thin gas-filled gap exists between the UO2 fuel pellet outer
    # surface and the inner surface of the Zircaloy-4 cladding tube.
    # Purpose: accommodate fuel pellet swelling (thermal expansion +
    # fission product accumulation) and fission gas release during
    # burnup. At beginning of life, it is filled with helium gas at
    # low pressure (~1-3 atm).
    gap = openmc.Material(name='Helium Gap')
    gap.set_density('g/cc', 0.001598)  # helium at ~2 atm, ~900K
    gap.add_nuclide('He4', 1.0)

    # ---- Zircaloy-4 Cladding ----
    # Zircaloy-4 (Zr-4) is the standard cladding material for PWR fuel
    # rods. It was selected for its combination of:
    #   - Low thermal neutron absorption cross section (Zr is nearly
    #     transparent to neutrons, minimizing parasitic absorption)
    #   - Good corrosion resistance in high-temperature water
    #   - Adequate mechanical strength at operating temperatures
    #   - Compatibility with UO2 fuel
    #
    # Composition (weight percent):
    #   Zr: 98.2% (balance)
    #   Sn: 1.5% (improves mechanical strength and corrosion resistance)
    #   Fe: 0.2% (improves corrosion resistance)
    #   Cr: 0.1% (improves corrosion resistance)
    #
    # Zircaloy-4 is also used for the guide tube and instrument tube walls.
    clad = openmc.Material(name='Zircaloy-4')
    clad.set_density('g/cc', 6.56)
    clad.add_element('Zr', 0.982, 'wo')  # zirconium (balance)
    clad.add_element('Sn', 0.015, 'wo')  # tin
    clad.add_element('Fe', 0.002, 'wo')  # iron
    clad.add_element('Cr', 0.001, 'wo')  # chromium

    # ---- Unborated Light Water Moderator ----
    # THIS IS THE KEY DIFFERENCE FROM CONVENTIONAL PWR MODELS.
    #
    # In a standard PWR (like the VERA benchmark), the moderator contains
    # dissolved boric acid (H3BO3) at concentrations of 500-1500 ppm boron.
    # The B-10 isotope (19.9% natural abundance) has a very large thermal
    # neutron absorption cross section (~3840 barns at 0.025 eV), making
    # soluble boron an effective means of reactivity control.
    #
    # The NuScale SMR operates with ZERO soluble boron. Instead, reactivity
    # control relies on:
    #   1. Burnable absorbers: Gd2O3 (gadolinium oxide) mixed into some
    #      fuel pins. Gd-155 and Gd-157 have enormous thermal absorption
    #      cross sections (~61,000 and ~254,000 barns) that suppress
    #      excess reactivity at beginning of life and gradually deplete.
    #   2. Control rods: Ag-In-Cd or B4C absorber rods inserted into
    #      the guide tubes for fine reactivity control and shutdown.
    #
    # Benefits of boron-free operation:
    #   - Eliminates risk of boron dilution accidents
    #   - Simplifies chemical volume control system (CVCS)
    #   - Eliminates boron precipitation concerns during LOCA
    #   - More negative moderator temperature coefficient (safety benefit)
    #
    # Moderator conditions:
    #   Temperature: 557 K (~284 degC)
    #   Pressure: 12.76 MPa
    #   Density: 0.7527 g/cc (subcooled water at these conditions)
    moderator = openmc.Material(name='Unborated Water')
    moderator.set_density('g/cc', 0.7527)
    moderator.add_nuclide('H1', 2.0)   # two hydrogen atoms per water molecule
    moderator.add_nuclide('O16', 1.0)  # one oxygen atom per water molecule
    # Add S(alpha,beta) thermal scattering data for hydrogen bound in water.
    # This is CRITICAL for correctly modeling neutron thermalization. Without
    # it, the free-gas scattering kernel would be used, which does not account
    # for the molecular binding effects in H2O. This would give an incorrect
    # thermal neutron spectrum and a significantly biased k-eff (typically
    # several hundred pcm or more).
    moderator.add_s_alpha_beta('c_H_in_H2O')

    return fuel, gap, clad, moderator


def _create_settings(particles, batches, inactive, half_extent):
    """
    Create OpenMC settings for the eigenvalue calculation.

    Parameters
    ----------
    particles : int
        Neutron histories per batch.
    batches : int
        Total batches (active + inactive).
    inactive : int
        Inactive batches for fission source convergence.
    half_extent : float
        Half-width of the source box in x and y [cm].

    Returns
    -------
    openmc.Settings
        Configured settings object.
    """
    settings = openmc.Settings()

    # Eigenvalue mode: solve the neutron transport eigenvalue problem
    # to find k-effective, the multiplication factor.
    settings.run_mode = 'eigenvalue'

    # Number of neutron histories simulated per batch.
    settings.particles = particles

    # Total number of batches. The first 'inactive' batches are discarded
    # to allow the fission source distribution to converge to the
    # fundamental mode. The remaining (batches - inactive) active batches
    # are used for tallying.
    settings.batches = batches
    settings.inactive = inactive

    # Temperature handling: use nearest available temperature in the
    # cross section library. The library may not have data at exactly
    # 900K or 557K, so we allow OpenMC to interpolate or use the
    # nearest available data set with a generous tolerance.
    settings.temperature = {'method': 'nearest', 'tolerance': 1000}

    # Initial fission source distribution: uniform box covering the
    # fuel region. OpenMC samples points uniformly in the box and
    # rejects any that don't land in fissile material.
    settings.source = [openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-half_extent, -half_extent, -0.5),
            upper_right=(+half_extent, +half_extent, +0.5)
        )
    )]

    return settings


def main():
    """Parse command-line arguments, build the model, and export/run."""
    parser = argparse.ArgumentParser(
        description='NuScale-like SMR Benchmark: Single Fuel Assembly (OpenMC)'
    )
    parser.add_argument('--particles', type=int, default=50000,
                        help='Neutron histories per batch (default: 50000)')
    parser.add_argument('--batches', type=int, default=200,
                        help='Total active batches (default: 200)')
    parser.add_argument('--inactive', type=int, default=50,
                        help='Inactive batches for source convergence (default: 50)')
    parser.add_argument('--model', type=str, default='assembly',
                        choices=['pin-cell', 'assembly'],
                        help='Model type: pin-cell or assembly (default: assembly)')
    parser.add_argument('--run', action='store_true',
                        help='Run OpenMC after exporting XML files')
    args = parser.parse_args()

    # Build the requested model
    if args.model == 'pin-cell':
        print("Building NuScale-like SMR pin cell model...")
        model = build_pin_cell(
            particles=args.particles,
            batches=args.batches,
            inactive=args.inactive
        )
        print("  Single fuel pin with reflective BCs (infinite pin lattice)")
        print("  Key feature: UNBORATED moderator (NuScale boron-free design)")
    else:
        print("Building NuScale-like SMR 17x17 assembly model...")
        model = build_assembly(
            particles=args.particles,
            batches=args.batches,
            inactive=args.inactive
        )
        # Print assembly layout summary
        n_fuel = 264
        n_gt = len(GUIDE_TUBE_POSITIONS)
        n_it = 1
        print(f"  Layout: {n_fuel} fuel pins, {n_gt} guide tubes, "
              f"{n_it} instrument tube = {n_fuel + n_gt + n_it} positions")
        print(f"  Pin pitch: 1.2598 cm, Assembly pitch: 21.5036 cm")
        print(f"  Key feature: UNBORATED moderator (NuScale boron-free design)")

    # Export XML input files
    model.export_to_xml()
    print("Exported OpenMC XML input files.")

    # Optionally run the simulation
    if args.run:
        print(f"\nRunning OpenMC: {args.particles} particles/batch, "
              f"{args.batches} active + {args.inactive} inactive batches")
        model.run()


if __name__ == '__main__':
    main()
