#!/usr/bin/env python3
"""
BEAVRS Benchmark - Single 17x17 Fuel Assembly at Hot Zero Power
================================================================

This script builds an OpenMC model for a single fuel assembly from the BEAVRS
(Benchmark for Evaluation And Validation of Reactor Simulations) benchmark.

The BEAVRS benchmark is described in:

    N. Horelik, B. Herman, B. Forget, and K. Smith,
    "Benchmark for Evaluation and Validation of Reactor Simulations (BEAVRS),
     v1.0.1," Proc. Int. Conf. Mathematics and Computational Methods Applied
     to Nucl. Sci. & Eng., Sun Valley, Idaho, 2013.

    Full specification: MIT Computational Reactor Physics Group
    https://github.com/mit-crpg/BEAVRS

BEAVRS is based on REAL MEASURED DATA from an operating Westinghouse 4-loop
pressurized water reactor rated at 3411 MWth. Unlike many computational
benchmarks that compare code-to-code, BEAVRS provides actual plant
measurements including:
  - Hot zero power (HZP) critical boron concentrations
  - Boron letdown curves during cycle 1 and cycle 2
  - 3D in-core flux maps from 58 instrumented assemblies
  - Control rod bank worths
  - Isothermal temperature coefficients

This makes BEAVRS uniquely valuable for true validation (not just
verification) of reactor physics codes against reality.

Single Assembly Model (this script)
-------------------------------------
We model a single 17x17 fuel assembly with reflective boundary conditions on
all faces, simulating an infinite lattice of identical assemblies. This is
a standard first step in BEAVRS analysis and isolates the assembly-level
neutronics from full-core effects (leakage, inter-assembly heterogeneity).

The assembly contains three types of pin locations:
  - 264 fuel pins: UO2 pellet + He gap + Zircaloy-4 clad + moderator
  -  24 guide tubes: water-filled Zirc-4 tubes (control rods withdrawn)
  -   1 instrument tube: center position (8,8), houses in-core detectors

Together: 264 + 24 + 1 = 289 = 17 x 17 positions.

Hot Zero Power (HZP) conditions:
  - Fuel temperature:      600 K (BEAVRS HZP specification)
  - Moderator temperature: 565 K
  - Moderator density:     0.740 g/cc (at ~565 K)
  - Soluble boron:         975 ppm (measured critical boron at HZP)
  - No control rods inserted (all rods out, ARO)

Key differences from the VERA benchmark (CASL):
  - Fuel pellet radius:  0.39218 cm (BEAVRS) vs 0.4096 cm (VERA)
  - Clad inner radius:   0.40005 cm (BEAVRS) vs 0.418  cm (VERA)
  - Clad outer radius:   0.45720 cm (BEAVRS) vs 0.475  cm (VERA)
  - Pin pitch:           1.25984 cm (BEAVRS) vs 1.26   cm (VERA)
  - Assembly pitch:      21.50364 cm (BEAVRS) vs 21.50 cm (VERA)
  - Critical boron:      975 ppm (BEAVRS, measured) vs 1300 ppm (VERA)
  - Reference data:      real plant measurements vs computed reference

Reference Results:
  - HZP measured critical boron concentration: ~975 ppm
  - At 975 ppm, k-eff should be ~1.000 for full core (our infinite lattice
    result will differ due to absence of leakage and assembly heterogeneity)

Usage:
    python model.py [--particles N] [--batches N] [--inactive N]
                    [--enrichment {1.6,2.4,3.1}] [--boron-ppm N]
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
# The guide tubes are symmetrically placed in the assembly following the
# standard Westinghouse pattern. During power operation, 24 of these tubes
# house control rod fingers (in rodded assemblies) or are simply water-filled
# (in unrodded assemblies). The center position (8,8) holds the instrument
# thimble for in-core flux measurement.
#
# These positions are identical to those used in the VERA benchmark (CASL),
# as both are based on the standard Westinghouse 17x17 design.

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
# It houses in-core neutron flux detectors used for core monitoring.
# In the BEAVRS plant, 58 of the 193 assemblies are instrumented.
INSTRUMENT_TUBE_POSITION = (8, 8)


# =============================================================================
# UO2 FUEL COMPOSITIONS FOR DIFFERENT ENRICHMENTS
# =============================================================================
# BEAVRS Cycle 1 uses three fuel enrichments: 1.6, 2.4, and 3.1 wt% U-235.
# The default is 3.1 wt%, which is the most common enrichment in the core.
# Number densities are computed from standard UO2 relationships at the given
# enrichment and theoretical density of 10.257 g/cc.

# Precomputed atom densities [atoms/barn-cm] for UO2 at 10.257 g/cc
# at each enrichment level. These use standard molecular weight relationships
# for UO2 and natural oxygen-16.
FUEL_COMPOSITIONS = {
    1.6: {
        # 1.6 wt% U-235 enrichment (lowest enrichment in BEAVRS Cycle 1)
        # Used in central region of the core where flux/power is highest
        'U235': 3.7025e-04,
        'U238': 2.2407e-02,
        'O16':  4.5759e-02,
    },
    2.4: {
        # 2.4 wt% U-235 enrichment (intermediate enrichment)
        # Used in the intermediate radial region of the core
        'U235': 5.5616e-04,
        'U238': 2.2263e-02,
        'O16':  4.5759e-02,
    },
    3.1: {
        # 3.1 wt% U-235 enrichment (highest enrichment in BEAVRS Cycle 1)
        # Used in the peripheral region of the core to flatten the power
        'U235': 7.18132e-04,
        'U238': 2.21546e-02,
        'O16':  4.57642e-02,
    },
}


def build_model(particles=30000, batches=200, inactive=50, enrichment=3.1,
                boron_ppm=975, small_tallies=False):
    """
    Construct the OpenMC model for a BEAVRS single fuel assembly.

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch.
    batches : int
        Total number of batches (active + inactive).
    inactive : int
        Number of inactive (discarded) batches for source convergence.
    enrichment : float
        Fuel enrichment in weight percent U-235. Must be 1.6, 2.4, or 3.1.
    boron_ppm : float
        Soluble boron concentration in ppm. Default 975 ppm corresponds to
        the BEAVRS HZP measured critical boron concentration.
    small_tallies : bool
        If True, skip the pin power mesh tally (faster runtime).

    Returns
    -------
    openmc.Model
        The complete OpenMC model ready for export or execution.
    """

    # Validate enrichment
    if enrichment not in FUEL_COMPOSITIONS:
        raise ValueError(
            f"Enrichment {enrichment} not available. "
            f"Choose from: {list(FUEL_COMPOSITIONS.keys())}"
        )

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # Material compositions for the BEAVRS benchmark. The fuel uses explicit
    # number densities computed from the UO2 formula at 10.257 g/cc theoretical
    # density. The moderator uses the BEAVRS HZP specifications.

    # ---- UO2 Fuel ----
    # Density: 10.257 g/cc (95.5% of theoretical UO2 density of 10.97 g/cc)
    # This is the standard assumed pellet density for Westinghouse PWR fuel.
    fuel_comp = FUEL_COMPOSITIONS[enrichment]
    fuel = openmc.Material(name=f'UO2 Fuel {enrichment}%')
    fuel.set_density('sum')  # density determined by sum of number densities
    for nuclide, density in fuel_comp.items():
        fuel.add_nuclide(nuclide, density)

    # ---- Helium Gap ----
    # A thin gas-filled gap exists between the fuel pellet outer surface and
    # the cladding inner surface. At beginning of life, this gap is filled
    # with helium backfill gas at approximately 1 atmosphere. The gap serves
    # multiple purposes:
    #   1. Accommodates fuel pellet swelling during irradiation
    #   2. Provides space for fission gas release
    #   3. The helium provides better thermal conductivity than a vacuum
    # The gap thickness in BEAVRS is 0.40005 - 0.39218 = 0.00787 cm = 78.7 um
    gap = openmc.Material(name='Helium Gap')
    gap.set_density('g/cc', 0.001598)  # helium at ~1 atm, ~600K
    gap.add_nuclide('He4', 1.0)

    # ---- Zircaloy-4 Cladding ----
    # Zircaloy-4 is the standard cladding and structural tube material in
    # Westinghouse PWR fuel assemblies. Its key properties:
    #   - Very low thermal neutron absorption cross section (~0.18 barns)
    #   - Good corrosion resistance in high-temperature water
    #   - Adequate mechanical strength at operating temperatures
    #   - Maintains integrity under irradiation
    # Composition (by weight): ~98.2% Zr, 1.5% Sn, 0.2% Fe, 0.1% Cr
    # Density: 6.56 g/cc
    clad = openmc.Material(name='Zircaloy-4')
    clad.set_density('g/cc', 6.56)
    clad.add_element('Zr', 0.982, 'wo')  # zirconium (balance)
    clad.add_element('Sn', 0.015, 'wo')  # tin - solid solution strengthener
    clad.add_element('Fe', 0.002, 'wo')  # iron - precipitate former
    clad.add_element('Cr', 0.001, 'wo')  # chromium - corrosion resistance

    # ---- Borated Water Moderator ----
    # The moderator is light water with dissolved boric acid (H3BO3). Soluble
    # boron is the primary means of long-term reactivity control in a PWR.
    # Boron-10 (natural abundance ~19.9%) has a very large thermal neutron
    # absorption cross section (~3840 barns at 0.025 eV), making it an
    # effective neutron poison.
    #
    # BEAVRS HZP conditions:
    #   - Temperature: 565 K (~292 C, ~557 F)
    #   - Density: 0.740 g/cc (subcooled liquid at HZP pressure ~155 bar)
    #   - Boron: 975 ppm (measured critical concentration at HZP)
    #
    # NOTE: If the cross section library only has room-temperature (294K) data,
    # there will be a temperature bias in k-eff (Doppler effect). Using
    # 0.740 g/cc density with 294K cross sections mixes two temperature
    # regimes, but this is standard practice when HZP-temperature libraries
    # are unavailable.
    #
    # We compute the number densities from the density and boron concentration.
    # Water molecular weight: 18.015 g/mol
    # Avogadro's number: 6.02214e23 /mol
    # 1 barn-cm = 1e-24 cm^2 * cm = 1e-24 cm^3
    moderator_density = 0.740  # g/cc at 565 K

    # Compute water number density
    # n_H2O = rho * N_A / M_H2O (in atoms/barn-cm, need factor of 1e-24)
    N_A = 6.02214076e23  # Avogadro's number
    M_H2O = 18.015       # g/mol for H2O

    # Number density of water molecules [molecules/barn-cm]
    n_water = moderator_density * N_A / M_H2O * 1e-24  # molecules/barn-cm

    # Hydrogen and oxygen from water
    n_H = 2.0 * n_water   # 2 H atoms per water molecule
    n_O = 1.0 * n_water    # 1 O atom per water molecule

    # Boron number density from ppm concentration
    # ppm is mass of boron per mass of solution * 1e6
    # n_B = (boron_ppm * 1e-6) * rho_mod * N_A / M_B * 1e-24
    M_B = 10.811  # g/mol (natural boron atomic weight)
    n_B_total = (boron_ppm * 1e-6) * moderator_density * N_A / M_B * 1e-24

    # Natural boron isotopic fractions: 19.9% B-10, 80.1% B-11
    n_B10 = 0.199 * n_B_total
    n_B11 = 0.801 * n_B_total

    moderator = openmc.Material(name='Borated Water')
    moderator.set_density('sum')
    moderator.add_nuclide('H1', n_H)
    moderator.add_nuclide('O16', n_O)
    moderator.add_nuclide('B10', n_B10)
    moderator.add_nuclide('B11', n_B11)
    # S(alpha,beta) thermal scattering law for hydrogen bound in water.
    # This is CRITICAL for correct neutron thermalization. Without it,
    # OpenMC would use the free-gas scattering kernel for hydrogen, which
    # does not account for the molecular binding effects in water. The
    # bound kernel properly treats:
    #   - Translational motion of the water molecule
    #   - Rotational modes of H around the O-H bond
    #   - Vibrational modes of the O-H stretch
    # Using the free-gas kernel would significantly bias k-eff (~500+ pcm).
    moderator.add_s_alpha_beta('c_H_in_H2O')

    # Collect all materials
    materials = openmc.Materials([fuel, gap, clad, moderator])

    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # The BEAVRS geometry uses slightly different dimensions from VERA.
    # All dimensions are from the BEAVRS specification document.

    # --- Pin pitch ---
    # Center-to-center distance between adjacent pin positions.
    # BEAVRS: 1.25984 cm (vs VERA: 1.26 cm - a small but real difference)
    pin_pitch = 1.25984  # cm

    # --- Assembly pitch ---
    # Center-to-center distance between adjacent fuel assemblies in the core.
    # This is larger than 17 * pin_pitch to account for the inter-assembly
    # water gap. BEAVRS: 21.50364 cm
    # Lattice width: 17 * 1.25984 = 21.41728 cm
    # Water gap on each side: (21.50364 - 21.41728) / 2 = 0.04318 cm
    assembly_pitch = 21.50364  # cm

    # =====================================================================
    # FUEL PIN UNIVERSE
    # =====================================================================
    # The fuel pin is a set of concentric cylinders representing (from
    # center outward): UO2 fuel pellet, helium gap, Zircaloy-4 cladding,
    # and surrounding moderator water.
    #
    # BEAVRS fuel pin dimensions:
    #   Fuel pellet OR:  0.39218 cm (smaller than VERA's 0.4096 cm)
    #   Clad IR (gap):   0.40005 cm
    #   Clad OR:         0.45720 cm
    #   Gap thickness:   0.40005 - 0.39218 = 0.00787 cm = 78.7 um
    #   Clad thickness:  0.45720 - 0.40005 = 0.05715 cm = 571.5 um
    #
    # Cross-section (not to scale):
    #
    #     Moderator (borated water)
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

    # Fuel pellet outer radius
    fuel_or = openmc.ZCylinder(r=0.39218, name='Fuel pellet OR')
    # Cladding inner radius (defines outer boundary of helium gap)
    clad_ir = openmc.ZCylinder(r=0.40005, name='Fuel clad IR')
    # Cladding outer radius
    clad_or = openmc.ZCylinder(r=0.45720, name='Fuel clad OR')

    # Fuel pellet cell: inside the fuel outer radius
    fuel_cell = openmc.Cell(name='Fuel pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or
    fuel_cell.temperature = 600.0  # BEAVRS HZP fuel temperature [K]

    # Helium gap cell: annular region between fuel pellet and clad
    gap_cell = openmc.Cell(name='Helium gap')
    gap_cell.fill = gap
    gap_cell.region = +fuel_or & -clad_ir
    gap_cell.temperature = 600.0

    # Cladding cell: annular region of Zircaloy-4
    clad_cell = openmc.Cell(name='Fuel cladding')
    clad_cell.fill = clad
    clad_cell.region = +clad_ir & -clad_or
    clad_cell.temperature = 565.0  # clad temperature closer to coolant

    # Moderator cell: everything outside the cladding within this universe.
    # When placed in a lattice, OpenMC clips each universe to its lattice
    # cell boundaries automatically.
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
    # Guide tubes are larger-diameter Zircaloy-4 tubes that house control
    # rod fingers during operation. In this HZP all-rods-out configuration,
    # the guide tubes are simply filled with moderator water.
    #
    # BEAVRS guide tube dimensions:
    #   Inner radius: 0.56134 cm
    #   Outer radius: 0.60198 cm
    #   Wall thickness: 0.60198 - 0.56134 = 0.04064 cm = 406.4 um
    #
    # These are significantly larger than the fuel pin cladding OR
    # (0.45720 cm) to allow insertion of control rod absorber pins.

    gt_ir = openmc.ZCylinder(r=0.56134, name='Guide tube IR')
    gt_or = openmc.ZCylinder(r=0.60198, name='Guide tube OR')

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

    guide_tube_universe = openmc.Universe(
        name='Guide Tube',
        cells=[gt_inner_cell, gt_wall_cell, gt_outer_cell]
    )

    # =====================================================================
    # INSTRUMENT TUBE UNIVERSE
    # =====================================================================
    # The instrument tube at position (8,8) houses in-core neutron flux
    # detectors used for core power monitoring. In the BEAVRS plant, 58
    # of 193 fuel assemblies are instrumented with movable in-core
    # detector systems.
    #
    # BEAVRS instrument tube dimensions (slightly different from guide tubes):
    #   Inner radius: 0.55880 cm
    #   Outer radius: 0.60452 cm
    #   Wall thickness: 0.60452 - 0.55880 = 0.04572 cm = 457.2 um
    #
    # The instrument tube has a slightly smaller inner radius and slightly
    # larger outer radius than the guide tube, giving a thicker wall.
    # The detector itself is not modeled (negligible neutronic effect).

    it_ir = openmc.ZCylinder(r=0.55880, name='Instrument tube IR')
    it_or = openmc.ZCylinder(r=0.60452, name='Instrument tube OR')

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

    instrument_tube_universe = openmc.Universe(
        name='Instrument Tube',
        cells=[it_inner_cell, it_wall_cell, it_outer_cell]
    )

    # =====================================================================
    # 17x17 ASSEMBLY LATTICE
    # =====================================================================
    # Build the 17x17 lattice by populating a 2D array of universes.
    # Most positions get the fuel pin universe; the 24 guide tube positions
    # and 1 instrument tube position are overridden.
    #
    # The lattice is centered at the origin. Each element spans one
    # pin pitch (1.25984 cm x 1.25984 cm).

    lattice = openmc.RectLattice(name='17x17 BEAVRS Fuel Assembly Lattice')

    # Lattice pitch = pin-to-pin spacing
    lattice.pitch = (pin_pitch, pin_pitch)

    # Center the lattice at the origin
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

    # Outer universe: any point outside the lattice boundaries gets this.
    # The inter-assembly gap is filled with moderator water.
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
    # (21.50364 cm) to properly represent the inter-assembly water gap.
    #
    # Lattice extent: 17 * 1.25984 = 21.41728 cm
    # Assembly pitch: 21.50364 cm
    # Water gap each side: (21.50364 - 21.41728) / 2 = 0.04318 cm
    # This thin water gap between assemblies is automatically filled by
    # the lattice's outer universe.

    half_assembly = assembly_pitch / 2.0  # 10.75182 cm

    # Assembly boundary planes with reflective BCs.
    # Reflective boundaries simulate an infinite 2D lattice of identical
    # assemblies - appropriate for a single-assembly eigenvalue problem.
    left = openmc.XPlane(x0=-half_assembly, boundary_type='reflective',
                         name='Assembly left boundary')
    right = openmc.XPlane(x0=+half_assembly, boundary_type='reflective',
                          name='Assembly right boundary')
    bottom = openmc.YPlane(y0=-half_assembly, boundary_type='reflective',
                           name='Assembly bottom boundary')
    top = openmc.YPlane(y0=+half_assembly, boundary_type='reflective',
                        name='Assembly top boundary')

    # Axial boundaries - 2D problem approximated with a thin slab.
    # Reflective BCs in z simulate infinite axial extent (no axial leakage).
    z_min = openmc.ZPlane(z0=-0.5, boundary_type='reflective',
                          name='Axial bottom (reflective)')
    z_max = openmc.ZPlane(z0=+0.5, boundary_type='reflective',
                          name='Axial top (reflective)')

    # Create the main assembly cell containing the lattice
    assembly_cell = openmc.Cell(name='BEAVRS Fuel Assembly')
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

    # Eigenvalue mode: solve for the fundamental mode multiplication factor.
    settings.run_mode = 'eigenvalue'

    # Number of neutron histories per batch.
    # 30,000 is adequate for a fuel assembly (lower than VERA default since
    # BEAVRS assembly is very similar in complexity).
    settings.particles = particles

    # Total number of batches (active + inactive).
    settings.batches = batches

    # Inactive batches for fission source convergence.
    # The fission source must equilibrate across all 264 fuel pins before
    # we start accumulating statistics. 50 inactive batches is typically
    # sufficient for a single assembly.
    settings.inactive = inactive

    # Temperature settings: use nearest available temperature data.
    # The cross section library may not have data at exactly 565K or 600K,
    # so we tell OpenMC to use the nearest available temperature.
    settings.temperature = {'method': 'nearest', 'tolerance': 1000}

    # Initial source distribution: uniform box spanning the fuel region.
    # We cover the entire lattice so that all fuel pins receive initial
    # source particles. OpenMC will reject samples in non-fissile regions.
    lattice_half = 17 * pin_pitch / 2.0  # 10.70864 cm
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
        # Pin power distributions are a key validation quantity for
        # assembly calculations. In the BEAVRS benchmark, measured in-core
        # detector signals provide real validation data for flux
        # distributions (though at the full-core level, not single assembly).

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
        description='BEAVRS Benchmark: 17x17 PWR Fuel Assembly (OpenMC)'
    )
    parser.add_argument('--particles', type=int, default=30000,
                        help='Neutron histories per batch (default: 30000)')
    parser.add_argument('--batches', type=int, default=200,
                        help='Total active batches (default: 200)')
    parser.add_argument('--inactive', type=int, default=50,
                        help='Inactive batches for source convergence (default: 50)')
    parser.add_argument('--enrichment', type=float, default=3.1,
                        choices=[1.6, 2.4, 3.1],
                        help='Fuel enrichment in wt%% U-235 (default: 3.1)')
    parser.add_argument('--boron-ppm', type=float, default=975,
                        help='Soluble boron concentration in ppm (default: 975)')
    parser.add_argument('--small-tallies', action='store_true',
                        help='Skip pin power mesh tally (faster runtime)')
    parser.add_argument('--run', action='store_true',
                        help='Run OpenMC after exporting XML files')
    args = parser.parse_args()

    # Build the model with specified parameters
    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
        enrichment=args.enrichment,
        boron_ppm=args.boron_ppm,
        small_tallies=args.small_tallies
    )

    # Export XML input files
    model.export_to_xml()
    print("Exported OpenMC XML input files for BEAVRS single assembly.")
    print(f"  Fuel enrichment: {args.enrichment} wt% U-235")
    print(f"  Boron concentration: {args.boron_ppm} ppm")

    # Print assembly layout summary
    n_fuel = 264
    n_gt = len(GUIDE_TUBE_POSITIONS)
    n_it = 1
    print(f"  Assembly layout: {n_fuel} fuel pins, {n_gt} guide tubes, "
          f"{n_it} instrument tube = {n_fuel + n_gt + n_it} positions")

    # Optionally run the simulation
    if args.run:
        print(f"\nRunning OpenMC: {args.particles} particles/batch, "
              f"{args.batches - args.inactive} active + "
              f"{args.inactive} inactive batches")
        model.run()


if __name__ == '__main__':
    main()
