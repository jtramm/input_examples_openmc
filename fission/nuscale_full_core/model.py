#!/usr/bin/env python3
"""
McSAFER NuScale-like SMR Full-Core Benchmark (37 Fuel Assemblies)
==================================================================

This script builds an OpenMC model of the NuScale-like SMR full core as
defined in the Euratom McSAFER project benchmark:

    "Neutronics Benchmark of the NuScale-like SMR in the McSAFER Project"
    J. Bousquet et al., J. Nucl. Eng. 2025, 6(4), 44
    https://www.mdpi.com/2673-4362/6/4/44

This is the FULL CORE model with all 37 fuel assemblies, as opposed to
the single-assembly infinite-lattice model in ../nuscale_smr/.

NuScale Power Module Overview
------------------------------
The NuScale SMR is a 77 MWe (160 MWth per module in the original design,
later uprated to 250 MWth) integral pressurized water reactor. The entire
primary system is housed in a single integrated vessel submerged in a
below-grade pool.

Core Layout
-----------
The 37 fuel assemblies are arranged in a 7x7 grid with the four corners
removed, giving a roughly cylindrical core cross-section:

    Row 0:   _    _    3    3    3    _    _
    Row 1:   _    3    2    2    2    3    _
    Row 2:   3    2    1    1    1    2    3
    Row 3:   3    2    1    1    1    2    3
    Row 4:   3    2    1    1    1    2    3
    Row 5:   _    3    2    2    2    3    _
    Row 6:   _    _    3    3    3    _    _

    1 = 1.6 wt% U-235 (center region, 9 assemblies)
    2 = 2.4 wt% U-235 (middle ring, 12 assemblies)
    3 = 3.1 wt% U-235 (outer ring, 16 assemblies)
    _ = water (no assembly)

This gives exactly 37 fuel assemblies, making this the smallest full-core
reactor benchmark -- ideal for memory-constrained systems.

Assembly Types (simplified -- no burnable absorbers)
-----------------------------------------------------
All three assembly types use the same 17x17 Westinghouse layout with
264 fuel pins, 24 guide tubes, and 1 instrument tube. They differ only
in fuel enrichment:
  - Type 1: 1.6 wt% U-235 (lowest enrichment, center)
  - Type 2: 2.4 wt% U-235 (intermediate enrichment, middle ring)
  - Type 3: 3.1 wt% U-235 (highest enrichment, outer ring)

The original McSAFER benchmark includes 7 assembly types with varying
Gd2O3 burnable absorber pin counts. For simplicity, this model omits
gadolinium and uses enrichment-only differentiation.

Operating Conditions (Hot Zero Power)
--------------------------------------
  - Fuel temperature:      543 K (HZP -- same as moderator)
  - Moderator temperature: 543 K (~270 degC)
  - Moderator density:     0.725 g/cc (subcooled water at NuScale conditions)
  - Soluble boron:         0 ppm (BORON-FREE design)
  - Active fuel height:    200 cm
  - Core power:            160 MWth

Geometry
--------
  - Pin pitch:             1.26 cm
  - Assembly pitch:        21.50 cm
  - Fuel pellet OR:        0.4096 cm
  - Fuel clad IR:          0.418 cm
  - Fuel clad OR:          0.475 cm
  - Guide tube IR:         0.561 cm
  - Guide tube OR:         0.602 cm
  - Instrument tube IR:    0.559 cm
  - Instrument tube OR:    0.605 cm
  - Core barrel IR:        ~100 cm
  - Core barrel OR:        ~105 cm

Reference Results
------------------
  HZP, ARO, fresh fuel: k-eff ~ 1.01-1.03 (from Serpent reference)
  SCALE/KENO-VI agrees within 44 pcm

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# GUIDE TUBE AND INSTRUMENT TUBE POSITIONS (17x17 Westinghouse layout)
# =============================================================================
# These are the same for all assembly types. Positions given as (row, col)
# with 0-based indexing, where (0,0) is the top-left corner.

GUIDE_TUBE_POSITIONS = [
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5),          (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
]

INSTRUMENT_TUBE_POSITION = (8, 8)


# =============================================================================
# CORE LOADING PATTERN
# =============================================================================
# 7x7 grid with corners removed. Values:
#   0 = empty (water)
#   1 = Type 1 assembly (1.6 wt% U-235)
#   2 = Type 2 assembly (2.4 wt% U-235)
#   3 = Type 3 assembly (3.1 wt% U-235)
#
# This loading pattern places the lowest-enrichment fuel in the center
# (where the neutron flux is highest) and the highest-enrichment fuel on
# the periphery (where leakage reduces the flux). This "low-leakage"
# loading pattern flattens the radial power distribution.

CORE_MAP = [
    [0, 0, 3, 3, 3, 0, 0],
    [0, 3, 2, 2, 2, 3, 0],
    [3, 2, 1, 1, 1, 2, 3],
    [3, 2, 1, 1, 1, 2, 3],
    [3, 2, 1, 1, 1, 2, 3],
    [0, 3, 2, 2, 2, 3, 0],
    [0, 0, 3, 3, 3, 0, 0],
]

# Verify: count assemblies
_n_assemblies = sum(1 for row in CORE_MAP for val in row if val > 0)
assert _n_assemblies == 37, f"Expected 37 assemblies, got {_n_assemblies}"


# =============================================================================
# ENRICHMENT SPECIFICATIONS FOR EACH ASSEMBLY TYPE
# =============================================================================
# Number densities for UO2 at 10.257 g/cc at different enrichments.
# These are computed from the enrichment and UO2 molecular weight.
#
# At 10.257 g/cc, the total uranium number density is approximately
# 2.2373e-2 atoms/barn-cm, and the oxygen density is ~4.5764e-2.
# The U-235 fraction varies by enrichment; U-234 and U-236 are trace
# isotopes from the enrichment process (scaled proportionally).

FUEL_COMPOSITIONS = {
    # Type 1: 1.6 wt% U-235
    1: {
        'name': 'UO2 Fuel 1.6%',
        'nuclides': {
            'U234': 3.06854e-06,
            'U235': 3.70512e-04,
            'U236': 1.70327e-06,
            'U238': 2.24025e-02,
            'O16':  4.57642e-02,
        }
    },
    # Type 2: 2.4 wt% U-235
    2: {
        'name': 'UO2 Fuel 2.4%',
        'nuclides': {
            'U234': 4.63210e-06,
            'U235': 5.55730e-04,
            'U236': 2.55442e-06,
            'U238': 2.22173e-02,
            'O16':  4.57642e-02,
        }
    },
    # Type 3: 3.1 wt% U-235
    3: {
        'name': 'UO2 Fuel 3.1%',
        'nuclides': {
            'U234': 6.11864e-06,
            'U235': 7.18132e-04,
            'U236': 3.29861e-06,
            'U238': 2.21546e-02,
            'O16':  4.57642e-02,
        }
    },
}


def build_model(particles=20000, batches=200, inactive=50):
    """
    Construct the NuScale-like SMR full-core model with 37 fuel assemblies.

    The model builds three fuel enrichment variants, each as a full 17x17
    assembly universe, then arranges them in a 7x7 core lattice inside a
    cylindrical core barrel with a water downcomer.

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch.
    batches : int
        Total number of batches (active + inactive).
    inactive : int
        Number of inactive batches for fission source convergence.

    Returns
    -------
    openmc.Model
        The complete OpenMC full-core model.
    """

    # =========================================================================
    # MATERIALS
    # =========================================================================

    # --- Three fuel materials (one per enrichment zone) ---
    fuel_materials = {}
    for fuel_type, spec in FUEL_COMPOSITIONS.items():
        mat = openmc.Material(name=spec['name'])
        mat.set_density('sum')
        for nuclide, density in spec['nuclides'].items():
            mat.add_nuclide(nuclide, density)
        fuel_materials[fuel_type] = mat

    # --- Helium gap ---
    # Thin gas-filled gap between fuel pellet and cladding. Accommodates
    # fuel swelling and fission gas release during burnup.
    gap = openmc.Material(name='Helium Gap')
    gap.set_density('g/cc', 0.001598)
    gap.add_nuclide('He4', 1.0)

    # --- Zircaloy-4 cladding ---
    # Used for fuel rod cladding, guide tubes, and instrument tubes.
    # Selected for its low neutron absorption cross section and good
    # corrosion resistance in high-temperature water.
    clad = openmc.Material(name='Zircaloy-4')
    clad.set_density('g/cc', 6.55)
    clad.add_element('Zr', 0.982, 'wo')
    clad.add_element('Sn', 0.015, 'wo')
    clad.add_element('Fe', 0.002, 'wo')
    clad.add_element('Cr', 0.001, 'wo')

    # --- Unborated light water moderator/coolant ---
    # NuScale operates with ZERO soluble boron -- a key design feature.
    # Reactivity control relies on burnable absorbers (Gd2O3) and control
    # rods. This simplifies the chemical volume control system and
    # eliminates boron dilution accident scenarios.
    #
    # Density: 0.725 g/cc at ~543 K average moderator temperature.
    # NuScale operates at lower pressure (~12.76 MPa) and temperature
    # than a standard PWR.
    moderator = openmc.Material(name='Unborated Water')
    moderator.set_density('g/cc', 0.725)
    moderator.add_nuclide('H1', 2.0)
    moderator.add_nuclide('O16', 1.0)
    # S(alpha,beta) thermal scattering law for hydrogen bound in water.
    # CRITICAL for correct neutron thermalization -- without it, the
    # free-gas scattering kernel gives wrong thermal spectrum and biased
    # k-eff (typically hundreds of pcm error).
    moderator.add_s_alpha_beta('c_H_in_H2O')

    # --- SS304 core barrel ---
    # The core barrel is a thick cylindrical shell that separates the
    # fuel assemblies from the downcomer water. It provides structural
    # support and directs coolant flow.
    ss304 = openmc.Material(name='SS304')
    ss304.set_density('g/cc', 8.03)
    ss304.add_element('Fe', 0.695, 'wo')
    ss304.add_element('Cr', 0.190, 'wo')
    ss304.add_element('Ni', 0.095, 'wo')
    ss304.add_element('Mn', 0.020, 'wo')

    # Collect all materials
    all_materials = list(fuel_materials.values()) + [gap, clad, moderator, ss304]
    materials = openmc.Materials(all_materials)
    materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # PIN-LEVEL GEOMETRY
    # =========================================================================

    pin_pitch = 1.26          # cm, center-to-center fuel rod spacing
    assembly_pitch = 21.50    # cm, center-to-center assembly spacing
    active_height = 200.0     # cm, active fuel height

    # --- Radial surfaces for fuel pin ---
    fuel_or_surf = openmc.ZCylinder(r=0.4096, name='Fuel pellet OR')
    clad_ir_surf = openmc.ZCylinder(r=0.418, name='Fuel clad IR')
    clad_or_surf = openmc.ZCylinder(r=0.475, name='Fuel clad OR')

    # --- Radial surfaces for guide tube ---
    gt_ir_surf = openmc.ZCylinder(r=0.561, name='Guide tube IR')
    gt_or_surf = openmc.ZCylinder(r=0.602, name='Guide tube OR')

    # --- Radial surfaces for instrument tube ---
    it_ir_surf = openmc.ZCylinder(r=0.559, name='Instrument tube IR')
    it_or_surf = openmc.ZCylinder(r=0.605, name='Instrument tube OR')

    # =====================================================================
    # FUEL PIN UNIVERSES (one per enrichment type)
    # =====================================================================
    # Each fuel pin consists of concentric cylinders:
    #   1. UO2 fuel pellet (r = 0.4096 cm)
    #   2. Helium gap      (r = 0.418 cm)
    #   3. Zircaloy-4 clad (r = 0.475 cm)
    #   4. Moderator water (fills remaining lattice cell)

    fuel_pin_universes = {}
    for fuel_type, fuel_mat in fuel_materials.items():
        enrichment = {1: '1.6', 2: '2.4', 3: '3.1'}[fuel_type]

        fuel_cell = openmc.Cell(name=f'Fuel pellet {enrichment}%')
        fuel_cell.fill = fuel_mat
        fuel_cell.region = -fuel_or_surf
        fuel_cell.temperature = 543.0  # HZP fuel temperature [K]

        gap_cell = openmc.Cell(name=f'He gap {enrichment}%')
        gap_cell.fill = gap
        gap_cell.region = +fuel_or_surf & -clad_ir_surf
        gap_cell.temperature = 543.0

        clad_cell = openmc.Cell(name=f'Fuel clad {enrichment}%')
        clad_cell.fill = clad
        clad_cell.region = +clad_ir_surf & -clad_or_surf
        clad_cell.temperature = 543.0

        mod_cell = openmc.Cell(name=f'Fuel mod {enrichment}%')
        mod_cell.fill = moderator
        mod_cell.region = +clad_or_surf
        mod_cell.temperature = 543.0

        fuel_pin_universes[fuel_type] = openmc.Universe(
            name=f'Fuel Pin {enrichment}%',
            cells=[fuel_cell, gap_cell, clad_cell, mod_cell]
        )

    # =====================================================================
    # GUIDE TUBE UNIVERSE
    # =====================================================================
    # Water-filled Zircaloy-4 tube at each of 24 guide tube positions.
    # In this ARO (all rods out) model, no control rods are inserted.

    gt_inner = openmc.Cell(name='GT inner water')
    gt_inner.fill = moderator
    gt_inner.region = -gt_ir_surf
    gt_inner.temperature = 543.0

    gt_wall = openmc.Cell(name='GT wall')
    gt_wall.fill = clad
    gt_wall.region = +gt_ir_surf & -gt_or_surf
    gt_wall.temperature = 543.0

    gt_outer = openmc.Cell(name='GT outer water')
    gt_outer.fill = moderator
    gt_outer.region = +gt_or_surf
    gt_outer.temperature = 543.0

    guide_tube_universe = openmc.Universe(
        name='Guide Tube',
        cells=[gt_inner, gt_wall, gt_outer]
    )

    # =====================================================================
    # INSTRUMENT TUBE UNIVERSE
    # =====================================================================
    # Slightly different dimensions from guide tubes. Houses in-core
    # neutron flux detectors (not modeled; tube interior is water-filled).

    it_inner = openmc.Cell(name='IT inner water')
    it_inner.fill = moderator
    it_inner.region = -it_ir_surf
    it_inner.temperature = 543.0

    it_wall = openmc.Cell(name='IT wall')
    it_wall.fill = clad
    it_wall.region = +it_ir_surf & -it_or_surf
    it_wall.temperature = 543.0

    it_outer = openmc.Cell(name='IT outer water')
    it_outer.fill = moderator
    it_outer.region = +it_or_surf
    it_outer.temperature = 543.0

    instrument_tube_universe = openmc.Universe(
        name='Instrument Tube',
        cells=[it_inner, it_wall, it_outer]
    )

    # =====================================================================
    # 17x17 ASSEMBLY UNIVERSES (one per enrichment type)
    # =====================================================================
    # Each assembly type is a RectLattice with the appropriate fuel pin
    # universe at fuel positions, guide tubes at the 24 GT positions,
    # and the instrument tube at the center.

    # Outer universe for the assembly lattices: water fills any space
    # outside the 17x17 pin array but inside the assembly cell. This
    # represents the inter-assembly water gap.
    lattice_outer_universe = openmc.Universe(name='Assembly outer (water)')
    lattice_outer_cell = openmc.Cell(
        name='Inter-assembly water', fill=moderator
    )
    lattice_outer_cell.temperature = 543.0
    lattice_outer_universe.add_cell(lattice_outer_cell)

    assembly_universes = {}
    for fuel_type in [1, 2, 3]:
        enrichment = {1: '1.6', 2: '2.4', 3: '3.1'}[fuel_type]

        lattice = openmc.RectLattice(
            name=f'17x17 Assembly {enrichment}% U-235'
        )
        lattice.pitch = (pin_pitch, pin_pitch)
        lattice.lower_left = (-17 * pin_pitch / 2.0, -17 * pin_pitch / 2.0)

        # Build the 17x17 universe array
        universes = []
        for row in range(17):
            row_list = []
            for col in range(17):
                if (row, col) == INSTRUMENT_TUBE_POSITION:
                    row_list.append(instrument_tube_universe)
                elif (row, col) in GUIDE_TUBE_POSITIONS:
                    row_list.append(guide_tube_universe)
                else:
                    row_list.append(fuel_pin_universes[fuel_type])
            universes.append(row_list)
        lattice.universes = universes
        lattice.outer = lattice_outer_universe

        # Wrap the lattice in a universe with a cell
        assy_cell = openmc.Cell(
            name=f'Assembly {enrichment}%',
            fill=lattice
        )
        assy_universe = openmc.Universe(
            name=f'Assembly {enrichment}%',
            cells=[assy_cell]
        )
        assembly_universes[fuel_type] = assy_universe

    # =====================================================================
    # WATER UNIVERSE (for empty positions in the core lattice)
    # =====================================================================

    water_cell = openmc.Cell(name='Core water', fill=moderator)
    water_cell.temperature = 543.0
    water_universe = openmc.Universe(name='Water', cells=[water_cell])

    # =====================================================================
    # 7x7 CORE LATTICE
    # =====================================================================
    # The 37 fuel assemblies are placed in a 7x7 rectangular lattice.
    # Empty corner positions are filled with the water universe.

    core_lattice = openmc.RectLattice(name='7x7 NuScale Core')
    core_lattice.pitch = (assembly_pitch, assembly_pitch)
    core_lattice.lower_left = (
        -7 * assembly_pitch / 2.0,
        -7 * assembly_pitch / 2.0,
    )

    # Populate the 7x7 lattice from the CORE_MAP
    core_universes = []
    for row in CORE_MAP:
        row_list = []
        for val in row:
            if val == 0:
                row_list.append(water_universe)
            else:
                row_list.append(assembly_universes[val])
        core_universes.append(row_list)
    core_lattice.universes = core_universes

    # Outer universe for the core lattice: water outside the 7x7 grid
    core_lattice.outer = water_universe

    # =====================================================================
    # CORE BARREL AND REACTOR VESSEL
    # =====================================================================
    # The core sits inside a cylindrical core barrel made of SS304.
    # Outside the barrel is the downcomer (water), bounded by the
    # reactor pressure vessel (RPV). We model up to the RPV outer surface.

    core_barrel_ir = 100.0   # cm, core barrel inner radius
    core_barrel_or = 105.0   # cm, core barrel outer radius
    rpv_ir = 115.0           # cm, approximate RPV inner radius

    barrel_ir_surf = openmc.ZCylinder(r=core_barrel_ir, name='Core barrel IR')
    barrel_or_surf = openmc.ZCylinder(r=core_barrel_or, name='Core barrel OR')
    rpv_ir_surf = openmc.ZCylinder(
        r=rpv_ir, boundary_type='vacuum', name='RPV inner surface'
    )

    # Axial boundaries
    z_min = openmc.ZPlane(
        z0=-active_height / 2.0, boundary_type='reflective',
        name='Axial bottom (reflective)'
    )
    z_max = openmc.ZPlane(
        z0=+active_height / 2.0, boundary_type='reflective',
        name='Axial top (reflective)'
    )

    # --- Core region: inside the core barrel, filled with core lattice ---
    core_cell = openmc.Cell(name='Core (fuel assemblies + water)')
    core_cell.fill = core_lattice
    core_cell.region = -barrel_ir_surf & +z_min & -z_max

    # --- Core barrel wall: SS304 annulus ---
    barrel_cell = openmc.Cell(name='Core barrel (SS304)')
    barrel_cell.fill = ss304
    barrel_cell.region = +barrel_ir_surf & -barrel_or_surf & +z_min & -z_max
    barrel_cell.temperature = 543.0

    # --- Downcomer: water between barrel and RPV ---
    downcomer_cell = openmc.Cell(name='Downcomer (water)')
    downcomer_cell.fill = moderator
    downcomer_cell.region = +barrel_or_surf & -rpv_ir_surf & +z_min & -z_max
    downcomer_cell.temperature = 543.0

    # =====================================================================
    # ROOT UNIVERSE AND GEOMETRY
    # =====================================================================

    root_universe = openmc.Universe(
        name='Root',
        cells=[core_cell, barrel_cell, downcomer_cell]
    )
    geometry = openmc.Geometry(root_universe)

    # =========================================================================
    # SETTINGS
    # =========================================================================

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    # Temperature treatment: use nearest available temperature data.
    # The cross section library may not have data at exactly 543 K,
    # so we allow a generous tolerance for interpolation/nearest lookup.
    settings.temperature = {'method': 'nearest', 'tolerance': 1000}

    # Initial fission source: uniform box covering the core region.
    # Only positions that land in fissile material are accepted.
    # The core spans approximately +/- 3.5 * 21.50 = +/- 75.25 cm in XY.
    core_half_width = 3.5 * assembly_pitch  # 75.25 cm
    settings.source = [openmc.IndependentSource(
        space=openmc.stats.Box(
            lower_left=(-core_half_width, -core_half_width,
                        -active_height / 2.0),
            upper_right=(+core_half_width, +core_half_width,
                         +active_height / 2.0)
        )
    )]

    # =========================================================================
    # TALLIES
    # =========================================================================
    # Assembly-level fission rate tally on a 7x7 mesh aligned with the
    # core lattice. This gives relative assembly powers for comparison
    # with reference solutions.

    tallies = openmc.Tallies()

    # Assembly power mesh (7x7, one element per assembly position)
    assy_mesh = openmc.RegularMesh(name='Assembly power mesh')
    assy_mesh.dimension = (7, 7)
    assy_mesh.lower_left = (-core_half_width, -core_half_width)
    assy_mesh.upper_right = (+core_half_width, +core_half_width)

    assy_power_tally = openmc.Tally(name='Assembly Powers')
    assy_power_tally.filters = [openmc.MeshFilter(assy_mesh)]
    assy_power_tally.scores = ['fission']
    tallies.append(assy_power_tally)

    # =========================================================================
    # ASSEMBLE AND RETURN MODEL
    # =========================================================================

    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )
    return model


def main():
    """Parse command-line arguments, build the model, and export/run."""
    parser = argparse.ArgumentParser(
        description='McSAFER NuScale-like SMR Full-Core Benchmark (OpenMC)'
    )
    parser.add_argument(
        '--particles', type=int, default=20000,
        help='Neutron histories per batch (default: 20000)'
    )
    parser.add_argument(
        '--batches', type=int, default=200,
        help='Total batches (default: 200)'
    )
    parser.add_argument(
        '--inactive', type=int, default=50,
        help='Inactive batches for source convergence (default: 50)'
    )
    parser.add_argument(
        '--run', action='store_true',
        help='Run OpenMC after exporting XML files'
    )
    args = parser.parse_args()

    print("Building McSAFER NuScale-like SMR full-core model...")
    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
    )

    # Print core layout summary
    n_type = {1: 0, 2: 0, 3: 0}
    for row in CORE_MAP:
        for val in row:
            if val > 0:
                n_type[val] += 1
    print(f"  Core layout: 7x7 grid, 37 fuel assemblies")
    print(f"    Type 1 (1.6% U-235): {n_type[1]} assemblies (center)")
    print(f"    Type 2 (2.4% U-235): {n_type[2]} assemblies (middle ring)")
    print(f"    Type 3 (3.1% U-235): {n_type[3]} assemblies (outer ring)")
    print(f"  Each assembly: 17x17 (264 fuel + 24 GT + 1 IT)")
    print(f"  Total fuel pins: {37 * 264} = 37 x 264")
    print(f"  Active height: 200 cm")
    print(f"  BORON-FREE moderator (NuScale design feature)")
    print(f"  Conditions: HZP, ARO, fresh fuel")
    print(f"  Cross sections: ENDF/B-VIII.0")

    model.export_to_xml()
    print("\nExported OpenMC XML input files.")

    if args.run:
        print(f"\nRunning OpenMC: {args.particles} particles/batch, "
              f"{args.batches} total batches, {args.inactive} inactive")
        model.run()


if __name__ == '__main__':
    main()
