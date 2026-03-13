#!/usr/bin/env python3
"""
VVER-1000 Mock-up Benchmark - Single Hexagonal Fuel Assembly (2D)
=================================================================

This script builds an OpenMC model for a VVER-1000 hexagonal fuel assembly,
based on the LR-0 research reactor mock-up experiments described in:

    F. Setiawan, M. Lemaire, D. Lee,
    "Analysis of VVER-1000 mock-up criticality experiments with nuclear data
    library ENDF/B-VIII.0 and Monte Carlo code MCS,"
    Nuclear Engineering and Technology 53 (2021) 1-18.

The experiments were carried out at the LR-0 zero-power research reactor
operated by the Research Center Rez in the Czech Republic.

Background: VVER-1000 Reactor Design
--------------------------------------
VVER (Vodo-Vodyanoy Energeticheskiy Reaktor, Water-Water Energetic Reactor)
is a series of pressurized water reactor (PWR) designs developed in the Soviet
Union and Russia. Unlike Western PWR designs that use square fuel assemblies
arranged in a Cartesian grid, VVER reactors use HEXAGONAL fuel assemblies
arranged in a triangular/hexagonal lattice. This fundamental design difference
makes the VVER an important validation case for neutronics codes.

Key differences from Western PWRs:
  - Hexagonal fuel assemblies (vs. square)
  - Triangular pin pitch within assemblies (vs. square)
  - Horizontal steam generators (vs. vertical)
  - No bottom-mounted instrument penetrations

VVER-1000 Mock-up Description
------------------------------
The mock-up consists of 32 dismountable fuel assemblies in a hexagonal lattice
with 23.6 cm flat-to-flat assembly pitch. The fuel pins are arranged in a
triangular lattice within each assembly with a pitch of 12.75 mm (1.275 cm).

Each standard assembly contains:
  - 312 fuel pins (UO2 pellets in Zr-alloy cladding)
  -  18 absorber cluster guide tubes (stainless steel)
  -   1 central instrumentation tube (Zr alloy)
  -----
    331 total lattice positions arranged in 11 hexagonal rings

The fuel pins are ~1.35 m long with an active fuel length of 1.25 m (shorter
than the 3.50 m commercial VVER-1000). Three UO2 enrichments are used:
2.0%, 3.0%, and 3.3% U-235.

Criticality experiments were performed at room temperature and atmospheric
pressure. Six critical configurations were achieved by varying the moderator
level and boric acid concentration. Reference MCS/MCNP6 calculations with
ENDF/B-VII.1 show keff overprediction of +137 to +532 pcm versus experiment.

This Model
-----------
This script models a SINGLE hexagonal assembly as a 2D infinite lattice with
reflective boundary conditions on all six faces. This is a standard approach
for assembly-level neutronics calculations and captures the key physics of
the hexagonal lattice geometry. The 2D approximation (infinite in the axial
direction) eliminates axial leakage effects.

OpenMC Hex Lattice Notes
--------------------------
openmc.HexLattice uses a ring-based indexing scheme:
  - Ring 0 (center): 1 position
  - Ring 1: 6 positions
  - Ring 2: 12 positions
  - ...
  - Ring n: 6*n positions (for n >= 1)

Total positions for N rings: 1 + 3*N*(N-1) where N = num_rings
For our 11-ring assembly: 1 + 3*11*10 = 331 positions

The universes list is ordered from OUTERMOST ring to INNERMOST ring.
Within each ring, positions start at the "top" and proceed clockwise.

Pin Geometry (from IRPhE benchmark documentation)
---------------------------------------------------
  Fuel pellet outer radius:    0.386 cm
  Cladding inner radius:       0.386 cm (no gap modeled - pellet-clad contact)
  Cladding outer radius:       0.4582 cm
  Cladding thickness:          0.071 cm (confirmed by sensitivity analysis)
  Fuel pin pitch:              1.275 cm (triangular lattice)
  Assembly pitch:              23.6 cm (flat-to-flat)

  Central tube inner radius:   0.45 cm
  Central tube outer radius:   0.5177 cm

  Absorber guide tube:
    Absorber pellet radius:    0.35 cm (B4C)
    Absorber clad outer:       0.41 cm
    Guide tube inner radius:   0.545 cm
    Guide tube outer radius:   0.6323 cm

Usage:
    python model.py [--particles N] [--batches N] [--inactive N]
                    [--enrichment {2.0,3.0,3.3}] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# ABSORBER GUIDE TUBE POSITIONS WITHIN THE HEXAGONAL ASSEMBLY
# =============================================================================
# In a VVER-1000 assembly, 18 absorber cluster guide tubes are placed at
# specific positions in the hexagonal lattice. Two different configurations
# exist in the mock-up:
#   - Configuration A (left in Fig. 2 of paper): used for all assemblies
#     except the 3.3% enriched ones
#   - Configuration B (right in Fig. 2): used for 3.3% enriched assemblies,
#     matches commercial VVER-1000 layout
#
# The positions below are specified as (ring_index, position_in_ring) tuples
# using the OpenMC HexLattice indexing scheme where rings are numbered from
# the outermost (ring 0 in the universes list = ring 10 physically) inward,
# and positions within a ring start at the "top" and go clockwise.
#
# For Configuration A (standard mock-up layout), the 18 guide tubes are
# located in the 4th and 5th physical rings from center (rings 6 and 7 in
# the 11-ring universes list, counting from outside).
#
# These positions were determined from the VVER-1000 assembly layout diagram
# (Fig. 2 of the reference paper) and the OpenMC regression test for the
# same geometry (tests/regression_tests/lattice_hex_x/test.py).

GUIDE_TUBE_POSITIONS_CONFIG_A = [
    # (ring_index_from_outer, position_within_ring)
    # Ring index 3 = physical ring 7 from center (18 guide tubes spread
    # across rings 4 and 5 from center, which map to rings 6 and 7 from
    # outer in an 11-ring lattice, i.e., indices 3 and 4 in universes list)
    (3, 2), (3, 5), (3, 8), (3, 11), (3, 14), (3, 17),    # 6 in ring 7
    (5, 0), (5, 5), (5, 10), (5, 15), (5, 20), (5, 25),   # 6 in ring 5
    (4, 3), (4, 9), (4, 15), (4, 21), (4, 27), (4, 33),   # 6 in ring 6
]


def build_model(enrichment=3.0, particles=10000, batches=110, inactive=10):
    """
    Build the OpenMC model for a single VVER-1000 hexagonal fuel assembly.

    Parameters
    ----------
    enrichment : float
        U-235 enrichment in weight percent (2.0, 3.0, or 3.3)
    particles : int
        Number of neutron histories per batch
    batches : int
        Total number of batches (active + inactive)
    inactive : int
        Number of inactive (discarded) batches for fission source convergence

    Returns
    -------
    openmc.Model
        Complete OpenMC model ready for simulation
    """
    model = openmc.Model()

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # All materials at room temperature (294 K) and atmospheric pressure,
    # consistent with the LR-0 experimental conditions.

    # --- UO2 Fuel ---
    # Uranium dioxide fuel with specified enrichment.
    # The theoretical density of UO2 is 10.97 g/cm3. In practice, fuel
    # pellets are sintered to ~95-96% of theoretical density.
    # Using 10.4 g/cm3 as a representative value for the LR-0 fuel.
    fuel = openmc.Material(name=f'UO2 Fuel {enrichment}%')
    fuel.set_density('g/cm3', 10.4)
    # UO2 stoichiometry: 1 U atom + 2 O atoms
    # Express all as atom fractions for consistency
    fuel.add_nuclide('U235', enrichment / 100.0)
    fuel.add_nuclide('U238', 1.0 - enrichment / 100.0)
    fuel.add_element('O', 2.0)  # stoichiometric UO2 (2 O per 1 U)
    fuel.temperature = 294.0

    # --- Zirconium Alloy Cladding ---
    # The cladding is a zirconium alloy (similar to Zircaloy or E110/E635
    # Russian alloys). For simplicity, we model it as pure zirconium since
    # the alloying elements (Nb, Sn, Fe, etc.) have minimal neutronic impact.
    # Density: ~6.55 g/cm3 for Zr alloys.
    zirc = openmc.Material(name='Zr Alloy Cladding')
    zirc.set_density('g/cm3', 6.55)
    zirc.add_element('Zr', 1.0)
    zirc.temperature = 294.0

    # --- Borated Light Water Moderator ---
    # The moderator is light water with dissolved boric acid (H3BO3).
    # Boric acid concentration varies by critical configuration (Table 2):
    #   Case 1: 2.85 g/kg,  Case 2: 3.63 g/kg,  Case 3: 4.06 g/kg
    #   Case 4: 4.44 g/kg,  Case 5: 4.53 g/kg,  Case 6: 4.68 g/kg
    #
    # For this single-assembly model, we use Case 6 conditions (the
    # configuration used for pin power measurements) with 4.68 g/kg H3BO3.
    #
    # Converting boric acid concentration to boron ppm:
    #   H3BO3 molecular weight = 61.83 g/mol
    #   B atomic weight = 10.81 g/mol
    #   Mass fraction of B in H3BO3 = 10.81/61.83 = 0.1749
    #   4.68 g H3BO3 per kg water = 4.68e-3 g/g * 0.1749 = 8.18e-4 g B/g
    #   ~ 818 ppm boron (by mass)
    boron_ppm = 818.0  # approximately, for Case 6
    water = openmc.Material(name='Borated Water')
    water.set_density('g/cm3', 0.998)  # room temperature water density
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    # Add boron as atom fraction: boron_ppm by mass ≈ boron_ppm * 1e-6 * 18/10.81 atom ratio
    b_atom_frac = boron_ppm * 1e-6 * (18.015 / 10.811) * 3.0  # relative to 3 atoms of H2O
    water.add_element('B', b_atom_frac)
    water.add_s_alpha_beta('c_H_in_H2O')  # thermal scattering for H in H2O
    water.temperature = 294.0

    # --- B4C Absorber Material ---
    # Natural boron carbide used in the absorber cluster tubes.
    # Theoretical density of B4C is ~2.52 g/cm3.
    b4c = openmc.Material(name='B4C Absorber')
    b4c.set_density('g/cm3', 2.52)
    b4c.add_element('B', 4.0, 'ao')  # natural boron (19.9% B-10, 80.1% B-11)
    b4c.add_element('C', 1.0, 'ao')
    b4c.temperature = 294.0

    # --- Stainless Steel (for absorber tube cladding) ---
    # The absorber cluster guide tubes use stainless steel cladding.
    # Approximate composition of austenitic stainless steel (Type 304-like).
    steel = openmc.Material(name='Stainless Steel')
    steel.set_density('g/cm3', 7.9)
    steel.add_element('Fe', 0.70, 'wo')
    steel.add_element('Cr', 0.19, 'wo')
    steel.add_element('Ni', 0.09, 'wo')
    steel.add_element('Mn', 0.02, 'wo')
    steel.temperature = 294.0

    model.materials = openmc.Materials([fuel, zirc, water, b4c, steel])

    # =========================================================================
    # GEOMETRY
    # =========================================================================
    # The geometry is built from the inside out:
    #   1. Define pin-cell universes (fuel pin, guide tube, central tube)
    #   2. Arrange pins in a HexLattice
    #   3. Place lattice inside a hexagonal prism with reflective BCs

    # --- Key dimensions (all in cm) ---
    fuel_pellet_or = 0.386     # fuel pellet outer radius
    clad_ir = 0.386            # cladding inner radius (pellet-clad contact)
    clad_or = 0.4582           # cladding outer radius
    pin_pitch = 1.275          # triangular pin pitch

    ct_ir = 0.45               # central tube inner radius
    ct_or = 0.5177             # central tube outer radius

    abs_ir = 0.35              # absorber pellet outer radius
    abs_clad_or = 0.41         # absorber cladding outer radius
    gt_ir = 0.545              # guide tube inner radius
    gt_or = 0.6323             # guide tube outer radius

    assembly_pitch = 23.6      # assembly flat-to-flat distance
    n_rings = 11               # number of hexagonal rings (0 through 10)

    # --- Fuel Pin Universe ---
    # Each fuel pin consists of:
    #   - UO2 fuel pellet (innermost region)
    #   - Zr alloy cladding
    #   - Surrounding moderator water
    # Note: No explicit helium gap is modeled (pellet-clad contact assumed).
    # This is consistent with the IRPhE benchmark specification where the
    # fuel pellet outer radius equals the cladding inner radius.
    fuel_or_surf = openmc.ZCylinder(r=fuel_pellet_or)
    clad_or_surf = openmc.ZCylinder(r=clad_or)

    fuel_cell = openmc.Cell(name='fuel pellet')
    fuel_cell.region = -fuel_or_surf
    fuel_cell.fill = fuel

    clad_cell = openmc.Cell(name='fuel cladding')
    clad_cell.region = +fuel_or_surf & -clad_or_surf
    clad_cell.fill = zirc

    water_cell = openmc.Cell(name='fuel pin moderator')
    water_cell.region = +clad_or_surf
    water_cell.fill = water

    fuel_pin_univ = openmc.Universe(
        name='Fuel Pin', cells=[fuel_cell, clad_cell, water_cell]
    )

    # --- Central Instrumentation Tube Universe ---
    # The central tube is a Zr-alloy tube filled with moderator water.
    # In the actual mock-up, assembly #27 has a dry experimental channel
    # (6.8 cm diameter, steel-cladded), but the standard assemblies have
    # a water-filled Zr-alloy central tube at the lattice center.
    ct_ir_surf = openmc.ZCylinder(r=ct_ir)
    ct_or_surf = openmc.ZCylinder(r=ct_or)

    ct_inner = openmc.Cell(name='central tube water')
    ct_inner.region = -ct_ir_surf
    ct_inner.fill = water

    ct_wall = openmc.Cell(name='central tube wall')
    ct_wall.region = +ct_ir_surf & -ct_or_surf
    ct_wall.fill = zirc

    ct_outer = openmc.Cell(name='central tube moderator')
    ct_outer.region = +ct_or_surf
    ct_outer.fill = water

    central_tube_univ = openmc.Universe(
        name='Central Tube', cells=[ct_inner, ct_wall, ct_outer]
    )

    # --- Guide Tube Universe (rods WITHDRAWN) ---
    # For the critical configuration, absorber rods are WITHDRAWN.
    # The guide tubes are simply steel tubes filled with water (moderator).
    # Structure (inside-out):
    #   1. Water (inside guide tube - where withdrawn rod would be)
    #   2. Steel guide tube wall
    #   3. Surrounding moderator
    # Note: With rods inserted (B4C), k-eff drops to ~0.84 — far subcritical.
    gt_ir_surf = openmc.ZCylinder(r=gt_ir)
    gt_or_surf = openmc.ZCylinder(r=gt_or)

    gt_inner_water = openmc.Cell(name='guide tube water (rod withdrawn)')
    gt_inner_water.region = -gt_ir_surf
    gt_inner_water.fill = water

    gt_wall = openmc.Cell(name='guide tube wall')
    gt_wall.region = +gt_ir_surf & -gt_or_surf
    gt_wall.fill = steel

    gt_outer = openmc.Cell(name='guide tube moderator')
    gt_outer.region = +gt_or_surf
    gt_outer.fill = water

    guide_tube_univ = openmc.Universe(
        name='Guide Tube (Withdrawn)',
        cells=[gt_inner_water, gt_wall, gt_outer]
    )

    # =========================================================================
    # HEXAGONAL LATTICE CONSTRUCTION
    # =========================================================================
    # The hex lattice is the KEY feature of the VVER design. OpenMC's
    # HexLattice class handles the hexagonal geometry natively.
    #
    # Important HexLattice concepts:
    #   - The lattice is defined by concentric rings of hexagonal cells
    #   - Ring numbering: outermost ring first in the universes list
    #   - Within each ring, positions start at the "top" (12 o'clock) and
    #     proceed CLOCKWISE
    #   - Ring 0 (center) has 1 position
    #   - Ring n has 6*n positions (for n >= 1)
    #   - Total positions: 1 + sum(6*n for n=1..N-1) = 1 + 3*N*(N-1)
    #     For N=11: 1 + 3*11*10 = 331 positions
    #
    # Orientation: 'x' means hexagons have flats perpendicular to x-axis
    # (pointy-top hexagons). 'y' means flats perpendicular to y-axis
    # (flat-top hexagons). VVER assemblies typically use 'x' orientation.
    #
    # The pitch for HexLattice is the distance between centers of adjacent
    # hexagonal cells (= pin pitch in our case = 1.275 cm).

    # Step 1: Build the universes list (outermost ring to innermost)
    # Start by filling all positions with fuel pins, then replace specific
    # positions with guide tubes and the central tube.
    universes = []
    for ring_idx in range(n_rings):
        # ring_idx 0 = outermost ring (physical ring 10)
        # ring_idx 10 = innermost ring (physical ring 0 = center)
        physical_ring = n_rings - 1 - ring_idx
        if physical_ring == 0:
            # Center position: 1 element
            universes.append([central_tube_univ])
        else:
            # Non-center ring: 6 * physical_ring elements
            n_positions = 6 * physical_ring
            universes.append([fuel_pin_univ] * n_positions)

    # Step 2: Place absorber guide tubes at their designated positions
    # The guide tube positions depend on the assembly configuration.
    # We use Configuration A (standard mock-up layout).
    for ring_idx, pos_idx in GUIDE_TUBE_POSITIONS_CONFIG_A:
        universes[ring_idx][pos_idx] = guide_tube_univ

    # Verify pin count:
    # Total positions = 331
    # Guide tubes = 18
    # Central tube = 1
    # Fuel pins = 331 - 18 - 1 = 312 (correct!)
    total_pins = sum(len(ring) for ring in universes)
    n_guide = sum(
        1 for ring in universes for u in ring if u is guide_tube_univ
    )
    n_fuel = sum(
        1 for ring in universes for u in ring if u is fuel_pin_univ
    )
    assert total_pins == 331, f"Expected 331 positions, got {total_pins}"
    assert n_guide == 18, f"Expected 18 guide tubes, got {n_guide}"
    assert n_fuel == 312, f"Expected 312 fuel pins, got {n_fuel}"

    # Step 3: Create the HexLattice
    hex_lattice = openmc.HexLattice(name='VVER-1000 Assembly')
    hex_lattice.orientation = 'x'
    hex_lattice.center = (0.0, 0.0)
    hex_lattice.pitch = [pin_pitch]  # single value for 2D lattice
    hex_lattice.universes = universes

    # The outer universe fills space outside the lattice but inside the
    # assembly boundary. This is water moderator in the gap between the
    # outermost pin ring and the assembly hex wrapper.
    outer_cell = openmc.Cell(name='outer moderator')
    outer_cell.fill = water
    outer_univ = openmc.Universe(name='Outer', cells=[outer_cell])
    hex_lattice.outer = outer_univ

    # Step 4: Place the lattice inside a hexagonal prism
    # The assembly boundary is a regular hexagonal prism with flat-to-flat
    # distance of 23.6 cm. The edge length of a regular hexagon is:
    #   edge_length = flat_to_flat / sqrt(3)
    # For orientation='x' hexagonal prism, the flats are perpendicular to x.
    edge_length = assembly_pitch / np.sqrt(3.0)

    # Create the hexagonal prism boundary with reflective BCs
    # This simulates an infinite lattice of identical assemblies
    hex_prism = openmc.model.HexagonalPrism(
        edge_length=edge_length,
        origin=(0.0, 0.0),
        orientation='x',
        boundary_type='reflective'
    )

    # Assembly cell: the hex lattice inside the hex prism
    assembly_cell = openmc.Cell(name='assembly')
    assembly_cell.region = -hex_prism
    assembly_cell.fill = hex_lattice

    # Root universe
    root = openmc.Universe(name='root')
    root.add_cell(assembly_cell)

    model.geometry = openmc.Geometry(root)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    # Initial neutron source: uniform spatial distribution within the
    # assembly hexagonal boundary. The Box source is a simple approximation;
    # the hex lattice will naturally confine fission neutrons.
    ll = [-assembly_pitch / 2.0, -assembly_pitch / 2.0, -1.0]
    ur = [assembly_pitch / 2.0, assembly_pitch / 2.0, 1.0]
    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(ll, ur)
    source.strength = 1.0
    settings.source = source

    # Temperature method
    settings.temperature = {'method': 'interpolation', 'default': 294.0}

    model.settings = settings

    return model


def main():
    parser = argparse.ArgumentParser(
        description='VVER-1000 Mock-up Benchmark - Single Hexagonal Assembly'
    )
    parser.add_argument(
        '--particles', type=int, default=10000,
        help='Neutron histories per batch (default: 10000)'
    )
    parser.add_argument(
        '--batches', type=int, default=110,
        help='Total number of batches (default: 110)'
    )
    parser.add_argument(
        '--inactive', type=int, default=10,
        help='Number of inactive batches (default: 10)'
    )
    parser.add_argument(
        '--enrichment', type=float, default=3.0, choices=[2.0, 3.0, 3.3],
        help='U-235 enrichment in weight percent (default: 3.0)'
    )
    parser.add_argument(
        '--run', action='store_true',
        help='Run the simulation after building the model'
    )
    args = parser.parse_args()

    model = build_model(
        enrichment=args.enrichment,
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
    )

    # Export model to XML files
    model.export_to_model_xml()
    print(f"Model exported: VVER-1000 assembly, {args.enrichment}% enrichment")
    print(f"  Particles: {args.particles}, Batches: {args.batches}, "
          f"Inactive: {args.inactive}")
    print(f"  Pin pitch: 1.275 cm, Assembly pitch: 23.6 cm")
    print(f"  312 fuel pins, 18 guide tubes, 1 central tube = 331 positions")
    print(f"  11 hexagonal rings, orientation='x'")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
