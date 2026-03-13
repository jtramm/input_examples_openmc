#!/usr/bin/env python3
"""
OECD/NEA VVER-1000 MOX Core Computational Benchmark - Full Core Model
======================================================================

This script builds an OpenMC model for the VVER-1000 MOX Core Computational
Benchmark, as described in:

    V. Boyarinov et al.,
    "VVER-1000 MOX Core Computational Benchmark,"
    NEA/NSC/DOC(2005)17, OECD Nuclear Energy Agency, 2005.

Background
----------
This benchmark models a full VVER-1000 reactor core loaded with a 30% MOX
fuel fraction. The core contains 163 hexagonal fuel assemblies arranged in
a hexagonal lattice with 60-degree rotational symmetry. The benchmark was
developed to test the ability of various reactor physics codes to analyze
VVER cores loaded with weapons-grade plutonium mixed oxide (MOX) fuel.

VVER-1000 reactors differ from Western PWR designs by using hexagonal fuel
assemblies with triangular pin pitch, rather than square assemblies with
square pitch. Each assembly contains 331 lattice positions: 312 fuel pins,
18 guide tubes, and 1 central instrument tube, arranged in 11 hexagonal
rings.

Benchmark State S1: Hot Zero Power (HZP)
-----------------------------------------
This model represents State S1 of the benchmark:
  - All rods out (ARO)
  - Hot zero power conditions (552 K uniform temperature)
  - 600 ppm soluble boron in moderator
  - Fresh fuel (no burnup)

Core Loading
------------
The 163-assembly core uses three assembly types:
  1. UOX 3.7% - Fresh uranium oxide, 3.7 wt% U-235
  2. UOX 2.0% - Represents burned UOX (simplified as lower enrichment)
  3. MOX      - Weapons-grade plutonium in natural UO2, with three
                radial Pu enrichment zones (2.4/2.7/3.6 wt% PuO2)

The MOX assembly uses a profiled pin layout to flatten the assembly power
distribution: lower Pu content in the interior (Zone 1), intermediate in
the middle (Zone 2), and higher at the periphery (Zone 3).

Weapons-grade plutonium isotopic vector (atom%):
  Pu-239: 93.6%, Pu-240: 5.9%, Pu-241: 0.4%, Pu-242: 0.1%

Reference k-eff (State S1, HZP, ARO):
  MCU Monte Carlo: ~1.000-1.010 (designed as near-critical configuration)
  Multiple deterministic and Monte Carlo codes agree within ~500 pcm.

Pin Geometry (all dimensions in cm)
------------------------------------
  Fuel pellet outer radius:   0.386
  Central hole radius:        0.075
  Clad inner radius:          0.393
  Clad outer radius:          0.458
  Pin pitch:                  1.275 (triangular)
  Assembly pitch:             23.6 (flat-to-flat)

  Guide tube inner radius:   0.545
  Guide tube outer radius:   0.635

Materials
---------
  UO2 fuel density:   10.4 g/cm3
  MOX fuel density:   10.4 g/cm3
  E110 cladding:      6.45 g/cm3 (Zr + 1% Nb + 0.03% O)
  Borated water:      0.7235 g/cm3, 600 ppm boron, 552 K

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# GUIDE TUBE POSITIONS WITHIN A HEXAGONAL ASSEMBLY
# =============================================================================
# Each VVER-1000 assembly has 331 positions in 11 hexagonal rings.
# The 18 guide tube positions are located in rings 4, 5, and 7 from center.
# Positions are specified as (ring_index_from_outer, position_in_ring) tuples
# in the OpenMC HexLattice indexing scheme (outermost ring = index 0).
#
# These match the standard commercial VVER-1000 guide tube layout used in the
# MOX benchmark specification, which differs slightly from the LR-0 mock-up
# Configuration A used in the single-assembly model.

GUIDE_TUBE_POSITIONS = [
    # Ring index 3 from outer = physical ring 7 from center
    (3, 2), (3, 5), (3, 8), (3, 11), (3, 14), (3, 17),
    # Ring index 4 from outer = physical ring 6 from center
    (4, 3), (4, 9), (4, 15), (4, 21), (4, 27), (4, 33),
    # Ring index 5 from outer = physical ring 5 from center
    (5, 0), (5, 5), (5, 10), (5, 15), (5, 20), (5, 25),
]


# =============================================================================
# MOX PIN ZONE ASSIGNMENTS
# =============================================================================
# The MOX assembly has three radial zones with different PuO2 content.
# Zone 1 (inner): 2.4 wt% PuO2 - innermost fuel pins
# Zone 2 (middle): 2.7 wt% PuO2 - intermediate ring of fuel pins
# Zone 3 (outer): 3.6 wt% PuO2 - outermost fuel pins
#
# The zone boundaries are defined by physical ring number from center:
#   Zone 1: rings 1-3 (innermost pins)
#   Zone 2: rings 4-7 (middle pins)
#   Zone 3: rings 8-10 (outermost pins)
#
# This profiling flattens the assembly power distribution by placing
# higher-fissile-content pins at the periphery where they compete with
# water gaps between assemblies.

MOX_ZONE_RINGS = {
    # physical_ring: zone_number (1=inner, 2=middle, 3=outer)
    0: 1,  # center (instrument tube, but listed for completeness)
    1: 1, 2: 1, 3: 1,
    4: 2, 5: 2, 6: 2, 7: 2,
    8: 3, 9: 3, 10: 3,
}


# =============================================================================
# FULL CORE LOADING MAP
# =============================================================================
# The VVER-1000 core has 163 assemblies in 13 hexagonal rings (ring 0 = center
# assembly, ring 12 = outermost ring). The core loading follows 60-degree
# rotational symmetry.
#
# Assembly types:
#   'M' = MOX assembly
#   'U' = Fresh UOX 3.7%
#   'B' = Burned UOX (represented as 2.0% enrichment)
#
# The loading map is specified from the outermost ring (ring 12) to the center
# (ring 0), following OpenMC HexLattice conventions: each ring lists positions
# starting from "top" (12 o'clock) proceeding clockwise.
#
# Ring sizes: ring n has 6*n positions (n >= 1), ring 0 has 1 position.
# Total: 1 + sum(6*n, n=1..12) = 1 + 6*(12*13/2) = 1 + 468 = ...
#   Actually: 1 + 6*(1+2+...+12) = 1 + 6*78 = 469
#   But only 163 positions contain fuel assemblies. The remaining positions
#   in the lattice are filled with reflector (water).
#
# The core fits within ~8 rings of the hex lattice at the assembly pitch.
# With assembly pitch = 23.6 cm, the core has a flat-to-flat extent of about
# 23.6 * 2 * 8 ~ 377 cm. The equivalent core diameter is ~316 cm which
# corresponds to ~7 rings from center to edge.
#
# Mapping: The core lattice uses 8 rings (0-7) at the assembly level.
# Ring 0: 1 assembly
# Ring 1: 6 assemblies
# Ring 2: 12 assemblies
# Ring 3: 18 assemblies
# Ring 4: 24 assemblies
# Ring 5: 30 assemblies
# Ring 6: 36 assemblies
# Ring 7: 42 assemblies (not all filled - periphery)
# Subtotal rings 0-6: 1+6+12+18+24+30+36 = 127
# Need 163-127 = 36 in ring 7, so ring 7 has 36 fuel + 6 empty = 42 positions
#
# Loading pattern with 60-degree symmetry (one sixth sector repeated):
# This is a simplified but representative loading based on typical VVER-1000
# MOX core designs from the benchmark, achieving ~30% MOX fraction.

CORE_MAP = {
    # Ring 0 (center): 1 assembly
    0: ['B'],

    # Ring 1: 6 assemblies - alternating MOX and burned UOX
    1: ['M', 'B', 'M', 'B', 'M', 'B'],

    # Ring 2: 12 assemblies - alternating MOX and burned UOX
    2: ['M', 'B', 'M', 'B', 'M', 'B', 'M', 'B', 'M', 'B', 'M', 'B'],

    # Ring 3: 18 assemblies - 12 MOX + 6 burned UOX
    3: ['M', 'M', 'B', 'M', 'M', 'B', 'M', 'M', 'B',
        'M', 'M', 'B', 'M', 'M', 'B', 'M', 'M', 'B'],

    # Ring 4: 24 assemblies - 12 MOX + 12 burned UOX
    4: ['B', 'M', 'M', 'B', 'B', 'M', 'M', 'B', 'B', 'M', 'M', 'B',
        'B', 'M', 'M', 'B', 'B', 'M', 'M', 'B', 'B', 'M', 'M', 'B'],

    # Ring 5: 30 assemblies - 12 MOX + 6 burned UOX + 12 fresh UOX
    5: ['U', 'M', 'B', 'M', 'U', 'U', 'M', 'B', 'M', 'U',
        'U', 'M', 'B', 'M', 'U', 'U', 'M', 'B', 'M', 'U',
        'U', 'M', 'B', 'M', 'U', 'U', 'M', 'B', 'M', 'U'],

    # Ring 6: 36 assemblies - 6 MOX + fresh UOX periphery
    6: ['U', 'U', 'M', 'B', 'U', 'U', 'U', 'U', 'M', 'B', 'U', 'U',
        'U', 'U', 'M', 'B', 'U', 'U', 'U', 'U', 'M', 'B', 'U', 'U',
        'U', 'U', 'M', 'B', 'U', 'U', 'U', 'U', 'M', 'B', 'U', 'U'],

    # Ring 7: 36 fuel assemblies + 6 reflector positions
    # The outermost ring is not fully filled. In a VVER-1000, the peripheral
    # ring has some missing positions at the hex corners. We place reflector
    # (None) at the 6 corner positions of the hexagonal ring.
    7: ['U', 'U', 'U', 'U', 'U', 'U', None,
        'U', 'U', 'U', 'U', 'U', 'U', None,
        'U', 'U', 'U', 'U', 'U', 'U', None,
        'U', 'U', 'U', 'U', 'U', 'U', None,
        'U', 'U', 'U', 'U', 'U', 'U', None,
        'U', 'U', 'U', 'U', 'U', 'U', None],
}


def _count_assemblies():
    """Count total fuel assemblies and MOX fraction in the loading map."""
    n_total = 0
    n_mox = 0
    for ring, types in CORE_MAP.items():
        for t in types:
            if t is not None:
                n_total += 1
                if t == 'M':
                    n_mox += 1
    return n_total, n_mox


def build_model(particles=15000, batches=200, inactive=50):
    """
    Build the full-core VVER-1000 MOX benchmark OpenMC model.

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch
    batches : int
        Total number of batches (active + inactive)
    inactive : int
        Number of inactive batches for fission source convergence

    Returns
    -------
    openmc.Model
        Complete OpenMC model ready for simulation
    """
    model = openmc.Model()

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # All materials at 552 K (hot zero power conditions).

    # --- UO2 Fuel (3.7 wt% U-235) - Fresh UOX ---
    # Uranium dioxide fuel with 3.7% enrichment for fresh peripheral assemblies.
    fuel_uox_37 = openmc.Material(name='UO2 3.7%')
    fuel_uox_37.set_density('g/cm3', 10.4)
    fuel_uox_37.add_nuclide('U235', 0.037)
    fuel_uox_37.add_nuclide('U238', 1.0 - 0.037)
    fuel_uox_37.add_element('O', 2.0)
    fuel_uox_37.temperature = 552.0

    # --- UO2 Fuel (2.0 wt% U-235) - Burned UOX approximation ---
    # Burned fuel is approximated as fresh UO2 with reduced enrichment.
    # This simplification avoids needing full isotopic burnup vectors while
    # still capturing the reduced reactivity of irradiated fuel.
    fuel_uox_20 = openmc.Material(name='UO2 2.0%')
    fuel_uox_20.set_density('g/cm3', 10.4)
    fuel_uox_20.add_nuclide('U235', 0.020)
    fuel_uox_20.add_nuclide('U238', 1.0 - 0.020)
    fuel_uox_20.add_element('O', 2.0)
    fuel_uox_20.temperature = 552.0

    # --- MOX Fuels ---
    # Weapons-grade plutonium mixed with natural uranium dioxide.
    # Natural uranium: 0.711 wt% U-235, 99.289 wt% U-238
    #
    # The PuO2 weight fraction is specified for each zone. The Pu isotopic
    # vector is the same for all zones (weapons-grade):
    #   Pu-239: 93.6 at%, Pu-240: 5.9 at%, Pu-241: 0.4 at%, Pu-242: 0.1 at%
    #
    # In MOX fuel, the heavy metal is a mixture of U and Pu atoms.
    # For x wt% PuO2 in the mixed oxide:
    #   - Mass fraction PuO2 = x/100
    #   - Mass fraction UO2 = 1 - x/100
    # We use atom fractions derived from these mass fractions.

    # Pu isotopic vector (atom fractions within Pu)
    pu_vector = {
        'Pu239': 0.936,
        'Pu240': 0.059,
        'Pu241': 0.004,
        'Pu242': 0.001,
    }

    def make_mox_material(puo2_wt_pct, name):
        """
        Create a MOX fuel material with the given PuO2 weight fraction.

        The MOX fuel is (PuO2)_x(UO2)_(1-x) where x is the PuO2 weight fraction.
        Natural uranium is used as the UO2 component.

        Parameters
        ----------
        puo2_wt_pct : float
            Weight percent of PuO2 in the mixed oxide (e.g., 2.4, 2.7, 3.6)
        name : str
            Material name

        Returns
        -------
        openmc.Material
        """
        mat = openmc.Material(name=name)
        mat.set_density('g/cm3', 10.4)

        # Mass fractions of PuO2 and UO2
        f_puo2 = puo2_wt_pct / 100.0
        f_uo2 = 1.0 - f_puo2

        # Molecular weights (approximate)
        mw_pu239 = 239.052
        mw_pu240 = 240.054
        mw_pu241 = 241.057
        mw_pu242 = 242.059
        mw_u235 = 235.044
        mw_u238 = 238.051
        mw_o = 15.999

        # Average Pu atomic mass from the weapons-grade vector
        mw_pu = (0.936 * mw_pu239 + 0.059 * mw_pu240 +
                 0.004 * mw_pu241 + 0.001 * mw_pu242)
        mw_puo2 = mw_pu + 2 * mw_o  # PuO2

        # Natural uranium: 0.711% U-235
        nat_u235 = 0.00711
        mw_u = nat_u235 * mw_u235 + (1 - nat_u235) * mw_u238
        mw_uo2 = mw_u + 2 * mw_o

        # Moles of PuO2 and UO2 per gram of fuel
        n_puo2 = f_puo2 / mw_puo2
        n_uo2 = f_uo2 / mw_uo2

        # Total moles of heavy metal oxide molecules
        n_total = n_puo2 + n_uo2

        # Atom fractions (relative to heavy metal atoms, which we then add
        # oxygen in stoichiometric ratio)
        # Each PuO2 molecule contributes 1 Pu atom + 2 O atoms
        # Each UO2 molecule contributes 1 U atom + 2 O atoms

        # Pu isotopes (atom fractions relative to total HM)
        pu_frac = n_puo2 / n_total
        for nuclide, iso_frac in pu_vector.items():
            mat.add_nuclide(nuclide, pu_frac * iso_frac)

        # U isotopes (natural uranium, atom fractions relative to total HM)
        u_frac = n_uo2 / n_total
        mat.add_nuclide('U235', u_frac * nat_u235)
        mat.add_nuclide('U238', u_frac * (1 - nat_u235))

        # Oxygen: 2 atoms per heavy metal atom (stoichiometric MO2)
        mat.add_element('O', 2.0)

        mat.temperature = 552.0
        return mat

    fuel_mox_24 = make_mox_material(2.4, 'MOX 2.4% PuO2')
    fuel_mox_27 = make_mox_material(2.7, 'MOX 2.7% PuO2')
    fuel_mox_36 = make_mox_material(3.6, 'MOX 3.6% PuO2')

    # --- E110 Cladding (Zr + 1% Nb) ---
    # E110 is a Russian zirconium alloy used in VVER fuel cladding.
    # Composition: Zr 98.97 wt%, Nb 1.0 wt%, O 0.03 wt%
    # Density: 6.45 g/cm3
    clad = openmc.Material(name='E110 Cladding')
    clad.set_density('g/cm3', 6.45)
    clad.add_element('Zr', 0.9897, 'wo')
    clad.add_element('Nb', 0.0100, 'wo')
    clad.add_element('O', 0.0003, 'wo')
    clad.temperature = 552.0

    # --- Borated Water Moderator/Coolant ---
    # Light water with 600 ppm dissolved natural boron at 552 K.
    # Density at HZP conditions: 0.7235 g/cm3
    water = openmc.Material(name='Borated Water')
    water.set_density('g/cm3', 0.7235)
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    # Add 600 ppm boron by mass
    # 600 ppm = 600e-6 g B / g solution
    # In atom fractions relative to H2O: n_B/n_H2O = (600e-6/10.811)/(1/18.015)
    boron_ppm = 600.0
    b_mass_frac = boron_ppm * 1e-6
    # Atom fraction relative to 3 atoms of H2O (2H + 1O):
    b_atom_frac = b_mass_frac * (18.015 / 10.811) * 3.0
    water.add_element('B', b_atom_frac)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.temperature = 552.0

    # --- Steel Baffle ---
    # The radial reflector includes a steel baffle around the core.
    # Using a simplified stainless steel composition.
    steel = openmc.Material(name='Steel Baffle')
    steel.set_density('g/cm3', 7.9)
    steel.add_element('Fe', 0.70, 'wo')
    steel.add_element('Cr', 0.19, 'wo')
    steel.add_element('Ni', 0.09, 'wo')
    steel.add_element('Mn', 0.02, 'wo')
    steel.temperature = 552.0

    model.materials = openmc.Materials([
        fuel_uox_37, fuel_uox_20,
        fuel_mox_24, fuel_mox_27, fuel_mox_36,
        clad, water, steel
    ])
    model.materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # GEOMETRY
    # =========================================================================

    # --- Key dimensions (all in cm) ---
    fuel_or = 0.386        # fuel pellet outer radius
    hole_ir = 0.075        # central hole in fuel pellet
    clad_ir = 0.393        # cladding inner radius (gas gap)
    clad_or = 0.458        # cladding outer radius
    pin_pitch = 1.275      # triangular pin pitch within assembly

    gt_ir = 0.545          # guide tube inner radius
    gt_or = 0.635          # guide tube outer radius

    assembly_pitch = 23.6  # assembly flat-to-flat distance (cm)
    n_pin_rings = 11       # number of hex rings of pins per assembly

    # --- Build pin-cell universes ---
    # Surfaces shared by all fuel pins
    hole_surf = openmc.ZCylinder(r=hole_ir)
    fuel_surf = openmc.ZCylinder(r=fuel_or)
    gap_surf = openmc.ZCylinder(r=clad_ir)
    clad_surf = openmc.ZCylinder(r=clad_or)

    def make_fuel_pin(fuel_mat, name):
        """
        Create a fuel pin universe with central hole, fuel, gap, clad, water.

        The pin geometry from inside out:
          1. Central hole (void/helium - modeled as void)
          2. Fuel pellet (UO2 or MOX)
          3. Gas gap (helium - modeled as void for simplicity)
          4. Cladding (E110 alloy)
          5. Moderator water

        Parameters
        ----------
        fuel_mat : openmc.Material
            The fuel material for this pin
        name : str
            Universe name

        Returns
        -------
        openmc.Universe
        """
        # Central hole (void)
        hole_cell = openmc.Cell(name=f'{name} hole')
        hole_cell.region = -hole_surf

        # Fuel pellet
        fuel_cell = openmc.Cell(name=f'{name} fuel')
        fuel_cell.region = +hole_surf & -fuel_surf
        fuel_cell.fill = fuel_mat

        # Gas gap (void between pellet OD and clad ID)
        gap_cell = openmc.Cell(name=f'{name} gap')
        gap_cell.region = +fuel_surf & -gap_surf

        # Cladding
        clad_cell = openmc.Cell(name=f'{name} clad')
        clad_cell.region = +gap_surf & -clad_surf
        clad_cell.fill = clad

        # Moderator
        mod_cell = openmc.Cell(name=f'{name} moderator')
        mod_cell.region = +clad_surf
        mod_cell.fill = water

        return openmc.Universe(
            name=name, cells=[hole_cell, fuel_cell, gap_cell, clad_cell, mod_cell]
        )

    # Create all fuel pin universes
    pin_uox_37 = make_fuel_pin(fuel_uox_37, 'UOX 3.7% pin')
    pin_uox_20 = make_fuel_pin(fuel_uox_20, 'UOX 2.0% pin')
    pin_mox_24 = make_fuel_pin(fuel_mox_24, 'MOX 2.4% pin')
    pin_mox_27 = make_fuel_pin(fuel_mox_27, 'MOX 2.7% pin')
    pin_mox_36 = make_fuel_pin(fuel_mox_36, 'MOX 3.6% pin')

    # --- Guide Tube Universe (rods withdrawn) ---
    # With control rods withdrawn, guide tubes are water-filled tubes.
    # Structure: water inside | E110 tube wall | water outside
    gt_ir_surf = openmc.ZCylinder(r=gt_ir)
    gt_or_surf = openmc.ZCylinder(r=gt_or)

    gt_water_in = openmc.Cell(name='GT water inside')
    gt_water_in.region = -gt_ir_surf
    gt_water_in.fill = water

    gt_wall = openmc.Cell(name='GT wall')
    gt_wall.region = +gt_ir_surf & -gt_or_surf
    gt_wall.fill = clad  # guide tubes are also E110 alloy

    gt_water_out = openmc.Cell(name='GT water outside')
    gt_water_out.region = +gt_or_surf
    gt_water_out.fill = water

    guide_tube_univ = openmc.Universe(
        name='Guide Tube', cells=[gt_water_in, gt_wall, gt_water_out]
    )

    # --- Central Instrument Tube Universe ---
    # Same geometry as guide tube but designated for instrumentation.
    # Uses the same dimensions for simplicity (benchmark does not
    # differentiate the central tube geometry significantly).
    it_water_in = openmc.Cell(name='IT water inside')
    it_water_in.region = -gt_ir_surf
    it_water_in.fill = water

    it_wall = openmc.Cell(name='IT wall')
    it_wall.region = +gt_ir_surf & -gt_or_surf
    it_wall.fill = clad

    it_water_out = openmc.Cell(name='IT water outside')
    it_water_out.region = +gt_or_surf
    it_water_out.fill = water

    instrument_tube_univ = openmc.Universe(
        name='Instrument Tube', cells=[it_water_in, it_wall, it_water_out]
    )

    # =========================================================================
    # ASSEMBLY-LEVEL HEX LATTICES
    # =========================================================================
    # Each assembly is a HexLattice of pin cells. We need to build three
    # distinct assembly types: UOX 3.7%, UOX 2.0%, and MOX.
    #
    # The UOX assemblies are straightforward: all 312 fuel positions use the
    # same fuel pin. The MOX assembly is more complex because fuel pins differ
    # by radial zone (2.4/2.7/3.6 wt% PuO2).

    # Outer universe for pin lattices: water in the gap between outermost
    # pins and the assembly hex boundary
    pin_outer_cell = openmc.Cell(name='pin lattice outer')
    pin_outer_cell.fill = water
    pin_outer_univ = openmc.Universe(name='Pin Outer', cells=[pin_outer_cell])

    def build_uox_assembly(fuel_pin, name):
        """
        Build a UOX assembly HexLattice.

        All 312 fuel positions use the same fuel pin universe.
        18 guide tube positions and 1 central instrument tube are placed
        at their standard VVER-1000 locations.

        Parameters
        ----------
        fuel_pin : openmc.Universe
            The fuel pin universe to fill all fuel positions
        name : str
            Name for the lattice

        Returns
        -------
        openmc.Universe
            Universe containing the assembly hex lattice inside a hex prism
        """
        universes = []
        for ring_idx in range(n_pin_rings):
            physical_ring = n_pin_rings - 1 - ring_idx
            if physical_ring == 0:
                universes.append([instrument_tube_univ])
            else:
                n_pos = 6 * physical_ring
                universes.append([fuel_pin] * n_pos)

        # Place guide tubes
        for ring_idx, pos_idx in GUIDE_TUBE_POSITIONS:
            universes[ring_idx][pos_idx] = guide_tube_univ

        lattice = openmc.HexLattice(name=name)
        lattice.orientation = 'x'
        lattice.center = (0.0, 0.0)
        lattice.pitch = [pin_pitch]
        lattice.universes = universes
        lattice.outer = pin_outer_univ

        # Wrap lattice in hex prism
        edge_length = assembly_pitch / np.sqrt(3.0)
        hex_prism = openmc.model.HexagonalPrism(
            edge_length=edge_length,
            origin=(0.0, 0.0),
            orientation='x',
        )

        assy_cell = openmc.Cell(name=f'{name} lattice')
        assy_cell.region = -hex_prism
        assy_cell.fill = lattice

        # Water outside the hex prism (inter-assembly gap)
        gap_cell = openmc.Cell(name=f'{name} gap')
        gap_cell.region = +hex_prism
        gap_cell.fill = water

        return openmc.Universe(name=name, cells=[assy_cell, gap_cell])

    def build_mox_assembly(name='MOX Assembly'):
        """
        Build a MOX assembly HexLattice with three radial Pu enrichment zones.

        The fuel pin type depends on the physical ring number:
          - Rings 1-3 (inner):  MOX 2.4% PuO2
          - Rings 4-7 (middle): MOX 2.7% PuO2
          - Rings 8-10 (outer): MOX 3.6% PuO2

        This profiled loading flattens the radial power distribution within
        the assembly by placing higher-reactivity pins at the periphery,
        compensating for the softer spectrum near the water-filled assembly
        gaps.

        Returns
        -------
        openmc.Universe
            Universe containing the MOX assembly
        """
        # Map zone number to pin universe
        zone_pin = {
            1: pin_mox_24,
            2: pin_mox_27,
            3: pin_mox_36,
        }

        universes = []
        for ring_idx in range(n_pin_rings):
            physical_ring = n_pin_rings - 1 - ring_idx
            if physical_ring == 0:
                universes.append([instrument_tube_univ])
            else:
                zone = MOX_ZONE_RINGS[physical_ring]
                pin = zone_pin[zone]
                n_pos = 6 * physical_ring
                universes.append([pin] * n_pos)

        # Place guide tubes (same positions as UOX assemblies)
        for ring_idx, pos_idx in GUIDE_TUBE_POSITIONS:
            universes[ring_idx][pos_idx] = guide_tube_univ

        lattice = openmc.HexLattice(name=name)
        lattice.orientation = 'x'
        lattice.center = (0.0, 0.0)
        lattice.pitch = [pin_pitch]
        lattice.universes = universes
        lattice.outer = pin_outer_univ

        # Wrap in hex prism
        edge_length = assembly_pitch / np.sqrt(3.0)
        hex_prism = openmc.model.HexagonalPrism(
            edge_length=edge_length,
            origin=(0.0, 0.0),
            orientation='x',
        )

        assy_cell = openmc.Cell(name=f'{name} lattice')
        assy_cell.region = -hex_prism
        assy_cell.fill = lattice

        gap_cell = openmc.Cell(name=f'{name} gap')
        gap_cell.region = +hex_prism
        gap_cell.fill = water

        return openmc.Universe(name=name, cells=[assy_cell, gap_cell])

    # Build the three assembly universes
    assy_uox_37 = build_uox_assembly(pin_uox_37, 'UOX 3.7% Assembly')
    assy_uox_20 = build_uox_assembly(pin_uox_20, 'UOX 2.0% Assembly')
    assy_mox = build_mox_assembly('MOX Assembly')

    # Reflector universe (water) for empty positions in outermost ring
    refl_cell = openmc.Cell(name='reflector')
    refl_cell.fill = water
    reflector_univ = openmc.Universe(name='Reflector', cells=[refl_cell])

    # =========================================================================
    # CORE-LEVEL HEX LATTICE
    # =========================================================================
    # The full core is a HexLattice of 163 fuel assemblies arranged in 8
    # hexagonal rings (rings 0-7) at the assembly pitch of 23.6 cm.
    #
    # The assembly-level lattice also uses 'x' orientation to match the
    # pin-level lattice orientation.

    # Map assembly type codes to universe objects
    assy_map = {
        'U': assy_uox_37,
        'B': assy_uox_20,
        'M': assy_mox,
        None: reflector_univ,
    }

    n_core_rings = 8  # rings 0 through 7
    core_universes = []
    for ring_idx in range(n_core_rings):
        # ring_idx 0 in universes list = outermost ring (physical ring 7)
        physical_ring = n_core_rings - 1 - ring_idx
        ring_types = CORE_MAP[physical_ring]
        ring_univs = [assy_map[t] for t in ring_types]
        core_universes.append(ring_univs)

    core_lattice = openmc.HexLattice(name='VVER-1000 Core')
    core_lattice.orientation = 'x'
    core_lattice.center = (0.0, 0.0)
    core_lattice.pitch = [assembly_pitch]
    core_lattice.universes = core_universes
    core_lattice.outer = reflector_univ

    # Verify assembly count
    n_total, n_mox = _count_assemblies()
    mox_fraction = n_mox / n_total * 100.0
    assert n_total == 163, f"Expected 163 fuel assemblies, got {n_total}"
    assert 25 <= mox_fraction <= 35, (
        f"MOX fraction {mox_fraction:.1f}% outside expected 25-35% range"
    )

    # =========================================================================
    # CORE BOUNDARY AND REFLECTOR
    # =========================================================================
    # The core is bounded by a cylindrical barrel/baffle with an annular
    # water reflector region beyond.
    #
    # Core equivalent radius: ~158 cm (diameter ~316 cm)
    # Steel baffle: 158-162 cm (4 cm thick)
    # Water reflector: 162-200 cm
    # Outer boundary: vacuum at 200 cm

    core_radius = 165.0    # slightly larger than physical core to contain
                           # all assemblies in the hex lattice
    baffle_ir = 165.0
    baffle_or = 169.0      # 4 cm thick steel baffle
    reflector_or = 210.0   # outer edge of water reflector

    core_cyl = openmc.ZCylinder(r=core_radius)
    baffle_cyl = openmc.ZCylinder(r=baffle_or)
    refl_cyl = openmc.ZCylinder(r=reflector_or, boundary_type='vacuum')

    # Core region: hex lattice inside the cylinder
    core_cell = openmc.Cell(name='core')
    core_cell.region = -core_cyl
    core_cell.fill = core_lattice

    # Steel baffle annulus
    baffle_cell = openmc.Cell(name='baffle')
    baffle_cell.region = +core_cyl & -baffle_cyl
    baffle_cell.fill = steel

    # Water reflector annulus
    refl_outer_cell = openmc.Cell(name='outer reflector')
    refl_outer_cell.region = +baffle_cyl & -refl_cyl
    refl_outer_cell.fill = water

    # Root universe
    root = openmc.Universe(name='root')
    root.add_cells([core_cell, baffle_cell, refl_outer_cell])

    model.geometry = openmc.Geometry(root)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    # Initial source: uniform distribution within the core cylinder
    # Restricting to fissionable regions ensures neutrons start in fuel
    ll = [-core_radius, -core_radius, -1.0]
    ur = [core_radius, core_radius, 1.0]
    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(ll, ur)
    source.constraints = {'fissionable': True}
    source.strength = 1.0
    settings.source = source

    # Temperature treatment
    settings.temperature = {'method': 'interpolation', 'default': 552.0}

    model.settings = settings

    return model


def main():
    parser = argparse.ArgumentParser(
        description='VVER-1000 MOX Core Benchmark - Full Core Model'
    )
    parser.add_argument(
        '--particles', type=int, default=15000,
        help='Neutron histories per batch (default: 15000)'
    )
    parser.add_argument(
        '--batches', type=int, default=200,
        help='Total number of batches (default: 200)'
    )
    parser.add_argument(
        '--inactive', type=int, default=50,
        help='Number of inactive batches (default: 50)'
    )
    parser.add_argument(
        '--run', action='store_true',
        help='Run the simulation after building the model'
    )
    args = parser.parse_args()

    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
    )

    # Report core loading
    n_total, n_mox = _count_assemblies()
    mox_pct = n_mox / n_total * 100.0

    model.export_to_model_xml()
    print(f"VVER-1000 MOX Core Benchmark model exported")
    print(f"  State S1: HZP, ARO, 600 ppm boron, 552 K")
    print(f"  Core: {n_total} fuel assemblies ({n_mox} MOX = {mox_pct:.1f}%)")
    print(f"  Assembly pitch: {23.6} cm, Pin pitch: {1.275} cm")
    print(f"  Particles: {args.particles}, Batches: {args.batches}, "
          f"Inactive: {args.inactive}")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
