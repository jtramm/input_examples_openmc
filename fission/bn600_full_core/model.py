#!/usr/bin/env python3
"""
BN-600 Sodium-Cooled Fast Reactor - Full-Core Benchmark
========================================================

This script builds an OpenMC model for the BN-600 sodium-cooled fast breeder
reactor hybrid core, based on specifications from:

    IAEA-TECDOC-1623, "Benchmark Analyses on the Natural Circulation Test
    Performed During the PHENIX End-of-Life Experiments" and the companion
    fast reactor benchmark series.

The BN-600 is a real operating fast reactor at the Beloyarsk Nuclear Power
Plant in Russia, in service since 1980. It produces 1470 MWth / 600 MWe.

Background: BN-600 Reactor Design
-----------------------------------
BN-600 is a pool-type sodium-cooled fast breeder reactor. Unlike light water
reactors, fast reactors:
  - Use NO moderator (sodium coolant provides minimal moderation)
  - Operate with a FAST neutron spectrum (most fissions from >100 keV neutrons)
  - Use higher enrichment fuel (17-26% U-235 in BN-600's hybrid core)
  - Have a compact core (active fuel height ~100 cm vs ~400 cm in LWRs)
  - Include breeding blankets (depleted UO2 to produce Pu-239)

The core uses hexagonal fuel assemblies with wire-wrapped pins, a design
typical of sodium-cooled fast reactors worldwide (also used in EBR-II,
FFTF, Phenix, SuperPhenix, Monju, PFBR, CEFR, BN-350, BN-800).

Core Layout
-----------
The BN-600 hybrid core consists of 369 fuel assemblies in three enrichment
zones, surrounded by radial blanket assemblies:

  - LEZ (Low Enrichment Zone):    Inner core, 17.0 wt% U-235
  - MEZ (Medium Enrichment Zone): Middle core, 21.0 wt% U-235
  - HEZ (High Enrichment Zone):   Outer core, 26.0 wt% U-235
  - Radial Blanket:               Depleted UO2 (0.3 wt% U-235)

The three-zone enrichment flattens the radial power distribution, which is
critical for fast reactors where radial leakage is significant.

Above and below the active fuel are axial blanket regions (~30 cm each) of
depleted UO2, which breed Pu-239 from captured neutrons.

Assembly Design
---------------
Each fuel assembly contains 127 fuel pins in a hexagonal arrangement:
  - Pin pitch: 0.75 cm (tight pitch, typical of fast reactors)
  - Assembly flat-to-flat: 9.6 cm (including hex wrapper/duct)
  - Assembly pitch: 9.82 cm (center-to-center, with sodium gap)
  - Wrapper thickness: ~0.2 cm (stainless steel hex duct)

Pin dimensions:
  - Fuel pellet outer radius: 0.29 cm
  - Central hole radius: 0.10 cm (for fission gas / thermal expansion)
  - Clad inner radius: 0.30 cm (0.01 cm sodium-filled gap)
  - Clad outer radius: 0.34 cm
  - Wire wrap not modeled explicitly (smeared into coolant)

OpenMC Modeling Notes
---------------------
This model uses a two-level hex lattice approach:
  1. Pin-level HexLattice: 127 pins within each assembly universe
  2. Core-level HexLattice: assemblies arranged in the full core

Each assembly universe includes axial structure:
  - Lower axial blanket (depleted UO2)
  - Active fuel region (enriched UO2)
  - Upper axial blanket (depleted UO2)

Since this is a fast reactor, NO S(alpha,beta) thermal scattering libraries
are needed. The neutron spectrum is well above the thermal range.

Reference k-eff
----------------
With all control rods withdrawn (ARO), the BN-600 core has excess reactivity
of approximately 4-6% dk/k (k-eff ~ 1.04-1.06). Various codes in the
IAEA benchmark exercise agree within ~500-1000 pcm.

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import math

import numpy as np
import openmc


# =============================================================================
# CORE MAP CONFIGURATION
# =============================================================================
# BN-600 core layout: hexagonal rings from center outward.
# Each ring is assigned a zone type. The total number of assemblies in
# ring n (for n >= 1) is 6*n; ring 0 (center) has 1 assembly.
#
# Zone assignments (approximate, to match ~369 fuel assemblies total):
#   Rings 0-5:  LEZ  -> 1 + 6 + 12 + 18 + 24 + 30 = 91 assemblies
#   Rings 6-8:  MEZ  -> 36 + 42 + 48 = 126 assemblies
#   Rings 9-10: HEZ  -> 54 + 60 = 114 assemblies
#   Total fuel: 91 + 126 + 114 = 331 assemblies
#
# Note: The actual BN-600 has 369 fuel assemblies. Some positions in the
# outer rings contain control rod assemblies, which we model as fuel here
# (all rods out condition). The exact count depends on the specific loading
# pattern. We use 331 fuel + radial blanket for this benchmark.
#
# Ring 11-12: Radial blanket (depleted UO2 assemblies)
#   Ring 11: 66 assemblies, Ring 12: 72 assemblies -> 138 blanket assemblies

ZONE_LEZ = 'LEZ'
ZONE_MEZ = 'MEZ'
ZONE_HEZ = 'HEZ'
ZONE_BLANKET = 'BLANKET'

# Ring-to-zone mapping
RING_ZONES = {
    0: ZONE_LEZ,
    1: ZONE_LEZ,
    2: ZONE_LEZ,
    3: ZONE_LEZ,
    4: ZONE_LEZ,
    5: ZONE_LEZ,
    6: ZONE_MEZ,
    7: ZONE_MEZ,
    8: ZONE_MEZ,
    9: ZONE_HEZ,
    10: ZONE_HEZ,
    11: ZONE_BLANKET,
    12: ZONE_BLANKET,
}

NUM_CORE_RINGS = 13  # rings 0 through 12


def build_model(particles=15000, batches=200, inactive=50):
    """
    Build the full-core BN-600 fast reactor OpenMC model.

    This function constructs:
      1. Materials for three fuel enrichments, blanket, sodium, and steel
      2. Pin-cell universes for each enrichment zone and the blanket
      3. Assembly-level hex lattices (127 pins per assembly)
      4. Core-level hex lattice with proper enrichment zone loading
      5. Axial blanket regions above and below the active fuel
      6. Surrounding sodium reflector region

    Parameters
    ----------
    particles : int
        Number of neutron histories per batch (default: 15000)
    batches : int
        Total number of batches (default: 200)
    inactive : int
        Number of inactive batches for fission source convergence (default: 50)

    Returns
    -------
    openmc.Model
        Complete OpenMC model ready for simulation
    """
    model = openmc.Model()

    # =========================================================================
    # DIMENSIONS (all in cm)
    # =========================================================================
    # Pin geometry
    fuel_hole_r = 0.10       # central hole in fuel pellet
    fuel_or = 0.29           # fuel pellet outer radius
    clad_ir = 0.30           # cladding inner radius
    clad_or = 0.34           # cladding outer radius
    pin_pitch = 0.75         # pin-to-pin center distance

    # Assembly geometry
    assy_flat_to_flat = 9.6  # assembly flat-to-flat (including wrapper)
    assy_pitch = 9.82        # assembly center-to-center pitch
    wrapper_thickness = 0.2  # hex wrapper (duct) wall thickness

    # Inner flat-to-flat of the wrapper (fuel region boundary)
    inner_ftf = assy_flat_to_flat - 2.0 * wrapper_thickness  # 9.2 cm

    # Number of pin rings in assembly
    # 127 pins = 1 + 6 + 12 + 18 + 24 + 30 + 36 = 7 complete hex rings
    # (ring 0 center + rings 1-6 = 1 + 3*6*7 ... wait, let's verify)
    # Ring 0: 1, Ring 1: 6, Ring 2: 12, Ring 3: 18, Ring 4: 24, Ring 5: 30,
    # Ring 6: 36. Total = 1 + 6 + 12 + 18 + 24 + 30 + 36 = 127. Correct!
    n_pin_rings = 7

    # Axial dimensions
    active_height = 100.0       # active fuel height
    axial_blanket_height = 30.0  # upper and lower axial blanket each
    total_height = active_height + 2.0 * axial_blanket_height  # 160 cm

    # Axial plane positions (centered at z=0)
    z_bottom = -total_height / 2.0       # -80 cm
    z_fuel_bottom = -active_height / 2.0  # -50 cm
    z_fuel_top = active_height / 2.0      # +50 cm
    z_top = total_height / 2.0            # +80 cm

    # =========================================================================
    # MATERIALS
    # =========================================================================
    # All materials at 600 K (hot zero power approximation).
    # Fast reactor temperatures are higher than LWRs, but for eigenvalue
    # calculations an average temperature is acceptable.

    temperature = 600.0  # K, HZP approximation for all materials

    # --- UO2 Fuel: LEZ (17.0 wt% U-235) ---
    # UO2 density ~10.5 g/cm3 (slightly below theoretical 10.97 due to
    # porosity and central hole smearing). Fast reactor fuel pellets are
    # typically ~90-95% theoretical density.
    fuel_lez = openmc.Material(name='UO2 LEZ 17%')
    fuel_lez.set_density('g/cm3', 10.5)
    fuel_lez.add_nuclide('U235', 0.17, 'wo')
    fuel_lez.add_nuclide('U238', 0.83, 'wo')
    # Oxygen content: UO2 is ~11.85 wt% oxygen
    # But since we specified U by weight, we add O to make UO2 stoichiometry
    # More precisely: use add_element with percent_type='ao' for atom fractions
    # UO2: 1 U + 2 O atoms. U is already split into isotopes above, so we
    # re-do this properly using atom fractions.
    fuel_lez = openmc.Material(name='UO2 LEZ 17%')
    fuel_lez.set_density('g/cm3', 10.5)
    # Convert weight percent enrichment to atom fraction
    # For 17 wt% U-235 in U: N235/N_U = (0.17/235) / (0.17/235 + 0.83/238)
    enrich_lez = 0.17
    af235_lez = (enrich_lez / 235.0) / (enrich_lez / 235.0 + (1.0 - enrich_lez) / 238.0)
    fuel_lez.add_nuclide('U235', af235_lez, 'ao')
    fuel_lez.add_nuclide('U238', 1.0 - af235_lez, 'ao')
    fuel_lez.add_element('O', 2.0, 'ao')  # 2 O per U in UO2
    fuel_lez.temperature = temperature

    # --- UO2 Fuel: MEZ (21.0 wt% U-235) ---
    fuel_mez = openmc.Material(name='UO2 MEZ 21%')
    fuel_mez.set_density('g/cm3', 10.5)
    enrich_mez = 0.21
    af235_mez = (enrich_mez / 235.0) / (enrich_mez / 235.0 + (1.0 - enrich_mez) / 238.0)
    fuel_mez.add_nuclide('U235', af235_mez, 'ao')
    fuel_mez.add_nuclide('U238', 1.0 - af235_mez, 'ao')
    fuel_mez.add_element('O', 2.0, 'ao')
    fuel_mez.temperature = temperature

    # --- UO2 Fuel: HEZ (26.0 wt% U-235) ---
    fuel_hez = openmc.Material(name='UO2 HEZ 26%')
    fuel_hez.set_density('g/cm3', 10.5)
    enrich_hez = 0.26
    af235_hez = (enrich_hez / 235.0) / (enrich_hez / 235.0 + (1.0 - enrich_hez) / 238.0)
    fuel_hez.add_nuclide('U235', af235_hez, 'ao')
    fuel_hez.add_nuclide('U238', 1.0 - af235_hez, 'ao')
    fuel_hez.add_element('O', 2.0, 'ao')
    fuel_hez.temperature = temperature

    # --- Depleted UO2 Blanket (0.3 wt% U-235) ---
    # Axial and radial blankets use depleted uranium oxide.
    # Slightly lower density than fuel (~10.0 g/cm3) due to different
    # fabrication route and potentially different porosity.
    blanket_mat = openmc.Material(name='UO2 Blanket 0.3%')
    blanket_mat.set_density('g/cm3', 10.0)
    enrich_bkt = 0.003
    af235_bkt = (enrich_bkt / 235.0) / (enrich_bkt / 235.0 + (1.0 - enrich_bkt) / 238.0)
    blanket_mat.add_nuclide('U235', af235_bkt, 'ao')
    blanket_mat.add_nuclide('U238', 1.0 - af235_bkt, 'ao')
    blanket_mat.add_element('O', 2.0, 'ao')
    blanket_mat.temperature = temperature

    # --- Stainless Steel Cladding and Wrapper ---
    # Cr16Ni15Mo3 type austenitic stainless steel, typical of Soviet/Russian
    # fast reactor cladding (similar to Western 316-type steel).
    # Composition: Fe 64%, Cr 16%, Ni 15%, Mo 3%, Mn 1.5%, Si 0.5%
    steel = openmc.Material(name='SS Cladding/Wrapper')
    steel.set_density('g/cm3', 7.9)
    steel.add_element('Fe', 0.640, 'wo')
    steel.add_element('Cr', 0.160, 'wo')
    steel.add_element('Ni', 0.150, 'wo')
    steel.add_element('Mo', 0.030, 'wo')
    steel.add_element('Mn', 0.015, 'wo')
    steel.add_element('Si', 0.005, 'wo')
    steel.temperature = temperature

    # --- Sodium Coolant ---
    # Liquid sodium at operating temperature (~670 K average).
    # Sodium is an excellent fast reactor coolant: high thermal conductivity,
    # low neutron moderation (heavy nucleus, A=23), low absorption.
    # Density at ~670 K: 0.845 g/cm3 (drops from 0.927 at melting point 371 K)
    sodium = openmc.Material(name='Sodium Coolant')
    sodium.set_density('g/cm3', 0.845)
    sodium.add_element('Na', 1.0, 'ao')
    sodium.temperature = temperature

    # Map zone names to fuel materials
    fuel_materials = {
        ZONE_LEZ: fuel_lez,
        ZONE_MEZ: fuel_mez,
        ZONE_HEZ: fuel_hez,
        ZONE_BLANKET: blanket_mat,
    }

    model.materials = openmc.Materials([
        fuel_lez, fuel_mez, fuel_hez, blanket_mat, steel, sodium
    ])
    model.materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # PIN-CELL UNIVERSES
    # =========================================================================
    # Each pin cell has concentric cylindrical regions:
    #   1. Central void/hole (if present in fuel pellet)
    #   2. Fuel pellet (UO2)
    #   3. Sodium-filled gap (between pellet and clad)
    #   4. Cladding (stainless steel)
    #   5. Coolant (sodium)
    #
    # For the blanket pins, the same geometry is used but with depleted UO2.
    # The central hole is present in fuel pellets but not always in blanket
    # pellets; for simplicity we include it in all pins.

    # Cylindrical surfaces shared by all pin types
    hole_surf = openmc.ZCylinder(r=fuel_hole_r)
    fuel_surf = openmc.ZCylinder(r=fuel_or)
    gap_surf = openmc.ZCylinder(r=clad_ir)
    clad_surf = openmc.ZCylinder(r=clad_or)

    # Axial boundary surfaces
    z_bot_surf = openmc.ZPlane(z0=z_bottom)
    z_fbot_surf = openmc.ZPlane(z0=z_fuel_bottom)
    z_ftop_surf = openmc.ZPlane(z0=z_fuel_top)
    z_top_surf = openmc.ZPlane(z0=z_top)

    def make_fuel_pin_universe(fuel_mat, zone_name):
        """
        Create a fuel pin universe with axial blanket regions.

        The pin has three axial sections:
          - Lower blanket (z_bottom to z_fuel_bottom): depleted UO2
          - Active fuel (z_fuel_bottom to z_fuel_top): enriched UO2
          - Upper blanket (z_fuel_top to z_top): depleted UO2

        Each axial section has the same radial structure:
        hole -> fuel/blanket -> gap -> clad -> sodium

        Parameters
        ----------
        fuel_mat : openmc.Material
            The fuel material for the active region
        zone_name : str
            Zone identifier for naming

        Returns
        -------
        openmc.Universe
            Pin-cell universe with axial structure
        """
        # Axial regions
        lower_blanket_region = +z_bot_surf & -z_fbot_surf
        active_region = +z_fbot_surf & -z_ftop_surf
        upper_blanket_region = +z_ftop_surf & -z_top_surf
        full_axial = +z_bot_surf & -z_top_surf

        cells = []

        # --- Active fuel region ---
        # Central hole (void/sodium-filled in reality; model as sodium)
        c = openmc.Cell(name=f'{zone_name} fuel hole')
        c.region = -hole_surf & active_region
        c.fill = sodium
        cells.append(c)

        # Fuel pellet
        c = openmc.Cell(name=f'{zone_name} fuel pellet')
        c.region = +hole_surf & -fuel_surf & active_region
        c.fill = fuel_mat
        cells.append(c)

        # Sodium gap
        c = openmc.Cell(name=f'{zone_name} fuel gap')
        c.region = +fuel_surf & -gap_surf & active_region
        c.fill = sodium
        cells.append(c)

        # Cladding
        c = openmc.Cell(name=f'{zone_name} fuel clad')
        c.region = +gap_surf & -clad_surf & active_region
        c.fill = steel
        cells.append(c)

        # Coolant around pin in active region
        c = openmc.Cell(name=f'{zone_name} fuel coolant')
        c.region = +clad_surf & active_region
        c.fill = sodium
        cells.append(c)

        # --- Lower axial blanket ---
        c = openmc.Cell(name=f'{zone_name} lower blanket hole')
        c.region = -hole_surf & lower_blanket_region
        c.fill = sodium
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} lower blanket pellet')
        c.region = +hole_surf & -fuel_surf & lower_blanket_region
        c.fill = blanket_mat
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} lower blanket gap')
        c.region = +fuel_surf & -gap_surf & lower_blanket_region
        c.fill = sodium
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} lower blanket clad')
        c.region = +gap_surf & -clad_surf & lower_blanket_region
        c.fill = steel
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} lower blanket coolant')
        c.region = +clad_surf & lower_blanket_region
        c.fill = sodium
        cells.append(c)

        # --- Upper axial blanket ---
        c = openmc.Cell(name=f'{zone_name} upper blanket hole')
        c.region = -hole_surf & upper_blanket_region
        c.fill = sodium
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} upper blanket pellet')
        c.region = +hole_surf & -fuel_surf & upper_blanket_region
        c.fill = blanket_mat
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} upper blanket gap')
        c.region = +fuel_surf & -gap_surf & upper_blanket_region
        c.fill = sodium
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} upper blanket clad')
        c.region = +gap_surf & -clad_surf & upper_blanket_region
        c.fill = steel
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} upper blanket coolant')
        c.region = +clad_surf & upper_blanket_region
        c.fill = sodium
        cells.append(c)

        return openmc.Universe(name=f'{zone_name} Fuel Pin', cells=cells)

    # Create pin universes for each zone
    pin_universes = {}
    for zone_name, fuel_mat in fuel_materials.items():
        pin_universes[zone_name] = make_fuel_pin_universe(fuel_mat, zone_name)

    # =========================================================================
    # ASSEMBLY UNIVERSES
    # =========================================================================
    # Each assembly is a hex lattice of 127 pins inside a hex wrapper (duct).
    # The wrapper is a stainless steel hex tube. Between the outermost pin
    # ring and the wrapper is sodium coolant, and between adjacent assembly
    # wrappers is an inter-assembly sodium gap.
    #
    # Structure (inside-out):
    #   1. Pin hex lattice (127 pins, 7 rings)
    #   2. Sodium between outer pins and wrapper
    #   3. Stainless steel wrapper tube
    #   4. Sodium inter-assembly gap

    def make_assembly_universe(zone_name, pin_univ):
        """
        Create an assembly universe with hex pin lattice and wrapper.

        Parameters
        ----------
        zone_name : str
            Zone identifier for naming
        pin_univ : openmc.Universe
            The pin-cell universe to fill the lattice

        Returns
        -------
        openmc.Universe
            Assembly universe
        """
        # --- Pin hex lattice ---
        pin_lattice = openmc.HexLattice(name=f'{zone_name} Pin Lattice')
        pin_lattice.orientation = 'y'  # flat-top hexagons
        pin_lattice.center = (0.0, 0.0, 0.0)
        pin_lattice.pitch = [pin_pitch]  # 2D lattice (axial handled by pin universe)

        # Build universes list: outermost ring first
        # Ring 0 (outermost, index 0 in list) = physical ring 6 -> 36 pins
        # ...
        # Ring 6 (innermost, index 6 in list) = physical ring 0 -> 1 pin (center)
        univs = []
        for ring_idx in range(n_pin_rings):
            physical_ring = n_pin_rings - 1 - ring_idx
            if physical_ring == 0:
                univs.append([pin_univ])
            else:
                n_pos = 6 * physical_ring
                univs.append([pin_univ] * n_pos)

        pin_lattice.universes = univs

        # Outer universe for the pin lattice: sodium coolant fills space
        # outside the pin lattice but inside the wrapper
        sodium_outer_cell = openmc.Cell(name=f'{zone_name} lattice outer Na')
        sodium_outer_cell.fill = sodium
        sodium_outer_univ = openmc.Universe(
            name=f'{zone_name} lattice outer',
            cells=[sodium_outer_cell]
        )
        pin_lattice.outer = sodium_outer_univ

        # --- Wrapper (hex duct) ---
        # The wrapper is a hex tube. We define inner and outer hex prisms.
        # Inner hex prism: contains the pin lattice + sodium gap
        # Outer hex prism: outer surface of the steel wrapper
        inner_edge = inner_ftf / math.sqrt(3.0)
        outer_edge = assy_flat_to_flat / math.sqrt(3.0)

        inner_hex = openmc.model.HexagonalPrism(
            edge_length=inner_edge,
            orientation='y',
            origin=(0.0, 0.0),
        )
        outer_hex = openmc.model.HexagonalPrism(
            edge_length=outer_edge,
            orientation='y',
            origin=(0.0, 0.0),
        )

        # Axial extent of the assembly
        full_axial = +z_bot_surf & -z_top_surf

        # Cell 1: Pin lattice region (inside inner hex, full axial extent)
        lattice_cell = openmc.Cell(name=f'{zone_name} pin lattice')
        lattice_cell.region = -inner_hex & full_axial
        lattice_cell.fill = pin_lattice

        # Cell 2: Wrapper wall (between inner and outer hex)
        wrapper_cell = openmc.Cell(name=f'{zone_name} wrapper')
        wrapper_cell.region = +inner_hex & -outer_hex & full_axial
        wrapper_cell.fill = steel

        # Cell 3: Sodium above assembly
        above_cell = openmc.Cell(name=f'{zone_name} above')
        above_cell.region = -outer_hex & +z_top_surf
        above_cell.fill = sodium

        # Cell 4: Sodium below assembly
        below_cell = openmc.Cell(name=f'{zone_name} below')
        below_cell.region = -outer_hex & -z_bot_surf
        below_cell.fill = sodium

        # Cell 5: Inter-assembly sodium (outside outer hex)
        # This is bounded by the core lattice cell, so we just fill
        # everything outside the wrapper with sodium
        gap_cell = openmc.Cell(name=f'{zone_name} inter-assy Na')
        gap_cell.region = +outer_hex
        gap_cell.fill = sodium

        return openmc.Universe(
            name=f'{zone_name} Assembly',
            cells=[lattice_cell, wrapper_cell, above_cell, below_cell,
                   gap_cell]
        )

    # Create assembly universes for each zone
    assembly_universes = {}
    for zone_name in [ZONE_LEZ, ZONE_MEZ, ZONE_HEZ, ZONE_BLANKET]:
        assembly_universes[zone_name] = make_assembly_universe(
            zone_name, pin_universes[zone_name]
        )

    # =========================================================================
    # CORE-LEVEL HEXAGONAL LATTICE
    # =========================================================================
    # The core is a large hex lattice of assembly universes.
    # Ring assignments follow RING_ZONES defined at the top of this file.
    #
    # OpenMC HexLattice universes: outermost ring first, clockwise from top.

    core_lattice = openmc.HexLattice(name='BN-600 Core')
    core_lattice.orientation = 'y'
    core_lattice.center = (0.0, 0.0, 0.0)
    core_lattice.pitch = [assy_pitch]

    # Build the core universes list
    core_univs = []
    for ring_idx in range(NUM_CORE_RINGS):
        # ring_idx 0 = outermost ring (physical ring 12)
        physical_ring = NUM_CORE_RINGS - 1 - ring_idx
        zone = RING_ZONES[physical_ring]
        assy_univ = assembly_universes[zone]

        if physical_ring == 0:
            core_univs.append([assy_univ])
        else:
            n_pos = 6 * physical_ring
            core_univs.append([assy_univ] * n_pos)

    core_lattice.universes = core_univs

    # Outer universe for the core lattice: sodium reflector
    # Space outside the outermost ring of assemblies
    reflector_cell = openmc.Cell(name='Sodium reflector')
    reflector_cell.fill = sodium
    reflector_univ = openmc.Universe(name='Reflector', cells=[reflector_cell])
    core_lattice.outer = reflector_univ

    # Count assemblies per zone for verification
    n_lez = sum(max(6 * (NUM_CORE_RINGS - 1 - ri), 1)
                for ri in range(NUM_CORE_RINGS)
                if RING_ZONES[NUM_CORE_RINGS - 1 - ri] == ZONE_LEZ)
    n_mez = sum(6 * (NUM_CORE_RINGS - 1 - ri)
                for ri in range(NUM_CORE_RINGS)
                if RING_ZONES[NUM_CORE_RINGS - 1 - ri] == ZONE_MEZ)
    n_hez = sum(6 * (NUM_CORE_RINGS - 1 - ri)
                for ri in range(NUM_CORE_RINGS)
                if RING_ZONES[NUM_CORE_RINGS - 1 - ri] == ZONE_HEZ)
    n_bkt = sum(6 * (NUM_CORE_RINGS - 1 - ri)
                for ri in range(NUM_CORE_RINGS)
                if RING_ZONES[NUM_CORE_RINGS - 1 - ri] == ZONE_BLANKET)
    n_fuel_total = n_lez + n_mez + n_hez

    # =========================================================================
    # ROOT UNIVERSE AND BOUNDARY CONDITIONS
    # =========================================================================
    # The core lattice is placed inside a cylindrical boundary with vacuum BCs.
    # The cylinder radius is chosen to encompass the outermost blanket ring
    # plus some sodium reflector.

    # Radial boundary: outermost assembly center is at radius
    # (NUM_CORE_RINGS - 1) * pitch from center. Add half a flat-to-flat plus
    # a margin for the sodium reflector.
    core_radius = (NUM_CORE_RINGS - 0.5) * assy_pitch + 5.0  # extra 5 cm sodium
    radial_boundary = openmc.ZCylinder(
        r=core_radius, boundary_type='vacuum'
    )

    # Axial boundaries: generous extent beyond the assembly
    axial_margin = 20.0  # 20 cm of sodium above/below assembly ends
    z_lower_bound = openmc.ZPlane(
        z0=z_bottom - axial_margin, boundary_type='vacuum'
    )
    z_upper_bound = openmc.ZPlane(
        z0=z_top + axial_margin, boundary_type='vacuum'
    )

    # Core cell: hex lattice inside the cylindrical boundary
    core_cell = openmc.Cell(name='Core')
    core_cell.region = -radial_boundary & +z_lower_bound & -z_upper_bound
    core_cell.fill = core_lattice

    root = openmc.Universe(name='Root')
    root.add_cell(core_cell)

    model.geometry = openmc.Geometry(root)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    # Initial neutron source: uniform box source in the active core region.
    # The box spans the approximate radial extent of the fuel zones and the
    # active fuel height. This is a rough starting distribution; the inactive
    # batches will converge the fission source to the true distribution.
    source_radius = 10 * assy_pitch  # approximate HEZ outer radius
    ll = [-source_radius, -source_radius, z_fuel_bottom]
    ur = [source_radius, source_radius, z_fuel_top]
    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(ll, ur)
    source.constraints = {'fissionable': True}
    source.strength = 1.0
    settings.source = source

    # Temperature treatment
    # For fast reactors, temperature affects Doppler broadening of resonances
    # (particularly U-238 capture resonances which are important for the
    # negative Doppler feedback coefficient).
    settings.temperature = {'method': 'interpolation', 'default': temperature}

    model.settings = settings

    # =========================================================================
    # PRINT SUMMARY
    # =========================================================================
    print("BN-600 Full-Core Fast Reactor Model")
    print("=" * 60)
    print(f"  Core layout: {NUM_CORE_RINGS} hexagonal rings")
    print(f"  LEZ assemblies (17.0% U-235): {n_lez}")
    print(f"  MEZ assemblies (21.0% U-235): {n_mez}")
    print(f"  HEZ assemblies (26.0% U-235): {n_hez}")
    print(f"  Total fuel assemblies:         {n_fuel_total}")
    print(f"  Radial blanket assemblies:     {n_bkt}")
    print(f"  Pins per assembly:             127")
    print(f"  Pin pitch:                     {pin_pitch} cm")
    print(f"  Assembly pitch:                {assy_pitch} cm")
    print(f"  Active fuel height:            {active_height} cm")
    print(f"  Axial blanket (each):          {axial_blanket_height} cm")
    print(f"  Total axial extent:            {total_height} cm")
    print(f"  Particles/batch:               {particles}")
    print(f"  Batches (total/inactive):      {batches}/{inactive}")
    print(f"  Temperature:                   {temperature} K (HZP)")
    print("=" * 60)

    return model


def main():
    parser = argparse.ArgumentParser(
        description='BN-600 Sodium-Cooled Fast Reactor Full-Core Benchmark'
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

    model.export_to_model_xml()
    print("\nModel exported to model.xml")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
