#!/usr/bin/env python3
"""
BN-600 Sodium-Cooled Fast Reactor - Hybrid Core Benchmark
==========================================================

IAEA-TECDOC-1623: "Benchmark Analyses on the Natural Circulation Test
Performed During the PHENIX End-of-Life Experiments" / BN-600 Hybrid Core
Benchmark Analysis.

BN-600 is a pool-type sodium-cooled fast breeder reactor at Beloyarsk NPP,
Russia. 1470 MWth / 600 MWe, in service since 1980.

The hybrid core configuration contains UO2 fuel at three enrichments plus
a MOX zone, surrounded by steel shielding assemblies (no radial blanket).

Core Layout (Ring-Based Approximation)
--------------------------------------
  Rings 0-5:  LEZ  (Low Enrichment Zone, 17 wt% U-235)    91 assemblies
  Ring  6:    MEZ  (Medium Enrichment Zone, 21 wt% U-235)  36 assemblies
  Ring  7:    MOX  (Mixed Oxide, ~3.4% fissile Pu)         42 assemblies
  Rings 8-9:  HEZ  (High Enrichment Zone, 26 wt% U-235)  102 assemblies
  Ring  10:   SSA1 (Steel Shielding, 1st row)               60 assemblies
  Ring  11:   SSA2 (Steel Shielding, outer rows)             66 assemblies
  Ring  12:   Radial Reflector (steel + sodium)              72 assemblies

  Note: The actual benchmark specifies 369 fuel assemblies in a complex
  loading pattern with control rod positions. This ring-based approximation
  gives 271 fuel assemblies. The compositions and dimensions match the
  benchmark specification (Table 3.55, heterogeneous model).

Assembly Design
---------------
  127 fuel pins per assembly (7 hex rings: 0-6)
  Pin pitch: 0.795 cm (triangular lattice)
  Wire wrap spacer omitted; cladding density increased to compensate
  Hex wrapper: 9.6 cm flat-to-flat external, 0.2 cm wall thickness

Pin Dimensions (from benchmark, cold dimensions at 20C)
--------------------------------------------------------
  Fuel pellet outer radius: 0.305 cm (no central hole, no gap)
  Cladding inner radius:    0.305 cm (pellet = clad IR)
  Cladding outer radius:    0.345 cm

Axial Structure (ASYMMETRIC blankets)
--------------------------------------
  Lower axial reflector:  20.0 cm  (steel + sodium)
  Lower axial blanket:    35.2 cm  (depleted UO2)
  Active fuel:           104.4 cm  (enriched UO2 or MOX)
  Upper axial blanket:    30.2 cm  (depleted UO2)
  Upper axial reflector:  20.0 cm  (steel + sodium)
  Total:                 209.8 cm

Temperatures
------------
  Fuel isotopes (U, Pu, O, FP):  1500 K
  Structure and coolant:           600 K

Number Densities
----------------
  Heterogeneous (pin-level) number densities from IAEA-TECDOC-1623
  Table 3.55 for equilibrium-cycle compositions (includes burnup products:
  Pu-239/240/241, U-236, lumped fission product modeled as Mo95).

Reference k-eff
----------------
  Benchmark codes give k-eff ~ 1.003-1.010 for the reference state
  (SHR at mid-plane). For all-rods-out, k-eff ~ 1.04-1.06.

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import math

import numpy as np
import openmc


# =============================================================================
# ZONE DEFINITIONS
# =============================================================================
ZONE_LEZ = 'LEZ'
ZONE_MEZ = 'MEZ'
ZONE_MOX = 'MOX'
ZONE_HEZ = 'HEZ'
ZONE_SSA1 = 'SSA1'
ZONE_SSA2 = 'SSA2'
ZONE_RR = 'RR'

# Ring-to-zone mapping based on RZ model radial zones from IAEA-TECDOC-1623
# LEZ: 0-49.6 cm (rings 0-5), MEZ: 53-64 cm (ring 6),
# MOX: 64-77 cm (ring 7), HEZ: 77-87 cm (rings 8-9),
# SSA1: 87-96 cm (ring 10), SSA2: 96-122 cm (ring 11),
# RR: beyond (ring 12)
RING_ZONES = {
    0: ZONE_LEZ,
    1: ZONE_LEZ,
    2: ZONE_LEZ,
    3: ZONE_LEZ,
    4: ZONE_LEZ,
    5: ZONE_LEZ,
    6: ZONE_MEZ,
    7: ZONE_MOX,
    8: ZONE_HEZ,
    9: ZONE_HEZ,
    10: ZONE_SSA1,
    11: ZONE_SSA2,
    12: ZONE_RR,
}

NUM_CORE_RINGS = 13  # rings 0 through 12
FUEL_ZONES = {ZONE_LEZ, ZONE_MEZ, ZONE_MOX, ZONE_HEZ}


def build_model(particles=15000, batches=200, inactive=50):
    """Build the BN-600 hybrid core benchmark model."""
    model = openmc.Model()

    # =========================================================================
    # DIMENSIONS (all in cm, cold dimensions at 20C per benchmark spec)
    # =========================================================================
    # Pin geometry (benchmark: pellet diameter = clad inner diameter, no gap)
    fuel_or = 0.305          # fuel pellet outer radius (6.1 mm diameter)
    clad_ir = 0.305          # cladding inner radius = pellet OR (no gap)
    clad_or = 0.345          # cladding outer radius (6.9 mm diameter)
    pin_pitch = 0.795        # pin-to-pin center distance (7.95 mm)

    # Assembly geometry
    assy_ftf = 9.6           # wrapper external flat-to-flat (96 mm)
    wrapper_thick = 0.2      # wrapper wall thickness (2 mm)
    inner_ftf = assy_ftf - 2.0 * wrapper_thick  # 9.2 cm
    assy_pitch = 9.902       # assembly center-to-center (99.02 mm)

    # Pin lattice: 127 pins = 7 hex rings (0-6)
    n_pin_rings = 7

    # Axial dimensions (ASYMMETRIC blankets per TECDOC-1623 Fig. 3.2)
    active_height = 104.4     # active fuel height
    lower_blanket_h = 35.2    # lower axial blanket
    upper_blanket_h = 30.2    # upper axial blanket
    lower_refl_h = 20.0       # lower steel/sodium reflector
    upper_refl_h = 20.0       # upper steel/sodium reflector

    total_height = (lower_refl_h + lower_blanket_h + active_height
                    + upper_blanket_h + upper_refl_h)  # 209.8 cm

    # Axial plane positions (z=0 at bottom of lower reflector)
    z_bot = 0.0
    z_lb_bot = lower_refl_h                           # 20.0
    z_fuel_bot = z_lb_bot + lower_blanket_h            # 55.2
    z_fuel_top = z_fuel_bot + active_height            # 159.6
    z_ub_top = z_fuel_top + upper_blanket_h            # 189.8
    z_top = z_ub_top + upper_refl_h                    # 209.8

    # =========================================================================
    # MATERIALS (heterogeneous number densities from TECDOC-1623 Table 3.55)
    # =========================================================================
    # All number densities in units of 10^24 atoms/cm^3 (= atoms/barn-cm)
    # Lumped fission product modeled as Mo95

    fuel_temp = 1500.0  # K - fuel isotopes
    struct_temp = 600.0  # K - structure and coolant

    # --- LEZ Fuel Pellet (Composition 1 - LEZ1, 17 wt% enrichment) ---
    fuel_lez = openmc.Material(name='LEZ Fuel')
    fuel_lez.set_density('sum')
    fuel_lez.add_nuclide('U235', 2.674E-03)
    fuel_lez.add_nuclide('U236', 1.218E-04)
    fuel_lez.add_nuclide('U238', 1.527E-02)
    fuel_lez.add_nuclide('Pu239', 3.184E-04)
    fuel_lez.add_nuclide('Pu240', 1.002E-05)
    fuel_lez.add_nuclide('Pu241', 2.633E-07)
    fuel_lez.add_nuclide('Mo95', 5.552E-04)  # lumped FP
    fuel_lez.add_nuclide('O16', 3.793E-02)
    fuel_lez.temperature = fuel_temp

    # --- MEZ Fuel Pellet (Composition 3, 21 wt% enrichment) ---
    fuel_mez = openmc.Material(name='MEZ Fuel')
    fuel_mez.set_density('sum')
    fuel_mez.add_nuclide('U235', 3.397E-03)
    fuel_mez.add_nuclide('U236', 1.284E-04)
    fuel_mez.add_nuclide('U238', 1.459E-02)
    fuel_mez.add_nuclide('Pu239', 2.599E-04)
    fuel_mez.add_nuclide('Pu240', 6.584E-06)
    fuel_mez.add_nuclide('Pu241', 1.349E-07)
    fuel_mez.add_nuclide('Mo95', 5.788E-04)  # lumped FP
    fuel_mez.add_nuclide('O16', 3.793E-02)
    fuel_mez.temperature = fuel_temp

    # --- MOX Fuel Pellet (Composition 4, ~3.4% fissile Pu) ---
    fuel_mox = openmc.Material(name='MOX Fuel')
    fuel_mox.set_density('sum')
    fuel_mox.add_nuclide('U235', 3.857E-05)
    fuel_mox.add_nuclide('U236', 1.328E-06)
    fuel_mox.add_nuclide('U238', 1.459E-02)
    fuel_mox.add_nuclide('Pu239', 3.441E-03)
    fuel_mox.add_nuclide('Pu240', 3.072E-04)
    fuel_mox.add_nuclide('Pu241', 1.961E-05)
    fuel_mox.add_nuclide('Pu242', 1.570E-06)
    fuel_mox.add_nuclide('Mo95', 5.099E-04)  # lumped FP
    fuel_mox.add_nuclide('O16', 3.781E-02)
    fuel_mox.temperature = fuel_temp

    # --- HEZ Fuel Pellet (Composition 5, 26 wt% enrichment) ---
    fuel_hez = openmc.Material(name='HEZ Fuel')
    fuel_hez.set_density('sum')
    fuel_hez.add_nuclide('U235', 4.420E-03)
    fuel_hez.add_nuclide('U236', 1.180E-04)
    fuel_hez.add_nuclide('U238', 1.376E-02)
    fuel_hez.add_nuclide('Pu239', 1.736E-04)
    fuel_hez.add_nuclide('Pu240', 3.136E-06)
    fuel_hez.add_nuclide('Pu241', 8.700E-08)
    fuel_hez.add_nuclide('Mo95', 4.861E-04)  # lumped FP
    fuel_hez.add_nuclide('O16', 3.793E-02)
    fuel_hez.temperature = fuel_temp

    # --- Depleted UO2 Blanket Pellet ---
    # Approximately 0.3 wt% U-235, same pellet density as fissile fuel
    # Derived from homogeneous blanket compositions (Table 3.2)
    fuel_blanket = openmc.Material(name='Blanket UO2')
    fuel_blanket.set_density('sum')
    fuel_blanket.add_nuclide('U235', 5.70E-05)
    fuel_blanket.add_nuclide('U238', 1.894E-02)
    fuel_blanket.add_nuclide('O16', 3.800E-02)
    fuel_blanket.temperature = fuel_temp

    # --- Steel Cladding/Wrapper (Composition 26, boosted for wire wrap) ---
    # CrNiMo austenitic steel; density increased to compensate for omitted
    # wire wrap spacer per benchmark specification.
    steel = openmc.Material(name='Steel Clad/Wrapper')
    steel.set_density('sum')
    steel.add_nuclide('Fe56', 7.064E-02)
    steel.add_nuclide('Cr52', 1.096E-02)
    steel.add_nuclide('Ni58', 1.580E-04)
    steel.add_nuclide('Mo98', 9.238E-04)
    steel.temperature = struct_temp

    # --- Sodium Coolant (Composition 28) ---
    sodium = openmc.Material(name='Sodium')
    sodium.set_density('sum')
    sodium.add_nuclide('Na23', 2.074E-02)
    sodium.temperature = struct_temp

    # --- SSA1 Steel Shielding, 1st row (Composition 31, homogenized) ---
    ssa1_mat = openmc.Material(name='SSA1 Shielding')
    ssa1_mat.set_density('sum')
    ssa1_mat.add_nuclide('Na23', 5.638E-03)
    ssa1_mat.add_nuclide('Fe56', 5.252E-02)
    ssa1_mat.add_nuclide('Cr52', 7.636E-03)
    ssa1_mat.add_nuclide('Ni58', 7.447E-04)
    ssa1_mat.add_nuclide('Mo98', 3.589E-04)
    ssa1_mat.temperature = struct_temp

    # --- SSA2 Steel Shielding, outer rows (Composition 32, homogenized) ---
    ssa2_mat = openmc.Material(name='SSA2 Shielding')
    ssa2_mat.set_density('sum')
    ssa2_mat.add_nuclide('Na23', 5.875E-03)
    ssa2_mat.add_nuclide('Fe56', 5.179E-02)
    ssa2_mat.add_nuclide('Cr52', 7.530E-03)
    ssa2_mat.add_nuclide('Ni58', 7.332E-04)
    ssa2_mat.add_nuclide('Mo98', 3.550E-04)
    ssa2_mat.temperature = struct_temp

    # --- Radial Reflector (Composition 33, homogenized) ---
    rr_mat = openmc.Material(name='Radial Reflector')
    rr_mat.set_density('sum')
    rr_mat.add_nuclide('Na23', 4.860E-03)
    rr_mat.add_nuclide('Fe56', 4.630E-02)
    rr_mat.add_nuclide('Cr52', 1.340E-02)
    rr_mat.add_nuclide('Ni58', 6.280E-03)
    rr_mat.temperature = struct_temp

    # --- Axial Reflector (Composition 34, steel + sodium mix) ---
    axial_refl = openmc.Material(name='Axial Reflector')
    axial_refl.set_density('sum')
    axial_refl.add_nuclide('Fe56', 1.729E-02)
    axial_refl.add_nuclide('Na23', 1.368E-02)
    axial_refl.add_nuclide('Cr52', 2.344E-03)
    axial_refl.add_nuclide('Ni58', 4.948E-05)
    axial_refl.add_nuclide('Mo98', 1.487E-04)
    axial_refl.temperature = struct_temp

    # Fuel material map
    fuel_materials = {
        ZONE_LEZ: fuel_lez,
        ZONE_MEZ: fuel_mez,
        ZONE_MOX: fuel_mox,
        ZONE_HEZ: fuel_hez,
    }

    model.materials = openmc.Materials([
        fuel_lez, fuel_mez, fuel_mox, fuel_hez, fuel_blanket,
        steel, sodium, ssa1_mat, ssa2_mat, rr_mat, axial_refl,
    ])
    model.materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # SURFACES
    # =========================================================================
    fuel_surf = openmc.ZCylinder(r=fuel_or)
    clad_surf = openmc.ZCylinder(r=clad_or)

    z_bot_surf = openmc.ZPlane(z0=z_bot)
    z_lb_bot_surf = openmc.ZPlane(z0=z_lb_bot)
    z_fbot_surf = openmc.ZPlane(z0=z_fuel_bot)
    z_ftop_surf = openmc.ZPlane(z0=z_fuel_top)
    z_ub_top_surf = openmc.ZPlane(z0=z_ub_top)
    z_top_surf = openmc.ZPlane(z0=z_top)

    # =========================================================================
    # PIN-CELL UNIVERSES
    # =========================================================================
    # Pin structure: fuel pellet -> cladding (no gap) -> sodium coolant
    # No central hole per benchmark specification.

    def make_fuel_pin(fuel_mat, zone_name):
        """Create a fuel pin universe with axial blankets and reflectors."""
        # Axial regions
        lower_refl_reg = +z_bot_surf & -z_lb_bot_surf
        lower_blank_reg = +z_lb_bot_surf & -z_fbot_surf
        active_reg = +z_fbot_surf & -z_ftop_surf
        upper_blank_reg = +z_ftop_surf & -z_ub_top_surf
        upper_refl_reg = +z_ub_top_surf & -z_top_surf

        cells = []

        # --- Active fuel region ---
        c = openmc.Cell(name=f'{zone_name} fuel')
        c.region = -fuel_surf & active_reg
        c.fill = fuel_mat
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} clad active')
        c.region = +fuel_surf & -clad_surf & active_reg
        c.fill = steel
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} Na active')
        c.region = +clad_surf & active_reg
        c.fill = sodium
        cells.append(c)

        # --- Lower blanket ---
        c = openmc.Cell(name=f'{zone_name} LB fuel')
        c.region = -fuel_surf & lower_blank_reg
        c.fill = fuel_blanket
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} LB clad')
        c.region = +fuel_surf & -clad_surf & lower_blank_reg
        c.fill = steel
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} LB Na')
        c.region = +clad_surf & lower_blank_reg
        c.fill = sodium
        cells.append(c)

        # --- Upper blanket ---
        c = openmc.Cell(name=f'{zone_name} UB fuel')
        c.region = -fuel_surf & upper_blank_reg
        c.fill = fuel_blanket
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} UB clad')
        c.region = +fuel_surf & -clad_surf & upper_blank_reg
        c.fill = steel
        cells.append(c)

        c = openmc.Cell(name=f'{zone_name} UB Na')
        c.region = +clad_surf & upper_blank_reg
        c.fill = sodium
        cells.append(c)

        # --- Lower axial reflector (homogenized steel+sodium) ---
        c = openmc.Cell(name=f'{zone_name} lower refl')
        c.region = lower_refl_reg
        c.fill = axial_refl
        cells.append(c)

        # --- Upper axial reflector (homogenized steel+sodium) ---
        c = openmc.Cell(name=f'{zone_name} upper refl')
        c.region = upper_refl_reg
        c.fill = axial_refl
        cells.append(c)

        return openmc.Universe(name=f'{zone_name} Pin', cells=cells)

    pin_universes = {}
    for zone_name, fuel_mat in fuel_materials.items():
        pin_universes[zone_name] = make_fuel_pin(fuel_mat, zone_name)

    # =========================================================================
    # ASSEMBLY UNIVERSES
    # =========================================================================

    def make_fuel_assembly(zone_name, pin_univ):
        """Create a fuel assembly with hex pin lattice inside wrapper."""
        # Pin hex lattice
        lat = openmc.HexLattice(name=f'{zone_name} Pins')
        lat.orientation = 'y'
        lat.center = (0.0, 0.0)
        lat.pitch = [pin_pitch]

        # Build universes: outermost ring first in OpenMC convention
        univs = []
        for ring_idx in range(n_pin_rings):
            phys_ring = n_pin_rings - 1 - ring_idx
            if phys_ring == 0:
                univs.append([pin_univ])
            else:
                univs.append([pin_univ] * (6 * phys_ring))
        lat.universes = univs

        # Outer universe for lattice (sodium between pins and wrapper)
        na_outer = openmc.Universe(name=f'{zone_name} lat outer',
                                   cells=[openmc.Cell(fill=sodium)])
        lat.outer = na_outer

        # Wrapper hex prisms
        inner_edge = inner_ftf / math.sqrt(3.0)
        outer_edge = assy_ftf / math.sqrt(3.0)

        inner_hex = openmc.model.HexagonalPrism(
            edge_length=inner_edge, orientation='y', origin=(0., 0.))
        outer_hex = openmc.model.HexagonalPrism(
            edge_length=outer_edge, orientation='y', origin=(0., 0.))

        full_z = +z_bot_surf & -z_top_surf

        # Pin lattice inside wrapper
        c1 = openmc.Cell(name=f'{zone_name} lattice')
        c1.region = -inner_hex & full_z
        c1.fill = lat

        # Wrapper wall
        c2 = openmc.Cell(name=f'{zone_name} wrapper')
        c2.region = +inner_hex & -outer_hex & full_z
        c2.fill = steel

        # Inter-assembly sodium
        c3 = openmc.Cell(name=f'{zone_name} gap')
        c3.region = +outer_hex | ~full_z
        c3.fill = sodium

        return openmc.Universe(name=f'{zone_name} Assy',
                               cells=[c1, c2, c3])

    def make_homog_assembly(zone_name, material):
        """Create a homogenized assembly (SSA, reflector)."""
        outer_edge = assy_ftf / math.sqrt(3.0)
        outer_hex = openmc.model.HexagonalPrism(
            edge_length=outer_edge, orientation='y', origin=(0., 0.))

        full_z = +z_bot_surf & -z_top_surf

        # Homogenized interior
        c1 = openmc.Cell(name=f'{zone_name} interior')
        c1.region = -outer_hex & full_z
        c1.fill = material

        # Inter-assembly sodium (outside wrapper + above/below)
        c2 = openmc.Cell(name=f'{zone_name} gap')
        c2.region = +outer_hex | ~full_z
        c2.fill = sodium

        return openmc.Universe(name=f'{zone_name} Assy',
                               cells=[c1, c2])

    # Map zone materials for non-fuel assemblies
    homog_materials = {
        ZONE_SSA1: ssa1_mat,
        ZONE_SSA2: ssa2_mat,
        ZONE_RR: rr_mat,
    }

    assembly_universes = {}
    for zone_name in FUEL_ZONES:
        assembly_universes[zone_name] = make_fuel_assembly(
            zone_name, pin_universes[zone_name])
    for zone_name, mat in homog_materials.items():
        assembly_universes[zone_name] = make_homog_assembly(zone_name, mat)

    # =========================================================================
    # CORE-LEVEL HEX LATTICE
    # =========================================================================
    core_lattice = openmc.HexLattice(name='BN-600 Core')
    core_lattice.orientation = 'y'
    core_lattice.center = (0.0, 0.0)
    core_lattice.pitch = [assy_pitch]

    core_univs = []
    zone_counts = {}
    for ring_idx in range(NUM_CORE_RINGS):
        phys_ring = NUM_CORE_RINGS - 1 - ring_idx
        zone = RING_ZONES[phys_ring]
        assy_univ = assembly_universes[zone]
        zone_counts[zone] = zone_counts.get(zone, 0)

        if phys_ring == 0:
            core_univs.append([assy_univ])
            zone_counts[zone] += 1
        else:
            n_pos = 6 * phys_ring
            core_univs.append([assy_univ] * n_pos)
            zone_counts[zone] += n_pos

    core_lattice.universes = core_univs

    # Outer universe: sodium beyond outermost ring
    refl_cell = openmc.Cell(name='Sodium reflector outer')
    refl_cell.fill = sodium
    core_lattice.outer = openmc.Universe(name='Outer Na',
                                         cells=[refl_cell])

    n_fuel = sum(v for k, v in zone_counts.items() if k in FUEL_ZONES)

    # =========================================================================
    # ROOT UNIVERSE
    # =========================================================================
    core_radius = (NUM_CORE_RINGS - 0.5) * assy_pitch + 10.0
    radial_bc = openmc.ZCylinder(r=core_radius, boundary_type='vacuum')
    z_lo_bc = openmc.ZPlane(z0=z_bot - 10.0, boundary_type='vacuum')
    z_hi_bc = openmc.ZPlane(z0=z_top + 10.0, boundary_type='vacuum')

    core_cell = openmc.Cell(name='Core')
    core_cell.region = -radial_bc & +z_lo_bc & -z_hi_bc
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

    # Source in active fuel region
    src_r = 9 * assy_pitch
    z_mid = (z_fuel_bot + z_fuel_top) / 2.0
    z_half = active_height / 2.0
    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(
        [-src_r, -src_r, z_mid - z_half],
        [src_r, src_r, z_mid + z_half])
    source.constraints = {'fissionable': True}
    settings.source = source

    settings.temperature = {
        'method': 'interpolation',
        'default': struct_temp,
    }

    model.settings = settings

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("BN-600 Hybrid Core Benchmark Model (IAEA-TECDOC-1623)")
    print("=" * 60)
    for z in [ZONE_LEZ, ZONE_MEZ, ZONE_MOX, ZONE_HEZ,
              ZONE_SSA1, ZONE_SSA2, ZONE_RR]:
        print(f"  {z:6s}: {zone_counts.get(z, 0):4d} assemblies")
    print(f"  Total fuel assemblies:  {n_fuel}")
    print(f"  Pin pitch:              {pin_pitch} cm")
    print(f"  Assembly pitch:         {assy_pitch} cm")
    print(f"  Active height:          {active_height} cm")
    print(f"  Lower/Upper blanket:    {lower_blanket_h}/{upper_blanket_h} cm")
    print(f"  Fuel temperature:       {fuel_temp} K")
    print(f"  Structure temperature:  {struct_temp} K")
    print(f"  Particles/batch:        {particles}")
    print(f"  Batches (total/inact):  {batches}/{inactive}")
    print("=" * 60)

    return model


def main():
    parser = argparse.ArgumentParser(
        description='BN-600 Hybrid Core Benchmark (IAEA-TECDOC-1623)')
    parser.add_argument('--particles', type=int, default=15000)
    parser.add_argument('--batches', type=int, default=200)
    parser.add_argument('--inactive', type=int, default=50)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()

    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
    )
    model.export_to_model_xml()
    print("\nModel exported to XML files")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
