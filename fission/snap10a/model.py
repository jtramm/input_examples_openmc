#!/usr/bin/env python3
"""
=============================================================================
SNAP-10A/2 Space Reactor Benchmark - SCA-4B Case 8490a
=============================================================================

Source: ORNL/TM-2005/54, Volume I
        "Experimental Criticality Benchmarks for SNAP 10A/2 Reactor Cores"
        A. W. Krass and K. L. Goluoglu, April 2005

Overview
--------
The SNAP-10A (Systems for Nuclear Auxiliary Power) was a compact space reactor
developed by Atomics International (AI) in the 1960s. It was the first (and
only) nuclear fission reactor to be operated in space, launched aboard the
SNAPSHOT satellite on April 3, 1965. The reactor produced 500 watts of
electrical power using thermoelectric converters.

The SCA-4B (SNAP Critical Assembly 4B) experimental program tested critical
configurations of SNAP 10A/2 cores at the AI facility. These experiments
evaluated water immersion scenarios for criticality safety analysis.

This Model: Case 8490a
-----------------------
This is the simplest benchmark configuration from Phase I of the SCA-4B
program. It consists of:
  - 28 SCA-4 fuel elements in the first 28 lattice positions
  - No Lucite rods (vacant positions 29-37 are water-filled)
  - No internal beryllium reflector inserts
  - Core vessel flooded with water
  - Full water reflection via surrounding water tanks
  - Upper tank water level: 9.0 inches (22.86 cm)

The assembly was subcritical by approximately 0.02 dollars.

Reference MCNP5 Result
-----------------------
  k_exp  = 0.99984 (estimated experimental, -0.02$ subcritical)
  k_calc = 1.00491 +/- 0.00069 (MCNP5, ENDF/B-VI)
  k_calc/k_exp = 1.0051 (normalized)

Fuel Description
----------------
The fuel is a uranium-zirconium hydride (U-ZrH) consisting of approximately
10 wt% uranium (enriched to >= 93 wt% U-235) and 90 wt% zirconium hydride.
The hydrogen content provides moderation within the fuel itself. Each fuel
element contains approximately 128.5 g of U-235. The fuel rod is clad in
Hastelloy N (a nickel-molybdenum alloy). The gap between fuel and cladding
is filled with a hydrogen barrier coating containing samarium oxide (Sm2O3)
as a burnable poison.

Core Geometry
-------------
The reactor core vessel is a cylindrical shell of 316 stainless steel,
closed at the bottom with a removable aluminum cover plate on top. The
vessel holds up to 37 fuel elements on a triangular pitch of 1.260 inches
(3.2004 cm). For water reflection, the vessel is mounted inside nested
cylindrical stainless steel water tanks (upper tank, lower tank, and a
control cap tank below).

All dimensions are from ORNL/TM-2005/54, Appendix C.
=============================================================================
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments for number of particles, batches, and inactive
# =============================================================================
parser = argparse.ArgumentParser(
    description="SNAP-10A/2 SCA-4B Case 8490a benchmark model"
)
parser.add_argument(
    "--particles", type=int, default=10000,
    help="Number of neutron histories per batch (default: 10000)"
)
parser.add_argument(
    "--batches", type=int, default=150,
    help="Total number of batches (default: 150)"
)
parser.add_argument(
    "--inactive", type=int, default=50,
    help="Number of inactive (discarded) batches (default: 50)"
)
args = parser.parse_args()


# =============================================================================
# MATERIALS
# =============================================================================
# All number densities taken directly from ORNL/TM-2005/54, Appendix B.
# Units: atom/barn-cm

# ---- Material 1: U-ZrH fuel (SCA-4 average) ----
# Density: 6.036 g/cm3 (weighted average of all 37 fuel rod assemblies)
# Composition: ~10 wt% U (93.3% enriched), ~90 wt% Zr hydride
# Hydrogen content: 6.40e22 H atoms/cm3 (revised downward from 6.483e22)
# Total mass per 37 rods: 51,722.59 g (4,760.38 g U-235, 343.64 g U-238,
#   917.79 g H, 45,700.78 g Zr)
# NTotal = 0.1007323 atom/barn-cm
fuel = openmc.Material(name="U-ZrH fuel (SCA-4)")
fuel.set_density("atom/b-cm", 0.1007323)
fuel.add_nuclide("U235", 0.0014234)   # 235-U from enriched uranium
fuel.add_nuclide("U238", 0.0001015)   # 238-U remainder
fuel.add_nuclide("Zr90", 0.0352074 * 0.5145)   # Natural Zr isotopic fractions
fuel.add_nuclide("Zr91", 0.0352074 * 0.1122)
fuel.add_nuclide("Zr92", 0.0352074 * 0.1715)
fuel.add_nuclide("Zr94", 0.0352074 * 0.1738)
fuel.add_nuclide("Zr96", 0.0352074 * 0.0280)
fuel.add_nuclide("H1", 0.0640000)     # Hydrogen in ZrH
# S(alpha,beta) thermal scattering for H in ZrH and Zr in ZrH
fuel.add_s_alpha_beta("c_H_in_ZrH")
fuel.add_s_alpha_beta("c_Zr_in_ZrH")

# ---- Material 5: Hydrogen barrier coating (Sm2O3 + Al2O3 + SiO2) ----
# This coating fills the radial gap between the fuel rod and cladding tube.
# It serves as both a hydrogen diffusion barrier and contains samarium oxide
# as a burnable neutron absorber. The SCA-4 fuel uses the ca.1961 AI spec
# of 6.6 mg Sm2O3/inch of fuel.
# Total coat weight for 37 elements: 220.57 g (3.00 g Sm2O3, rest is
# simplified as 50/50 Al2O3 and SiO2 by weight)
# NTotal = 0.0254232 atom/barn-cm
coating = openmc.Material(name="H-barrier coating (Sm2O3/Al2O3/SiO2)")
coating.set_density("atom/b-cm", 0.0254232)
coating.add_nuclide("O16", 0.0161151)     # Oxygen from all three oxides
coating.add_nuclide("Al27", 0.0050217)    # Aluminum from Al2O3
coating.add_nuclide("Si28", 0.0042609)    # Silicon from SiO2 (treated as nat Si)
# Samarium isotopes (only those with significant absorption cross sections):
# 144Sm and 148Sm omitted per ORNL report (low absorption, no MCNP5 data)
coating.add_nuclide("Sm147", 0.0000061)   # 15.0 at% of natural Sm
coating.add_nuclide("Sm149", 0.0000056)   # 13.8 at% (highest absorber)
coating.add_nuclide("Sm150", 0.0000030)   # 7.4 at%
coating.add_nuclide("Sm152", 0.0000108)   # 26.7 at%

# ---- Material 6: Type 1100 Aluminum (grid plates and cover plate) ----
# Density: 2.7 g/cm3
# Composition: 100% Al by weight
# NTotal = 0.0602626 atom/barn-cm
aluminum = openmc.Material(name="Type 1100 Aluminum")
aluminum.set_density("atom/b-cm", 0.0602626)
aluminum.add_nuclide("Al27", 0.0602626)

# ---- Material 7: Light water (moderator and reflector) ----
# Density: 0.9982 g/cm3 at room temperature
# NTotal = 0.1001037 atom/barn-cm
water = openmc.Material(name="Light water")
water.set_density("atom/b-cm", 0.1001037)
water.add_nuclide("H1", 0.0667358)
water.add_nuclide("O16", 0.0333679)
water.add_s_alpha_beta("c_H_in_H2O")  # Thermal scattering for H in water

# ---- Material 8: Hastelloy N (fuel element cladding) ----
# Density: 8.86 g/cm3
# Composition: 71 wt% Ni, 7 wt% Cr, 16 wt% Mo, 5 wt% Fe, 1 wt% Si
# This is a nickel-molybdenum alloy with excellent high-temperature
# corrosion resistance, developed specifically for molten salt and
# high-temperature nuclear applications.
# NTotal = 0.0872946 atom/barn-cm
hastelloy_n = openmc.Material(name="Hastelloy N")
hastelloy_n.set_density("atom/b-cm", 0.0872946)
# Nickel isotopes (71 wt% Ni)
hastelloy_n.add_nuclide("Ni58", 0.0440590)
hastelloy_n.add_nuclide("Ni60", 0.0168440)
hastelloy_n.add_nuclide("Ni61", 0.0007293)
hastelloy_n.add_nuclide("Ni62", 0.0023169)
hastelloy_n.add_nuclide("Ni64", 0.0005873)
# Chromium isotopes (7 wt% Cr)
hastelloy_n.add_nuclide("Cr50", 0.0003121)
hastelloy_n.add_nuclide("Cr52", 0.0060187)
hastelloy_n.add_nuclide("Cr53", 0.0006824)
hastelloy_n.add_nuclide("Cr54", 0.0001699)
# Molybdenum (16 wt% Mo) - using natural Mo
hastelloy_n.add_element("Mo", 0.0088982)
# Iron isotopes (5 wt% Fe)
hastelloy_n.add_nuclide("Fe54", 0.0002818)
hastelloy_n.add_nuclide("Fe56", 0.0043815)
hastelloy_n.add_nuclide("Fe57", 0.0001003)
hastelloy_n.add_nuclide("Fe58", 0.0000134)
# Silicon (1 wt% Si) - using natural Si
hastelloy_n.add_element("Si", 0.0018998)

# ---- Material 9: SS316 stainless steel (core vessel, tank walls) ----
# Density: 8.03 g/cm3 (calculated from isotopic distribution)
# Composition: 65.5 wt% Fe, 17 wt% Cr, 12 wt% Ni, 2.5 wt% Mo,
#              2 wt% Mn, 1 wt% Si
# NTotal = 0.0871549 atom/barn-cm
ss316 = openmc.Material(name="SS316 stainless steel")
ss316.set_density("atom/b-cm", 0.0871549)
# Iron isotopes (65.5 wt% Fe)
ss316.add_nuclide("Fe54", 0.0033463)
ss316.add_nuclide("Fe56", 0.0520202)
ss316.add_nuclide("Fe57", 0.0011910)
ss316.add_nuclide("Fe58", 0.0001588)
# Chromium isotopes (17 wt% Cr)
ss316.add_nuclide("Cr50", 0.0006870)
ss316.add_nuclide("Cr52", 0.0132476)
ss316.add_nuclide("Cr53", 0.0015020)
ss316.add_nuclide("Cr54", 0.0003739)
# Nickel isotopes (12 wt% Ni)
ss316.add_nuclide("Ni58", 0.0067490)
ss316.add_nuclide("Ni60", 0.0025802)
ss316.add_nuclide("Ni61", 0.0001117)
ss316.add_nuclide("Ni62", 0.0003549)
ss316.add_nuclide("Ni64", 0.0000900)
# Molybdenum (2.5 wt% Mo)
ss316.add_element("Mo", 0.0012601)
# Manganese (2 wt% Mn)
ss316.add_nuclide("Mn55", 0.0017604)
# Silicon (1 wt% Si)
ss316.add_element("Si", 0.0017218)

# Collect all materials
materials = openmc.Materials([fuel, coating, aluminum, water, hastelloy_n, ss316])
materials.export_to_xml()


# =============================================================================
# GEOMETRY
# =============================================================================
# All dimensions from ORNL/TM-2005/54, Appendix C.
# The coordinate system is centered on the midplane of the fuel elements,
# with the z-axis along the vertical axis of the core vessel.

# ---- Fuel element dimensions ----
# From Appendix C, Section C-1.1 (Fig. C-1):
#   Fuel rod diameter: 1.212 in. (3.07848 cm) -> radius = 1.53924 cm
#   Cladding tube ID: 1.230 in. (3.1242 cm) -> IR = 1.5621 cm
#     (= OD - 2*wall = 1.250 - 2*0.010 = 1.230 in.)
#   Cladding tube OD: 1.250 in. (3.175 cm) -> OR = 1.5875 cm
#   Fuel rod length: 12.225 in. (31.0515 cm) -> half-length = 15.52575 cm
#   Fuel element length: 12.450 in. (31.623 cm) -> half-length = 15.8115 cm
#   Wall thickness: 0.010 in. (0.0254 cm)
#   Radial gap: (1.5621 - 1.53924) = 0.02286 cm (filled with H-barrier coating)

fuel_rod_radius = 1.53924       # cm, radius of UZrH fuel meat
clad_inner_radius = 1.5621      # cm, inner radius of Hastelloy N tube
clad_outer_radius = 1.5875      # cm, outer radius of Hastelloy N tube
fuel_half_length = 15.52575     # cm, half-length of active fuel region
element_half_length = 15.8115   # cm, half-length of total fuel element (incl. end caps)

# ---- Core vessel dimensions ----
# From Appendix C, Section C-2 (Fig. C-4):
#   Vessel ID: 8.900 in. (22.606 cm) -> IR = 11.303 cm
#   Vessel wall thickness: 0.032 in. at sides
#   Vessel OD: 8.965 in. + 0.062 in. wall -> OR = 11.38555 cm
#     Note: OR = (8.965/2 + 0.031) * 2.54 = 11.38555 cm
#   Flange OR: (8.965/2 + 2) * 2.54 = 16.46555 cm
vessel_ir = 11.303              # cm, inner radius of core vessel
vessel_or = 11.38555            # cm, outer radius of core vessel
flange_or = 16.46555            # cm, outer radius of vessel flange

# ---- Grid plate dimensions ----
# Upper grid plate: 0.375 in. thick aluminum
# Lower grid plate: 0.312 in. thick aluminum
# Both have diameter matching vessel ID (11.303 cm radius)
upper_grid_top = 16.764         # cm, top of upper grid plate
upper_grid_bot = 15.8115        # cm, bottom of upper grid plate (= top of fuel elements)
lower_grid_top = -15.8115       # cm, top of lower grid plate (= bottom of fuel elements)
lower_grid_bot = -16.60652      # cm, bottom of lower grid plate

# ---- Vessel vertical extents ----
# From Fig. C-4 of Appendix C:
cover_flange_top = 18.71726     # cm, top of cover plate flange
cover_indent_top = 18.161       # cm, top of cover plate indent
cover_interface = 18.08226      # cm, vessel cover/flange interface
flange_bottom = 17.92478        # cm, bottom of vessel flange
vessel_bottom_top = -16.764     # cm, top of vessel bottom plate (= bottom of lower grid)
vessel_bottom_bot = -16.84274   # cm, bottom of vessel bottom plate

# ---- Upper water tank ----
# Provides radial reflection over the upper 1/3 of the core vessel plus
# up to 9 inches of axial reflection above the vessel.
# Inner radius: vessel OR + 6 inches radial water gap
upper_tank_ir = 26.62555        # cm
upper_tank_or = 26.78303        # cm, (wall thickness 0.062 in.)
upper_tank_floor_bot = 6.01726  # cm, bottom of upper tank floor plate
upper_tank_floor_top = 6.17474  # cm, top of upper tank floor plate
upper_tank_water_top = 29.03474 # cm, water level at 9 inches in upper tank
upper_tank_top = 31.57474       # cm, top edge of upper tank

# ---- Lower water tank ----
# Provides radial reflection over the lower 2/3 of the core vessel.
# Inner radius: vessel OR + 6.5 inches radial water gap
lower_tank_ir = 27.89555        # cm
lower_tank_or = 28.05303        # cm, (wall thickness 0.062 in.)
lower_tank_floor_top = -17.95526  # cm, top of lower tank floor plate
lower_tank_floor_bot = -18.11274  # cm, bottom of lower tank floor plate

# ---- Control cap tank ----
# Small tank below the core providing axial reflection from below.
# Inner radius is the average of the lower tank IR and vessel OR.
control_cap_ir = 19.64055       # cm
control_cap_or = 19.79803       # cm
control_cap_top_plate_top = -18.36674     # cm
control_cap_top_plate_bot = -18.52422     # cm
control_cap_floor_top = -33.76422         # cm
control_cap_floor_bot = -33.9217          # cm

# ---- Outer boundary ----
# Spherical boundary to contain the entire geometry
outer_radius = 50.0             # cm


# =============================================================================
# FUEL ELEMENT POSITIONS
# =============================================================================
# The 37 lattice positions are arranged in a triangular pitch of 1.260 in.
# (3.2004 cm). For case 8490a, positions 1-28 contain fuel elements and
# positions 29-37 are empty (filled with water since the core is flooded).
#
# Position coordinates from ORNL/TM-2005/54, Appendix C, Table C-1.
# Position 1 is at the center-top of the array. The MCNP5 input places
# position 1 at (x=0, y=9.6012). All other positions are relative to
# position 1. The absolute coordinates are listed below.
#
# The array forms a roughly hexagonal pattern characteristic of the
# triangular pitch arrangement within the cylindrical vessel.

# Absolute (x, y) coordinates of all 37 fuel element positions (cm)
# Position 1 is at (0, 9.6012); others offset from Table C-1
fuel_positions = {
    # Position: (x_cm, y_cm)
    # --- Ring 1 (center of array, 1 element) ---
    1:  ( 0.0000,   9.6012),
    # --- Ring 2 (6 elements around center) ---
    2:  ( 0.0000,   6.4008),   # directly below pos 1
    3:  (-2.7716,   8.0010),   # lower-left of pos 1
    4:  ( 2.7716,   8.0010),   # lower-right of pos 1
    5:  ( 0.0000,   3.2004),
    6:  (-2.7716,   4.8006),
    7:  (-5.5432,   6.4008),
    8:  ( 2.7716,   4.8006),
    9:  ( 5.5432,   6.4008),
    # --- Ring 3 (expanding outward) ---
    10: ( 0.0000,   0.0000),
    11: (-2.7716,   1.6002),
    12: (-5.5432,   3.2004),
    13: (-8.3149,   4.8006),
    14: ( 2.7716,   1.6002),
    15: ( 5.5432,   3.2004),
    16: ( 8.3149,   4.8006),
    # --- Ring 4 ---
    17: ( 0.0000,  -3.2004),
    18: (-2.7716,  -1.6002),
    19: (-5.5432,   0.0000),
    20: (-8.3249,   1.6002),   # Note: x=-8.3249 per MCNP input (slight asymmetry)
    21: ( 2.7716,  -1.6002),
    22: ( 5.5432,   0.0000),
    23: ( 8.3149,   1.6002),
    # --- Ring 5 ---
    24: ( 0.0000,  -6.4008),
    25: (-2.7716,  -4.8006),
    26: (-5.5432,  -3.2004),
    27: (-8.3149,  -1.6002),
    28: ( 2.7716,  -4.8006),
    # --- Positions 29-37 (empty/water-filled for case 8490a) ---
    29: ( 5.5432,  -3.2004),
    30: ( 8.3149,  -1.6002),
    31: ( 0.0000,  -9.6012),
    32: (-2.7716,  -8.0010),
    33: (-5.5432,  -6.4008),   # Note: positions 33-37 from MCNP input
    34: (-8.3149,  -4.8006),
    35: ( 2.7716,  -8.0010),
    36: ( 5.5432,  -6.4008),
    37: ( 8.3149,  -4.8006),
}


# =============================================================================
# SURFACES AND CELLS
# =============================================================================

# ---- Planar surfaces for axial boundaries ----
# These define the vertical extents of fuel, grid plates, vessel, and tanks.
fuel_top_plane = openmc.ZPlane(z0=fuel_half_length)          # Top of active fuel
fuel_bot_plane = openmc.ZPlane(z0=-fuel_half_length)         # Bottom of active fuel
elem_top_plane = openmc.ZPlane(z0=element_half_length)       # Top of fuel element
elem_bot_plane = openmc.ZPlane(z0=-element_half_length)      # Bottom of fuel element
ugrid_top_plane = openmc.ZPlane(z0=upper_grid_top)           # Top of upper grid plate
lgrid_bot_plane = openmc.ZPlane(z0=lower_grid_bot)           # Bottom of lower grid plate
cover_top_plane = openmc.ZPlane(z0=cover_flange_top)         # Top of cover plate flange
cover_indent_plane = openmc.ZPlane(z0=cover_indent_top)      # Top of cover plate indent
cover_iface_plane = openmc.ZPlane(z0=cover_interface)        # Cover/flange interface
flange_bot_plane = openmc.ZPlane(z0=flange_bottom)           # Bottom of vessel flange
vessel_bot_top_plane = openmc.ZPlane(z0=vessel_bottom_top)   # Top of vessel bottom
vessel_bot_bot_plane = openmc.ZPlane(z0=vessel_bottom_bot)   # Bottom of vessel bottom

# Upper tank planes
utank_floor_bot_plane = openmc.ZPlane(z0=upper_tank_floor_bot)
utank_floor_top_plane = openmc.ZPlane(z0=upper_tank_floor_top)
utank_water_plane = openmc.ZPlane(z0=upper_tank_water_top)   # Water level @ 9 in.
utank_top_plane = openmc.ZPlane(z0=upper_tank_top)           # Top edge of tank

# Lower tank planes
ltank_floor_top_plane = openmc.ZPlane(z0=lower_tank_floor_top)
ltank_floor_bot_plane = openmc.ZPlane(z0=lower_tank_floor_bot)

# Control cap tank planes
ccap_top_plate_top_plane = openmc.ZPlane(z0=control_cap_top_plate_top)
ccap_top_plate_bot_plane = openmc.ZPlane(z0=control_cap_top_plate_bot)
ccap_floor_top_plane = openmc.ZPlane(z0=control_cap_floor_top)
ccap_floor_bot_plane = openmc.ZPlane(z0=control_cap_floor_bot)

# ---- Cylindrical surfaces ----
vessel_inner_cyl = openmc.ZCylinder(r=vessel_ir)             # Inner wall of core vessel
vessel_outer_cyl = openmc.ZCylinder(r=vessel_or)             # Outer wall of core vessel
flange_cyl = openmc.ZCylinder(r=flange_or)                   # Vessel flange outer radius
utank_inner_cyl = openmc.ZCylinder(r=upper_tank_ir)          # Upper tank inner wall
utank_outer_cyl = openmc.ZCylinder(r=upper_tank_or)          # Upper tank outer wall
ltank_inner_cyl = openmc.ZCylinder(r=lower_tank_ir)          # Lower tank inner wall
ltank_outer_cyl = openmc.ZCylinder(r=lower_tank_or)          # Lower tank outer wall
ccap_inner_cyl = openmc.ZCylinder(r=control_cap_ir)          # Control cap inner wall
ccap_outer_cyl = openmc.ZCylinder(r=control_cap_or)          # Control cap outer wall

# ---- Outer boundary sphere ----
outer_sphere = openmc.Sphere(r=outer_radius, boundary_type="vacuum")


# =============================================================================
# BUILD FUEL ELEMENT UNIVERSES
# =============================================================================
# Each fuel element is modeled as a set of concentric cylinders:
#   1. Inner fuel rod (U-ZrH) - radius 1.53924 cm
#   2. Hydrogen barrier coating (gap region) - from 1.53924 to 1.5621 cm
#   3. Hastelloy N cladding tube - from 1.5621 to 1.5875 cm
#
# The fuel rod extends axially from -15.52575 to +15.52575 cm.
# The cladding tube (with end caps) extends from -15.8115 to +15.8115 cm.
# End caps are solid Hastelloy N above/below the fuel rod.

# Create a single fuel element universe that will be placed at each position
fuel_element_universe = openmc.Universe(name="Fuel element")

# Fuel rod cylinder (centered at origin of the universe)
fuel_rod_cyl = openmc.ZCylinder(r=fuel_rod_radius)
# Inner cladding surface
clad_inner_cyl = openmc.ZCylinder(r=clad_inner_radius)
# Outer cladding surface
clad_outer_cyl = openmc.ZCylinder(r=clad_outer_radius)
# Axial extent of fuel rod
fuel_top = openmc.ZPlane(z0=fuel_half_length)
fuel_bot = openmc.ZPlane(z0=-fuel_half_length)
# Axial extent of entire element (including end caps)
elem_top = openmc.ZPlane(z0=element_half_length)
elem_bot = openmc.ZPlane(z0=-element_half_length)

# Cell 1: U-ZrH fuel meat (solid cylindrical rod)
fuel_cell = openmc.Cell(name="U-ZrH fuel rod")
fuel_cell.fill = fuel
fuel_cell.region = -fuel_rod_cyl & -fuel_top & +fuel_bot

# Cell 2: Hydrogen barrier/absorber coating (fills radial gap)
# This annular region between the fuel rod and cladding inner wall contains
# the Sm2O3/Al2O3/SiO2 coating mixture
gap_cell = openmc.Cell(name="H-barrier coating (radial gap)")
gap_cell.fill = coating
gap_cell.region = +fuel_rod_cyl & -clad_inner_cyl & -fuel_top & +fuel_bot

# Cell 3: Hastelloy N cladding tube and end caps
# The cladding is the annular tube plus solid end caps above/below fuel
clad_cell = openmc.Cell(name="Hastelloy N cladding")
clad_cell.fill = hastelloy_n
clad_cell.region = (+clad_inner_cyl & -clad_outer_cyl & -elem_top & +elem_bot) | \
                   (-clad_outer_cyl & (+fuel_top | -fuel_bot) & -elem_top & +elem_bot & -clad_inner_cyl) | \
                   (-fuel_rod_cyl & (+fuel_top | -fuel_bot) & -elem_top & +elem_bot)
# Simplified: everything inside the cladding OD, within element length,
# that is NOT the fuel or gap region
clad_cell.region = (-clad_outer_cyl & -elem_top & +elem_bot) & \
                   ~(-fuel_rod_cyl & -fuel_top & +fuel_bot) & \
                   ~(+fuel_rod_cyl & -clad_inner_cyl & -fuel_top & +fuel_bot)

# Cell 4: Water outside the cladding tube (within the universe bounding box)
# This is the water between fuel elements inside the vessel
water_around_cell = openmc.Cell(name="Water around fuel element")
water_around_cell.fill = water
water_around_cell.region = +clad_outer_cyl | +elem_top | -elem_bot

fuel_element_universe.add_cells([fuel_cell, gap_cell, clad_cell, water_around_cell])


# =============================================================================
# BUILD THE MAIN GEOMETRY (ROOT UNIVERSE)
# =============================================================================
root_universe = openmc.Universe(name="Root universe")
all_cells = []

# ---- Place fuel elements at positions 1-28 ----
# Each fuel element is placed as a cell filled with the fuel element universe,
# bounded by the cladding outer radius cylinder centered at the element's
# (x, y) position. The MCNP5 model uses the TRCL card for translations.
fuel_element_cells = []
for pos_num in range(1, 29):
    x, y = fuel_positions[pos_num]
    # Create a bounding cylinder for this fuel element at its (x,y) position
    bound_cyl = openmc.ZCylinder(x0=x, y0=y, r=clad_outer_radius)
    # Create a cell that contains the fuel element universe
    cell = openmc.Cell(
        name=f"Fuel element at position {pos_num}"
    )
    cell.fill = fuel_element_universe
    # The cell region is inside the bounding cylinder, within element axial extent
    cell.region = -bound_cyl & -elem_top_plane & +elem_bot_plane
    cell.translation = (x, y, 0.0)
    # Actually, we need to use the fuel_element_universe centered at origin
    # and translate it. The bounding region should use the translated cylinder.
    fuel_element_cells.append(cell)
    all_cells.append(cell)

# ---- Water inside vessel (between and around fuel elements) ----
# This is the water that fills the flooded core vessel interior, excluding
# the fuel element volumes. The region extends radially to the vessel ID
# and axially between the grid plates.
vessel_water_region = -vessel_inner_cyl & -elem_top_plane & +elem_bot_plane
for cell in fuel_element_cells:
    # Exclude each fuel element's bounding cylinder
    pos_num = fuel_element_cells.index(cell) + 1
    x, y = fuel_positions[pos_num]
    bound_cyl = openmc.ZCylinder(x0=x, y0=y, r=clad_outer_radius)
    vessel_water_region = vessel_water_region & +bound_cyl

vessel_water_cell = openmc.Cell(name="Water inside vessel (between fuel elements)")
vessel_water_cell.fill = water
vessel_water_cell.region = vessel_water_region
all_cells.append(vessel_water_cell)

# ---- Upper aluminum grid plate ----
# A solid aluminum disk at the top of the fuel element region.
# Diameter matches vessel ID (11.303 cm radius).
# Thickness: 0.375 in. (0.9525 cm), from elem_top to ugrid_top
upper_grid_cell = openmc.Cell(name="Upper aluminum grid plate")
upper_grid_cell.fill = aluminum
upper_grid_cell.region = -vessel_inner_cyl & -ugrid_top_plane & +elem_top_plane
all_cells.append(upper_grid_cell)

# ---- Lower aluminum grid plate ----
# Similar aluminum disk at the bottom of the fuel element region.
# Thickness: 0.312 in. (0.79248 cm), from lgrid_bot to elem_bot
lower_grid_cell = openmc.Cell(name="Lower aluminum grid plate")
lower_grid_cell.fill = aluminum
lower_grid_cell.region = -vessel_inner_cyl & -elem_bot_plane & +lgrid_bot_plane
all_cells.append(lower_grid_cell)

# ---- Water-filled vessel interior above upper grid plate ----
# Water between upper grid plate and the vessel cover/flange
water_above_grid = openmc.Cell(name="Water above upper grid plate (inside vessel)")
water_above_grid.fill = water
water_above_grid.region = -vessel_inner_cyl & -cover_iface_plane & +ugrid_top_plane
all_cells.append(water_above_grid)

# ---- Water-filled vessel interior below lower grid plate ----
# Water between lower grid plate and vessel bottom
water_below_grid = openmc.Cell(name="Water below lower grid plate (inside vessel)")
water_below_grid.fill = water
water_below_grid.region = -vessel_inner_cyl & -lgrid_bot_plane & +vessel_bot_top_plane
all_cells.append(water_below_grid)

# ---- SS316 core vessel walls ----
# The vessel consists of:
#   (a) Cylindrical side wall (vessel_ir to vessel_or)
#   (b) Bottom plate (vessel_bot_top to vessel_bot_bot)
#   (c) Flange and cover plate assembly at top
vessel_side_region = (+vessel_inner_cyl & -vessel_outer_cyl &
                      -flange_bot_plane & +vessel_bot_bot_plane)
vessel_bottom_region = (-vessel_outer_cyl & -vessel_bot_top_plane & +vessel_bot_bot_plane)
vessel_flange_region = (-flange_cyl & +vessel_inner_cyl &
                        -cover_iface_plane & +flange_bot_plane)
vessel_cover_region = ((-vessel_inner_cyl & -cover_indent_plane & +cover_iface_plane) |
                       (-flange_cyl & -cover_top_plane & +cover_iface_plane &
                        +cover_indent_plane))

# Combine all vessel steel regions
vessel_cell = openmc.Cell(name="SS316 core vessel (walls, bottom, flange)")
vessel_cell.fill = ss316
vessel_cell.region = vessel_side_region | vessel_bottom_region
all_cells.append(vessel_cell)

# ---- Aluminum cover plate ----
# The cover plate sits on the flange and seals the top of the vessel.
# It consists of a flat disk inside vessel IR + the flange lip.
cover_cell = openmc.Cell(name="Aluminum cover plate and flange")
cover_cell.fill = aluminum
cover_cell.region = ((-vessel_inner_cyl & -cover_top_plane & +cover_iface_plane) |
                     (-flange_cyl & +vessel_inner_cyl & -cover_top_plane & +flange_bot_plane))
all_cells.append(cover_cell)

# ---- Upper water tank ----
# The upper tank provides radial water reflection over the upper portion
# of the core vessel plus axial reflection above.
# Water-filled volume:
utank_water_cell = openmc.Cell(name="Water in upper tank")
utank_water_cell.fill = water
utank_water_cell.region = (
    # Radial water annulus between vessel and tank wall, above tank floor
    (-utank_inner_cyl & +vessel_outer_cyl & -utank_water_plane & +utank_floor_top_plane &
     -cover_top_plane) |
    # Water above the cover plate, inside tank
    (-utank_inner_cyl & -utank_water_plane & +cover_top_plane) |
    # Water above flange level, outside vessel but inside tank
    (-utank_inner_cyl & +flange_cyl & -cover_top_plane & +flange_bot_plane &
     -utank_water_plane & +utank_floor_top_plane)
)
all_cells.append(utank_water_cell)

# Void above water level in upper tank
utank_void_cell = openmc.Cell(name="Void above water in upper tank")
utank_void_cell.region = -utank_inner_cyl & -utank_top_plane & +utank_water_plane
all_cells.append(utank_void_cell)

# SS316 upper tank walls
utank_walls_cell = openmc.Cell(name="SS316 upper tank walls")
utank_walls_cell.fill = ss316
utank_walls_cell.region = (
    (+utank_inner_cyl & -utank_outer_cyl & -utank_top_plane & +utank_floor_bot_plane) |
    (-utank_outer_cyl & +utank_inner_cyl & -utank_floor_top_plane & +utank_floor_bot_plane) |
    (-utank_inner_cyl & -utank_floor_top_plane & +utank_floor_bot_plane &
     +vessel_outer_cyl)
)
# Simplify: side wall + floor plate
utank_walls_cell.region = (
    # Side walls
    (+utank_inner_cyl & -utank_outer_cyl & -utank_top_plane & +utank_floor_top_plane) |
    # Floor plate (annular, between vessel OD and tank ID)
    (-utank_inner_cyl & +vessel_outer_cyl & -utank_floor_top_plane & +utank_floor_bot_plane)
)
all_cells.append(utank_walls_cell)

# ---- Lower water tank ----
# Provides radial reflection over the lower 2/3 of the vessel.
ltank_water_cell = openmc.Cell(name="Water in lower tank")
ltank_water_cell.fill = water
ltank_water_cell.region = (
    (-ltank_inner_cyl & +vessel_outer_cyl &
     -utank_floor_bot_plane & +ltank_floor_top_plane) &
    ~(-utank_inner_cyl & +vessel_outer_cyl &
      -utank_floor_top_plane & +utank_floor_bot_plane)
)
# Simpler: water in the annulus around the vessel, in the lower tank region
# that doesn't overlap with the upper tank floor
ltank_water_cell.region = (
    # Water in lower tank annulus (below upper tank floor)
    (-ltank_inner_cyl & +vessel_outer_cyl &
     -utank_floor_bot_plane & +vessel_bot_bot_plane) &
    # Exclude the upper tank floor plate region
    ~(-utank_inner_cyl & +vessel_outer_cyl &
      -utank_floor_top_plane & +utank_floor_bot_plane) &
    # Exclude vessel wall itself (already defined)
    ~(+vessel_inner_cyl & -vessel_outer_cyl &
      -flange_bot_plane & +vessel_bot_bot_plane) &
    # Exclude vessel bottom plate
    ~(-vessel_outer_cyl & -vessel_bot_top_plane & +vessel_bot_bot_plane)
)
all_cells.append(ltank_water_cell)

# SS316 lower tank walls
ltank_walls_cell = openmc.Cell(name="SS316 lower tank walls")
ltank_walls_cell.fill = ss316
ltank_walls_cell.region = (
    # Side walls
    (+ltank_inner_cyl & -ltank_outer_cyl &
     -utank_floor_bot_plane & +ltank_floor_top_plane) |
    # Floor plate
    (-ltank_inner_cyl & -ltank_floor_top_plane & +ltank_floor_bot_plane)
)
all_cells.append(ltank_walls_cell)

# ---- Control cap tank ----
# Small tank below the core providing axial reflection from below.
ccap_water_cell = openmc.Cell(name="Water in control cap tank")
ccap_water_cell.fill = water
ccap_water_cell.region = (-ccap_inner_cyl &
                          -ccap_top_plate_bot_plane & +ccap_floor_top_plane)
all_cells.append(ccap_water_cell)

# Aluminum top plate of control cap
ccap_top_plate_cell = openmc.Cell(name="Aluminum top plate of control cap tank")
ccap_top_plate_cell.fill = aluminum
ccap_top_plate_cell.region = (-ccap_inner_cyl &
                              -ccap_top_plate_top_plane & +ccap_top_plate_bot_plane)
all_cells.append(ccap_top_plate_cell)

# SS316 control cap tank walls
ccap_walls_cell = openmc.Cell(name="SS316 control cap tank walls")
ccap_walls_cell.fill = ss316
ccap_walls_cell.region = (
    # Side walls
    (+ccap_inner_cyl & -ccap_outer_cyl &
     -ccap_top_plate_top_plane & +ccap_floor_top_plane) |
    # Floor plate
    (-ccap_inner_cyl & -ccap_floor_top_plane & +ccap_floor_bot_plane)
)
all_cells.append(ccap_walls_cell)

# ---- Void between tanks and vessel where no water ----
# Gap between lower tank floor and control cap top plate
gap_below_cell = openmc.Cell(name="Void between lower tank and control cap")
gap_below_cell.region = (
    (-ltank_inner_cyl & -ltank_floor_bot_plane & +ccap_top_plate_top_plane) &
    ~(-ccap_outer_cyl & -ccap_top_plate_top_plane & +ccap_floor_bot_plane)
)
all_cells.append(gap_below_cell)

# ---- Void outside all structures (within outer sphere) ----
# Everything inside the outer sphere that is not part of any defined cell
void_cell = openmc.Cell(name="Void outside all structures")
void_cell.region = (
    -outer_sphere &
    ~vessel_side_region &
    ~vessel_bottom_region
)
# This is complex to define exactly, so we use a simpler bounding approach:
# The outer void is everything inside the sphere but outside the outermost
# tank radii and above/below the tank axial extents.
# For simplicity, we'll use a catch-all approach.

# Let me rebuild the geometry more carefully using a simpler approach.
# I'll clear all_cells and start over with a cleaner structure.

all_cells.clear()

# =============================================================================
# SIMPLIFIED GEOMETRY BUILD (cleaner cell definitions)
# =============================================================================

# -- Place 28 fuel elements --
fuel_bounding_cyls = {}
for pos_num in range(1, 29):
    x, y = fuel_positions[pos_num]
    cyl = openmc.ZCylinder(x0=x, y0=y, r=clad_outer_radius)
    fuel_bounding_cyls[pos_num] = cyl
    cell = openmc.Cell(name=f"Fuel element position {pos_num}")
    cell.fill = fuel_element_universe
    cell.region = -cyl & -elem_top_plane & +elem_bot_plane
    cell.translation = (x, y, 0.0)
    all_cells.append(cell)

# -- Water inside vessel between fuel elements (core is flooded) --
# Region: inside vessel_inner_cyl, between grid plates, excluding all fuel elements
interior_region = -vessel_inner_cyl & -elem_top_plane & +elem_bot_plane
for pos_num in range(1, 29):
    interior_region = interior_region & +fuel_bounding_cyls[pos_num]
vessel_water = openmc.Cell(name="Water between fuel elements (flooded core)")
vessel_water.fill = water
vessel_water.region = interior_region
all_cells.append(vessel_water)

# -- Upper grid plate (aluminum) --
c_upper_grid = openmc.Cell(name="Upper aluminum grid plate (0.375 in. thick)")
c_upper_grid.fill = aluminum
c_upper_grid.region = -vessel_inner_cyl & +elem_top_plane & -ugrid_top_plane
all_cells.append(c_upper_grid)

# -- Lower grid plate (aluminum) --
c_lower_grid = openmc.Cell(name="Lower aluminum grid plate (0.312 in. thick)")
c_lower_grid.fill = aluminum
c_lower_grid.region = -vessel_inner_cyl & -elem_bot_plane & +lgrid_bot_plane
all_cells.append(c_lower_grid)

# -- Water above upper grid plate (inside vessel) --
c_water_above = openmc.Cell(name="Water above upper grid plate (inside vessel)")
c_water_above.fill = water
c_water_above.region = -vessel_inner_cyl & +ugrid_top_plane & -cover_iface_plane
all_cells.append(c_water_above)

# -- Water below lower grid plate (inside vessel) --
c_water_below = openmc.Cell(name="Water below lower grid plate (inside vessel)")
c_water_below.fill = water
c_water_below.region = -vessel_inner_cyl & -lgrid_bot_plane & +vessel_bot_top_plane
all_cells.append(c_water_below)

# -- SS316 vessel cylindrical side wall --
c_vessel_side = openmc.Cell(name="SS316 vessel side wall")
c_vessel_side.fill = ss316
c_vessel_side.region = (+vessel_inner_cyl & -vessel_outer_cyl &
                        +vessel_bot_bot_plane & -flange_bot_plane)
all_cells.append(c_vessel_side)

# -- SS316 vessel bottom plate --
c_vessel_bottom = openmc.Cell(name="SS316 vessel bottom plate")
c_vessel_bottom.fill = ss316
c_vessel_bottom.region = (-vessel_outer_cyl & +vessel_bot_bot_plane & -vessel_bot_top_plane)
all_cells.append(c_vessel_bottom)

# -- Aluminum cover plate + flange region --
# This is simplified as a single aluminum/SS316 region at the top of the vessel.
# The cover plate mates with the flange. Per the MCNP model, the cover plate
# is aluminum and the flange is part of the vessel (SS316).
c_cover_plate = openmc.Cell(name="Aluminum cover plate")
c_cover_plate.fill = aluminum
c_cover_plate.region = ((-vessel_inner_cyl & +cover_iface_plane & -cover_top_plane) |
                        (-flange_cyl & +vessel_inner_cyl & +cover_iface_plane & -cover_top_plane))
all_cells.append(c_cover_plate)

# -- SS316 vessel flange (below cover plate, around vessel) --
c_vessel_flange = openmc.Cell(name="SS316 vessel flange")
c_vessel_flange.fill = ss316
c_vessel_flange.region = (-flange_cyl & +vessel_inner_cyl &
                          +flange_bot_plane & -cover_iface_plane)
all_cells.append(c_vessel_flange)

# -- Upper tank water --
# Water in the annular space between vessel OD and upper tank ID,
# from the upper tank floor to the water level.
# Also water above the cover plate up to the water level.
c_utank_water = openmc.Cell(name="Water in upper tank (radial + axial reflection)")
c_utank_water.fill = water
c_utank_water.region = (
    (-utank_inner_cyl & +utank_floor_top_plane & -utank_water_plane) &
    ~(-vessel_outer_cyl & +vessel_bot_bot_plane & -flange_bot_plane) &
    ~(-flange_cyl & +flange_bot_plane & -cover_top_plane) &
    ~(-vessel_inner_cyl & -cover_iface_plane & +vessel_bot_top_plane)
)
all_cells.append(c_utank_water)

# -- Void above water in upper tank --
c_utank_void = openmc.Cell(name="Void above water in upper tank")
c_utank_void.region = -utank_inner_cyl & +utank_water_plane & -utank_top_plane
all_cells.append(c_utank_void)

# -- SS316 upper tank walls (side wall + floor plate) --
c_utank_walls = openmc.Cell(name="SS316 upper tank walls")
c_utank_walls.fill = ss316
c_utank_walls.region = (
    # Side wall
    (+utank_inner_cyl & -utank_outer_cyl & +utank_floor_top_plane & -utank_top_plane) |
    # Floor plate (annular between vessel OD and upper tank ID)
    (-utank_inner_cyl & +vessel_outer_cyl & +utank_floor_bot_plane & -utank_floor_top_plane)
)
all_cells.append(c_utank_walls)

# -- Lower tank water --
# Water in annular space around vessel in the lower region
c_ltank_water = openmc.Cell(name="Water in lower tank (radial reflection)")
c_ltank_water.fill = water
c_ltank_water.region = (
    (-ltank_inner_cyl & +vessel_outer_cyl &
     +ltank_floor_top_plane & -utank_floor_bot_plane) &
    # Exclude the upper tank floor plate
    ~(-utank_inner_cyl & +vessel_outer_cyl &
      +utank_floor_bot_plane & -utank_floor_top_plane) &
    # Exclude vessel side wall
    ~(+vessel_inner_cyl & -vessel_outer_cyl &
      +vessel_bot_bot_plane & -flange_bot_plane) &
    # Exclude vessel bottom
    ~(-vessel_outer_cyl & +vessel_bot_bot_plane & -vessel_bot_top_plane)
)
all_cells.append(c_ltank_water)

# -- SS316 lower tank walls --
c_ltank_walls = openmc.Cell(name="SS316 lower tank walls")
c_ltank_walls.fill = ss316
c_ltank_walls.region = (
    # Side wall
    (+ltank_inner_cyl & -ltank_outer_cyl &
     +ltank_floor_top_plane & -utank_floor_bot_plane) |
    # Floor plate
    (-ltank_inner_cyl & +ltank_floor_bot_plane & -ltank_floor_top_plane)
)
all_cells.append(c_ltank_walls)

# -- Control cap tank water --
c_ccap_water = openmc.Cell(name="Water in control cap tank (axial reflection)")
c_ccap_water.fill = water
c_ccap_water.region = (-ccap_inner_cyl &
                       +ccap_floor_top_plane & -ccap_top_plate_bot_plane)
all_cells.append(c_ccap_water)

# -- Control cap tank top plate (aluminum) --
c_ccap_top = openmc.Cell(name="Aluminum top plate of control cap tank")
c_ccap_top.fill = aluminum
c_ccap_top.region = (-ccap_inner_cyl &
                     +ccap_top_plate_bot_plane & -ccap_top_plate_top_plane)
all_cells.append(c_ccap_top)

# -- SS316 control cap tank walls --
c_ccap_walls = openmc.Cell(name="SS316 control cap tank walls")
c_ccap_walls.fill = ss316
c_ccap_walls.region = (
    # Side wall
    (+ccap_inner_cyl & -ccap_outer_cyl &
     +ccap_floor_top_plane & -ccap_top_plate_top_plane) |
    # Floor plate
    (-ccap_inner_cyl & +ccap_floor_bot_plane & -ccap_floor_top_plane)
)
all_cells.append(c_ccap_walls)

# -- Void regions to fill gaps --
# Between lower tank floor and control cap top plate
c_gap_ltank_ccap = openmc.Cell(name="Void between lower tank and control cap")
c_gap_ltank_ccap.region = (
    -outer_sphere &
    +ltank_floor_bot_plane & -ltank_floor_top_plane &
    +ltank_inner_cyl
)
# Actually need to handle the void around all structures more carefully.
# Let's define the overall void as everything inside the sphere but outside
# all defined structures.

# Collect all defined regions for the void complement
# For simplicity, define a single void cell as the complement of all material cells
c_void = openmc.Cell(name="Void (outside all structures, inside boundary)")
void_region = -outer_sphere
for c in all_cells:
    if c.region is not None:
        void_region = void_region & ~c.region
c_void.region = void_region
all_cells.append(c_void)

# Add all cells to root universe
root_universe.add_cells(all_cells)

# Create geometry
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# SETTINGS
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "eigenvalue"

# Criticality source settings
settings.particles = args.particles     # Neutrons per batch
settings.batches = args.batches         # Total batches
settings.inactive = args.inactive       # Inactive batches (for source convergence)

# Initial source: a point source at the center of the fuel array
# Position 10 is at (0, 0, 0), which is roughly the center of the array
settings.source = [openmc.IndependentSource(
    space=openmc.stats.Point(xyz=(0.0, 0.0, 0.0))
)]

# Temperature (room temperature, ~293.6 K)
settings.temperature = {"default": 293.6}

settings.export_to_xml()

# =============================================================================
# Summary
# =============================================================================
print("=" * 70)
print("SNAP-10A/2 SCA-4B Benchmark - Case 8490a")
print("=" * 70)
print(f"  Fuel elements:     28 (positions 1-28 of 37)")
print(f"  Be inserts:        None")
print(f"  Lucite rods:       None")
print(f"  Core flooding:     Water-flooded")
print(f"  Water reflection:  Full (9 in. in upper tank)")
print(f"  Particles/batch:   {args.particles}")
print(f"  Total batches:     {args.batches}")
print(f"  Inactive batches:  {args.inactive}")
print(f"  Reference k_exp:   0.99984 (-0.02$ subcritical)")
print(f"  Reference k_calc:  1.00491 +/- 0.00069 (MCNP5)")
print("=" * 70)
print()
print("XML files written: materials.xml, geometry.xml, settings.xml")
print("Run with: openmc")
