#!/usr/bin/env python3
"""
ST-FNSF Benchmark: Spherical Tokamak Fusion Nuclear Science Facility
=====================================================================

Device Overview
---------------
The Spherical Tokamak Fusion Nuclear Science Facility (ST-FNSF) is a compact
spherical tokamak concept proposed as a stepping stone between present-day
fusion experiments and a full demonstration power plant (DEMO). The design
was developed at the Princeton Plasma Physics Laboratory (PPPL) and is
described in detail by Menard et al. (2016).

Spherical tokamaks (STs) operate at very low aspect ratio (A = R0/a ~ 1.7),
which produces several advantages over conventional aspect-ratio tokamaks:
  - Higher beta (plasma pressure / magnetic pressure) for a given field
  - More compact devices for a given plasma volume
  - Natural elongation of the plasma cross-section
  - Strong natural divertor geometry

The ST-FNSF concept is designed to:
  - Produce 60-100 MW of fusion power in the FNSF operating phase
  - Achieve neutron wall loading of ~1.0 MW/m^2
  - Test tritium breeding blanket modules in a fusion neutron environment
  - Qualify materials and components for DEMO

The extremely compact geometry presents unique neutronics challenges:
  - The center column (R = 0 to ~70 cm) is very cramped, with minimal
    space for shielding between the plasma and the TF coil inner legs
  - No room for a breeding blanket on the inboard side
  - All tritium breeding must occur on the outboard side
  - Center column shielding is a critical design constraint

Device Parameters
-----------------
  Major radius:      R0 = 170 cm (1.7 m)
  Minor radius:      a  = 100 cm (1.0 m)
  Aspect ratio:      A  = 1.7
  Elongation:        kappa = 2.8
  Fusion power:      60-100 MW (FNSF phase)
  Neutron wall load: ~1.0 MW/m^2
  TF coils:          12 (30-degree sector)

Geometry (Simplified CSG Toroidal Sector)
-----------------------------------------
This model uses a 30-degree toroidal sector (360/12 TF coils) with reflective
boundary conditions on the sector planes. The geometry is asymmetric between
the inboard (small-R) and outboard (large-R) sides, which is the defining
characteristic of spherical tokamak neutronics.

  Inboard side (R < 170 cm):
    - Cu center post (TF coil inner legs):  R = 0 to 50 cm
    - WC + water neutron shield:            R = 50 to 70 cm
    - Plasma:                               R = 70 to 170 cm
    No breeding blanket on inboard -- insufficient space.

  Outboard side (R > 170 cm):
    - Plasma:                R = 170 to 270 cm (minor radius 100 cm)
    - F82H first wall:       3 cm thick
    - DCLL blanket:          60 cm thick (LiPb + He + RAFM structure)
    - WC + water shield:     30 cm thick
    - SS316 vacuum vessel:   15 cm thick
    - TF coil return legs:   R > 381 cm

Materials
---------
  Copper center post:    Pure Cu, 8.96 g/cm3
  Inboard shield:        80% WC + 20% water by volume, ~12.7 g/cm3
  F82H RAFM steel:       Fe-Cr-W-V, 7.89 g/cm3
  DCLL blanket:          25% RAFM + 65% LiPb + 10% He (void), ~7.1 g/cm3
  SS316 vacuum vessel:   Standard austenitic steel, 7.93 g/cm3
  LiPb breeder:          Pb-17Li with 90% Li-6 enrichment, 9.4 g/cm3

Reference
---------
  J. E. Menard et al., "Fusion Nuclear Science Facilities and Pilot Plants
  Based on the Spherical Tokamak," Nuclear Fusion, vol. 56, no. 10, p. 106023,
  2016. DOI: 10.1088/0029-5515/56/10/106023
"""

import argparse
import copy
import math

import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="ST-FNSF Benchmark -- OpenMC fixed-source model"
)
parser.add_argument(
    "--particles",
    type=int,
    default=50_000,
    help="Number of source particles per batch (default: 50,000)",
)
parser.add_argument(
    "--batches",
    type=int,
    default=100,
    help="Number of batches (default: 100)",
)
parser.add_argument(
    "--weight-windows",
    action="store_true",
    help="Generate weight windows using FW-CADIS before running",
)
args = parser.parse_args()


# =============================================================================
# Materials
# =============================================================================

# --- Copper center post (TF coil inner legs) ---
# The center column of a spherical tokamak is a single copper post that
# carries the current for all TF coils. In ST-FNSF, this extends from the
# axis (R=0) out to R~50 cm. Pure copper is used for its high electrical
# conductivity.
copper = openmc.Material(name="Copper center post")
copper.add_nuclide("Cu63", 0.6917, "ao")
copper.add_nuclide("Cu65", 0.3083, "ao")
copper.set_density("g/cm3", 8.96)

# --- WC + water inboard shield (80/20 vol%) ---
# The inboard neutron shield sits between the center post and the plasma
# (R = 50 to 70 cm). It must attenuate neutrons to protect the copper
# center post and TF coil superconductor from radiation damage. Tungsten
# carbide provides excellent fast neutron attenuation; water provides
# moderation for intermediate-energy neutrons.
#
# Volume fractions: 80% WC (density 15.63 g/cm3) + 20% H2O (1.0 g/cm3)
# Effective density: 0.80 * 15.63 + 0.20 * 1.0 = 12.7 g/cm3
inboard_shield = openmc.Material(name="WC + water inboard shield")
# WC component (80 vol%): W and C in 1:1 atomic ratio
# Molecular weight WC = 195.85 g/mol, density 15.63 g/cm3
# H2O component (20 vol%): standard water
# We homogenize by weight fractions:
# Mass WC = 0.80 * 15.63 = 12.504 g/cm3
# Mass H2O = 0.20 * 1.0 = 0.200 g/cm3
# Total = 12.704 g/cm3
# Weight fraction WC = 12.504 / 12.704 = 0.9843
# Weight fraction H2O = 0.200 / 12.704 = 0.0157
inboard_shield.add_element("W", 0.9343, "wo")   # tungsten in WC
inboard_shield.add_element("C", 0.0500, "wo")   # carbon in WC
inboard_shield.add_nuclide("H1", 0.0022, "wo")  # hydrogen in water
inboard_shield.add_nuclide("O16", 0.0135, "wo") # oxygen in water
inboard_shield.set_density("g/cm3", 12.7)

# --- F82H RAFM steel (Reduced Activation Ferritic/Martensitic) ---
# F82H is the reference structural steel for fusion blankets. It is a
# reduced-activation alloy where Mo and Nb (which produce long-lived
# activation products) are replaced by W and V/Ta.
# Composition (wt%): Fe 89%, Cr 7.5%, W 2%, V 0.2%, Ta 0.04%,
#                     Mn 0.1%, Si 0.1%, balance Fe
f82h = openmc.Material(name="F82H RAFM steel")
f82h.add_element("Fe", 0.8906, "wo")
f82h.add_element("Cr", 0.0750, "wo")
f82h.add_element("W",  0.0200, "wo")
f82h.add_element("V",  0.0020, "wo")
f82h.add_element("Ta", 0.0004, "wo")
f82h.add_element("Mn", 0.0010, "wo")
f82h.add_element("Si", 0.0010, "wo")
f82h.add_element("C",  0.0100, "wo")
f82h.set_density("g/cm3", 7.89)

# --- LiPb (Pb-17Li) breeder ---
# The eutectic lead-lithium alloy Pb-17Li is the primary tritium breeding
# material in the DCLL (Dual Coolant Lead-Lithium) blanket concept. The
# "17" refers to 17 atom% lithium. For enhanced tritium breeding, the
# lithium is enriched to 90% Li-6 (natural is ~7.5% Li-6).
#
# Composition: 17 at% Li (90% Li-6, 10% Li-7) + 83 at% Pb
# Density: 9.4 g/cm3 at operating temperature (~500 C)
lipb = openmc.Material(name="Pb-17Li breeder")
# Li: 17 at%, enriched to 90% Li-6
lipb.add_nuclide("Li6", 0.17 * 0.90, "ao")   # 15.3 at%
lipb.add_nuclide("Li7", 0.17 * 0.10, "ao")   # 1.7 at%
# Pb: 83 at% (natural isotopic composition)
lipb.add_element("Pb", 0.83, "ao")
lipb.set_density("g/cm3", 9.4)

# --- Homogenized DCLL blanket ---
# The Dual Coolant Lead-Lithium (DCLL) blanket uses LiPb as both breeder
# and coolant, with helium as a secondary coolant for the steel structure.
# Approximate volume fractions:
#   25% F82H RAFM structure
#   65% LiPb (Pb-17Li breeder/coolant)
#   10% He (void -- negligible mass)
#
# Effective density: 0.25 * 7.89 + 0.65 * 9.4 + 0.10 * 0.0 = 8.08 g/cm3
# (The user spec says ~7.1 g/cm3 which accounts for lower packing; we use
#  7.1 as specified.)
dcll_blanket = openmc.Material(name="DCLL blanket (homogenized)")
# F82H component (25 vol%, density 7.89 -> mass contribution 1.9725 g/cm3)
# LiPb component (65 vol%, density 9.4 -> mass contribution 6.11 g/cm3)
# Total effective mass = 1.9725 + 6.11 = 8.0825 g/cm3
# But we use the specified 7.1 g/cm3 to account for structure/void packing
# Weight fractions: F82H = 1.9725/8.0825 = 0.2441, LiPb = 6.11/8.0825 = 0.7559
# Expand F82H: Fe 0.2441*0.89 = 0.2172, Cr 0.2441*0.075 = 0.0183, etc.
# Expand LiPb: need atomic fractions -> convert to weight fractions
# For simplicity, build from atomic fractions of constituents:
dcll_blanket.add_element("Fe", 0.2172, "wo")
dcll_blanket.add_element("Cr", 0.0183, "wo")
dcll_blanket.add_element("W",  0.0049, "wo")
dcll_blanket.add_element("V",  0.0005, "wo")
dcll_blanket.add_element("Ta", 0.0001, "wo")
dcll_blanket.add_element("Mn", 0.0002, "wo")
dcll_blanket.add_element("Si", 0.0002, "wo")
dcll_blanket.add_element("C",  0.0024, "wo")
# LiPb portion (75.59 wt%)
dcll_blanket.add_nuclide("Li6", 0.7559 * 0.17 * 0.90 * 6.015 / 175.82, "wo")
dcll_blanket.add_nuclide("Li7", 0.7559 * 0.17 * 0.10 * 7.016 / 175.82, "wo")
dcll_blanket.add_element("Pb",  0.7559 * 0.83 * 207.2 / 175.82, "wo")
dcll_blanket.set_density("g/cm3", 7.1)

# --- WC + water outboard shield (same composition as inboard) ---
outboard_shield = openmc.Material(name="WC + water outboard shield")
outboard_shield.add_element("W", 0.9343, "wo")
outboard_shield.add_element("C", 0.0500, "wo")
outboard_shield.add_nuclide("H1", 0.0022, "wo")
outboard_shield.add_nuclide("O16", 0.0135, "wo")
outboard_shield.set_density("g/cm3", 12.7)

# --- SS316 vacuum vessel ---
# Standard austenitic stainless steel used for the vacuum vessel.
# Composition (wt%): Fe 65.395%, Cr 17%, Ni 12%, Mo 2.5%,
#                     Mn 2%, Si 1%, C 0.08%, N 0.025%
ss316 = openmc.Material(name="SS316 vacuum vessel")
ss316.add_element("Fe", 0.65395, "wo")
ss316.add_element("Cr", 0.17000, "wo")
ss316.add_element("Ni", 0.12000, "wo")
ss316.add_element("Mo", 0.02500, "wo")
ss316.add_element("Mn", 0.02000, "wo")
ss316.add_element("Si", 0.01000, "wo")
ss316.add_element("C",  0.00080, "wo")
ss316.add_element("N",  0.00025, "wo")
ss316.set_density("g/cm3", 7.93)

# --- Collect and export materials ---
materials = openmc.Materials([
    copper, inboard_shield, f82h, lipb, dcll_blanket,
    outboard_shield, ss316,
])
materials.cross_sections = "/data/endfb-viii.0-hdf5/cross_sections.xml"
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# The ST-FNSF geometry uses a 30-degree toroidal sector with reflective
# boundaries on the sector faces. The key feature of a spherical tokamak is
# the extreme asymmetry between the inboard and outboard sides.
#
# Radial layout (at midplane, z=0):
#   R = 0 to 50 cm:    Copper center post
#   R = 50 to 70 cm:   WC+water inboard shield
#   R = 70 to 270 cm:  Plasma (void) -- from R0-a to R0+a
#   R = 270 to 273 cm: F82H first wall (outboard only)
#   R = 273 to 333 cm: DCLL blanket (outboard only)
#   R = 333 to 363 cm: WC+water outboard shield
#   R = 363 to 378 cm: SS316 vacuum vessel
#   R > 378 cm:        Void to boundary

# --- Center column surfaces ---
# The center post is a vertical cylinder at the tokamak axis.
center_post_outer = openmc.ZCylinder(r=50.0, name="Center post outer surface")
inboard_shield_outer = openmc.ZCylinder(r=70.0, name="Inboard shield outer surface")

# --- Plasma boundary ---
# The plasma is bounded by a torus surface. For a circular cross-section
# plasma: ZTorus(a=R0, b=a, c=a) where R0 is major radius and a is minor radius.
# a=170 is the distance from the z-axis to the torus center,
# b=100 and c=100 define circular cross-section with minor radius 100 cm.
plasma_inner = openmc.ZTorus(a=170.0, b=100.0, c=100.0,
                             name="Plasma boundary")

# --- Concentric torus shells for outboard components ---
# These are concentric tori centered at R0=170 cm with increasing minor radii.
# Each shell represents a different component layer on the outboard side.
fw_outer = openmc.ZTorus(a=170.0, b=103.0, c=103.0,
                         name="First wall outer surface")
blanket_outer = openmc.ZTorus(a=170.0, b=163.0, c=163.0,
                              name="Blanket outer surface")
shield_outer = openmc.ZTorus(a=170.0, b=193.0, c=193.0,
                             name="Outboard shield outer surface")
vv_outer = openmc.ZTorus(a=170.0, b=211.0, c=211.0,
                         name="Vacuum vessel outer surface")

# --- Axial boundaries ---
z_top = openmc.ZPlane(z0=500.0, boundary_type="vacuum",
                      name="Top boundary")
z_bottom = openmc.ZPlane(z0=-500.0, boundary_type="vacuum",
                         name="Bottom boundary")

# --- Outer radial boundary ---
r_outer = openmc.ZCylinder(r=400.0, boundary_type="vacuum",
                           name="Outer radial boundary")

# --- Sector boundaries (30-degree toroidal sector) ---
# The sector spans from phi=0 to phi=30 degrees.
# First boundary: the y=0 plane (phi=0)
# Second boundary: a plane at phi=30 degrees from the x-axis
# Both are reflective to simulate the full 360 degrees with 12-fold symmetry.
sector_angle = 30.0  # degrees
phi = math.radians(sector_angle)

# y=0 plane: normal vector (0, 1, 0), passing through origin
plane_0 = openmc.YPlane(y0=0.0, boundary_type="reflective",
                        name="Sector boundary at phi=0")

# Plane at phi=30 degrees: normal vector (sin(phi), -cos(phi), 0)
# This plane contains the z-axis and makes angle phi with the x-axis.
plane_30 = openmc.Plane(
    a=math.sin(phi), b=-math.cos(phi), c=0.0, d=0.0,
    boundary_type="reflective",
    name="Sector boundary at phi=30 deg",
)

# --- Define the sector region ---
# The sector is the region between the two planes, above z_bottom,
# below z_top, and inside r_outer.
sector_region = +plane_0 & +plane_30 & +z_bottom & -z_top & -r_outer

# =============================================================================
# Cell definitions
# =============================================================================
# The geometry has a fundamentally asymmetric structure: on the inboard side
# (R < R0), only the center post and shield exist. On the outboard side
# (R > R0), the full blanket/shield/VV stack surrounds the plasma.
# The torus surfaces naturally handle this asymmetry since they close
# around the plasma on all sides, but the inboard is bounded by the
# center column cylinders.

# --- Cell 1: Copper center post ---
# Everything inside the center post cylinder, within the sector
center_post_cell = openmc.Cell(
    name="Copper center post",
    fill=copper,
    region=-center_post_outer & sector_region,
)

# --- Cell 2: Inboard neutron shield ---
# Annular region between center post and inboard shield boundary
inboard_shield_cell = openmc.Cell(
    name="Inboard neutron shield (WC+water)",
    fill=inboard_shield,
    region=+center_post_outer & -inboard_shield_outer & sector_region,
)

# --- Cell 3: Plasma (void) ---
# Inside the plasma torus, but outside the inboard shield cylinder.
# On the inboard side, the plasma extends from R=70 (inboard_shield_outer)
# to the torus surface. On the outboard side, it extends to R=270.
plasma_cell = openmc.Cell(
    name="Plasma (void)",
    region=-plasma_inner & +inboard_shield_outer & sector_region,
)

# --- Cell 4: First wall (F82H, outboard only) ---
# Between plasma surface and first wall outer surface.
# The cylinder cut ensures this only exists on the outboard side
# (for the inboard portion, the shield cylinder already bounds things).
fw_cell = openmc.Cell(
    name="First wall (F82H)",
    fill=f82h,
    region=+plasma_inner & -fw_outer & sector_region,
)

# --- Cell 5: DCLL blanket (outboard) ---
blanket_cell = openmc.Cell(
    name="DCLL blanket (outboard)",
    fill=dcll_blanket,
    region=+fw_outer & -blanket_outer & sector_region,
)

# --- Cell 6: Outboard shield (WC + water) ---
outboard_shield_cell = openmc.Cell(
    name="Outboard shield (WC+water)",
    fill=outboard_shield,
    region=+blanket_outer & -shield_outer & sector_region,
)

# --- Cell 7: Vacuum vessel (SS316) ---
vv_cell = openmc.Cell(
    name="Vacuum vessel (SS316)",
    fill=ss316,
    region=+shield_outer & -vv_outer & sector_region,
)

# --- Cell 8: Void (everything else in the sector) ---
# Outside the VV but inside the sector boundaries, and not in the
# center post or inboard shield regions.
void_cell = openmc.Cell(
    name="Void (outside VV)",
    region=(
        +vv_outer
        & +inboard_shield_outer
        & sector_region
    ),
)

# --- Cell 9: Inboard void ---
# The region between the inboard shield and the plasma on the inboard
# side is already covered by the plasma cell. But we need to handle
# the region outside the plasma torus on the inboard side (above/below
# the plasma) that is between the shield cylinder and the first wall torus.
# This is void (no material between inboard shield and plasma vertically).
inboard_void_cell = openmc.Cell(
    name="Inboard void (above/below plasma)",
    region=(
        +plasma_inner
        & +vv_outer
        & -inboard_shield_outer
        & sector_region
    ),
)

# Actually, the geometry needs careful handling. Let me restructure:
# The torus surfaces close around, so +plasma_inner means OUTSIDE the torus.
# Between the inboard shield cylinder and the VV torus outer surface,
# excluding the interior torus regions, we need a void cell.

# Rebuild cells more carefully:
# Region inside inboard shield cylinder is: center post + inboard shield
# Region inside plasma torus (but outside inboard shield cyl): plasma void
# Region between plasma torus and fw torus: first wall
# Region between fw torus and blanket torus: blanket
# Region between blanket torus and shield torus: outboard shield
# Region between shield torus and vv torus: vacuum vessel
# Everything else in the sector: void

# The "everything else" is:
#  (outside vv_outer torus OR inside inboard_shield_outer cylinder)
#  AND in sector
#  minus center post, minus inboard shield, minus plasma
# But center post and inboard shield are inside inboard_shield_outer,
# and plasma is inside plasma_inner, so:

# Void = in sector AND outside vv_outer AND outside inboard_shield_outer
#   PLUS: in sector AND inside inboard_shield_outer BUT that's already
#         covered by center_post and inboard_shield cells.
# The tricky region: between inboard shield cylinder outer surface and
# the plasma torus, but OUTSIDE the plasma torus (above/below the torus).
# This is: +inboard_shield_outer & +plasma_inner & -vv_outer & sector
# (outside inboard shield cyl, outside plasma torus, inside vv torus)

# Let's redefine the void regions properly:
# Region A: outside VV torus, outside inboard shield cyl -> void
# Region B: outside plasma torus, inside VV torus structure layers are
#           already covered by fw, blanket, shield, vv cells
# Region C: outside plasma torus, outside inboard shield cyl,
#           inside the VV torus but NOT in any shell -> this is the
#           inboard gap between shield cyl and torus shells

# Actually, the torus shells (fw, blanket, shield, vv) wrap around
# the full poloidal cross-section including inboard. But we only want
# blanket material on the outboard side. This is the subtlety of ST design.
#
# For this simplified CSG model, the concentric torus shells DO wrap
# around to the inboard side. On the inboard side, the torus shells
# overlap with the center column cylinders. The cylinder cells take
# priority because they are defined first with the cylinder boundary.
# The torus shell cells only occupy the region OUTSIDE the inboard
# shield cylinder, which is physically correct -- the blanket/shield
# structures only exist where there is space outside the center column.
#
# However, at the very top and bottom of the torus, the shells may
# extend inboard of the shield cylinder (at high |z|). In a real
# ST-FNSF, the blanket does not extend to the very top/bottom on the
# inboard side. For this simplified model, this is acceptable.

# Redefine cells to properly partition space:

# Clear previous cell definitions and rebuild
center_post_cell = openmc.Cell(
    name="Copper center post",
    fill=copper,
    region=-center_post_outer & sector_region,
)

inboard_shield_cell = openmc.Cell(
    name="Inboard neutron shield (WC+water)",
    fill=inboard_shield,
    region=+center_post_outer & -inboard_shield_outer & sector_region,
)

# Plasma: inside plasma torus AND outside inboard shield cylinder
plasma_cell = openmc.Cell(
    name="Plasma (void)",
    region=-plasma_inner & +inboard_shield_outer & sector_region,
)

# First wall: between plasma torus and fw torus, outside inboard shield cyl
fw_cell = openmc.Cell(
    name="First wall (F82H)",
    fill=f82h,
    region=+plasma_inner & -fw_outer & +inboard_shield_outer & sector_region,
)

# Blanket: between fw torus and blanket torus, outside inboard shield cyl
blanket_cell = openmc.Cell(
    name="DCLL blanket",
    fill=dcll_blanket,
    region=+fw_outer & -blanket_outer & +inboard_shield_outer & sector_region,
)

# Outboard shield: between blanket torus and shield torus, outside inboard shield cyl
outboard_shield_cell = openmc.Cell(
    name="Outboard shield (WC+water)",
    fill=outboard_shield,
    region=+blanket_outer & -shield_outer & +inboard_shield_outer & sector_region,
)

# Vacuum vessel: between shield torus and vv torus, outside inboard shield cyl
vv_cell = openmc.Cell(
    name="Vacuum vessel (SS316)",
    fill=ss316,
    region=+shield_outer & -vv_outer & +inboard_shield_outer & sector_region,
)

# Void: outside vv torus AND outside inboard shield cyl, in sector
void_cell = openmc.Cell(
    name="Void (outside VV, outboard)",
    region=+vv_outer & +inboard_shield_outer & sector_region,
)

# --- Build geometry ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[
        center_post_cell,
        inboard_shield_cell,
        plasma_cell,
        fw_cell,
        blanket_cell,
        outboard_shield_cell,
        vv_cell,
        void_cell,
    ],
)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The D-T fusion source is distributed throughout the plasma volume.
# For a spherical tokamak, the plasma is centered at R0=170 cm with
# minor radius a=100 cm, so it extends radially from R~70 to R~270 cm.
#
# We use a cylindrical independent source distribution:
#   R:   uniform from 100 to 240 cm (peaked region, not full extent)
#   phi: uniform from 0 to 30 degrees (one sector)
#   z:   uniform from -80 to +80 cm (within elongated plasma)
#
# Energy: 14.1 MeV (D-T fusion neutrons)
# Direction: isotropic

source = openmc.IndependentSource()
r_dist = openmc.stats.Uniform(a=100.0, b=240.0)
phi_dist = openmc.stats.Uniform(a=0.0, b=math.radians(30.0))
z_dist = openmc.stats.Uniform(a=-80.0, b=80.0)
source.space = openmc.stats.CylindricalIndependent(
    r=r_dist, phi=phi_dist, z=z_dist
)
source.energy = openmc.stats.Discrete([14.1e6], [1.0])
source.angle = openmc.stats.Isotropic()
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = source
settings.particles = args.particles
settings.batches = args.batches
settings.photon_transport = False
settings.temperature = {"default": 600.0}
settings.output = {"tallies": True}

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# Key tallies for ST-FNSF neutronics:
#
# 1. Tritium Breeding Ratio (TBR):
#    (n,Xt) reaction rate in the blanket. This is the single most important
#    figure of merit -- ST-FNSF needs TBR >= 1.0 for tritium self-sufficiency
#    but achieves only marginal breeding due to limited blanket coverage
#    (no inboard blanket).
#
# 2. Nuclear heating in all components:
#    Determines thermal-hydraulic requirements and helps identify hot spots.
#    Center column heating is especially critical because cooling is difficult.
#
# 3. Neutron flux at center column:
#    Fast neutron damage to the copper center post is a lifetime-limiting
#    factor. The compact inboard geometry means the center post sees
#    significant neutron flux despite the shield.
#
# 4. Fast neutron flux (E > 0.1 MeV) at the TF coil region:
#    Superconductor radiation damage limit is typically ~1e18 n/cm2.
#
# 5. Neutron current through inboard/outboard surfaces:
#    Quantifies the asymmetric neutron wall loading.

tallies = openmc.Tallies()

# --- Cell filters ---
blanket_filter = openmc.CellFilter([blanket_cell])
center_post_filter = openmc.CellFilter([center_post_cell])
inboard_shield_filter = openmc.CellFilter([inboard_shield_cell])
fw_filter = openmc.CellFilter([fw_cell])
outboard_shield_filter = openmc.CellFilter([outboard_shield_cell])
vv_filter = openmc.CellFilter([vv_cell])

# --- Tally 1: Tritium Breeding Ratio (TBR) ---
# Score (n,Xt) which counts all tritium-producing reactions including
# (n,t), (n,nt), (n,2nt), etc. This is the standard way to compute TBR.
tbr_tally = openmc.Tally(name="TBR")
tbr_tally.filters = [blanket_filter]
tbr_tally.scores = ["(n,Xt)"]
tallies.append(tbr_tally)

# --- Tally 2: Nuclear heating in all components ---
# Heating scores the energy deposited per source particle [eV/source].
heating_tally = openmc.Tally(name="nuclear_heating")
all_cells_filter = openmc.CellFilter([
    center_post_cell, inboard_shield_cell, fw_cell,
    blanket_cell, outboard_shield_cell, vv_cell,
])
heating_tally.filters = [all_cells_filter]
heating_tally.scores = ["heating"]
tallies.append(heating_tally)

# --- Tally 3: Neutron flux at center column ---
# The center post flux determines radiation damage rates (dpa) and
# nuclear heating in the copper conductor.
center_flux_tally = openmc.Tally(name="center_column_flux")
center_flux_tally.filters = [center_post_filter]
center_flux_tally.scores = ["flux"]
tallies.append(center_flux_tally)

# --- Tally 4: Neutron flux spectrum at center column ---
# Energy-resolved flux to characterize the neutron spectrum reaching
# the center post after penetrating the inboard shield.
energy_bins = np.logspace(np.log10(1.0e3), np.log10(15.0e6), 101)
energy_filter = openmc.EnergyFilter(energy_bins)

center_spectrum_tally = openmc.Tally(name="center_column_spectrum")
center_spectrum_tally.filters = [center_post_filter, energy_filter]
center_spectrum_tally.scores = ["flux"]
tallies.append(center_spectrum_tally)

# --- Tally 5: Fast flux (E > 0.1 MeV) at VV (proxy for TF coil) ---
# The vacuum vessel is the last material layer before the TF coils.
# Fast flux here indicates what reaches the superconductor.
fast_energy_filter = openmc.EnergyFilter([0.1e6, 15.0e6])
fast_flux_tally = openmc.Tally(name="fast_flux_at_vv")
fast_flux_tally.filters = [vv_filter, fast_energy_filter]
fast_flux_tally.scores = ["flux"]
tallies.append(fast_flux_tally)

# --- Tally 6: Inboard shield heating ---
# Heating in the inboard shield is a critical design parameter because
# the compact geometry concentrates nuclear heating.
inboard_heating_tally = openmc.Tally(name="inboard_shield_heating")
inboard_heating_tally.filters = [inboard_shield_filter]
inboard_heating_tally.scores = ["heating"]
tallies.append(inboard_heating_tally)

# --- Tally 7: Blanket neutron flux spectrum ---
# Characterizes the neutron spectrum in the breeding zone.
blanket_spectrum_tally = openmc.Tally(name="blanket_spectrum")
blanket_spectrum_tally.filters = [blanket_filter, energy_filter]
blanket_spectrum_tally.scores = ["flux"]
tallies.append(blanket_spectrum_tally)

tallies.export_to_xml()


# =============================================================================
# Weight window generation (optional)
# =============================================================================
if args.weight_windows:
    print("\n" + "=" * 70)
    print("Generating weight windows using FW-CADIS...")
    print("=" * 70)

    model = openmc.Model(
        geometry=geometry, materials=materials, settings=settings, tallies=tallies
    )

    model_rr = copy.deepcopy(model)
    model_rr.convert_to_multigroup(
        method="material_wise", groups="CASMO-4", nparticles=3000
    )
    model_rr.convert_to_random_ray()

    bbox = model_rr.geometry.bounding_box
    nx = max(1, int((bbox.upper_right[0] - bbox.lower_left[0]) / 10.0))
    ny = max(1, int((bbox.upper_right[1] - bbox.lower_left[1]) / 10.0))
    nz = max(1, int((bbox.upper_right[2] - bbox.lower_left[2]) / 10.0))

    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (nx, ny, nz)
    ww_mesh.lower_left = list(bbox.lower_left)
    ww_mesh.upper_right = list(bbox.upper_right)

    root = model_rr.geometry.root_universe
    model_rr.settings.random_ray["source_region_meshes"] = [(ww_mesh, [root])]

    wwg = openmc.WeightWindowGenerator(
        method="fw_cadis",
        mesh=ww_mesh,
        max_realizations=model_rr.settings.batches,
    )
    model_rr.settings.weight_window_generators = wwg
    model_rr.export_to_xml()
    openmc.run()

    print("Weight windows generated successfully.")


# =============================================================================
# Summary
# =============================================================================
print("=" * 70)
print("ST-FNSF Benchmark -- OpenMC Model")
print("=" * 70)
print(f"  Major radius R0:   170.0 cm")
print(f"  Minor radius a:    100.0 cm")
print(f"  Aspect ratio A:    1.7")
print(f"  Sector angle:      {sector_angle} degrees (12 TF coils)")
print(f"  Components:")
print(f"    Center post:     R = 0 to 50 cm (Cu)")
print(f"    Inboard shield:  R = 50 to 70 cm (WC+water)")
print(f"    Plasma:          R = 70 to 270 cm (void)")
print(f"    First wall:      +3 cm (F82H)")
print(f"    DCLL blanket:    +60 cm (LiPb + RAFM + He)")
print(f"    Outboard shield: +30 cm (WC+water)")
print(f"    Vacuum vessel:   +18 cm (SS316)")
print(f"  Source:            14.1 MeV D-T, R = 100-240 cm, z = -80 to 80 cm")
print(f"  Particles:         {args.particles:,} per batch x {args.batches} batches")
print(f"  Temperature:       600 K")
print(f"  Cross sections:    ENDF/B-VIII.0")
print(f"  Tallies:           TBR, heating (6 components), flux spectra")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
