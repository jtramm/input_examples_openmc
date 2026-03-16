#!/usr/bin/env python3
"""
CFETR 22.5-Degree Toroidal Sector Fusion Benchmark
=====================================================

Facility
--------
The China Fusion Engineering Test Reactor (CFETR) is a proposed tokamak
designed by the Chinese fusion community to bridge the gap between ITER and a
demonstration fusion power plant (DEMO). CFETR aims to demonstrate steady-state
operation, electricity generation, and tritium self-sufficiency. The conceptual
design features a major radius of 5.6 m, a plasma minor radius of 1.3 m, and
a fusion power ranging from 200 MW (Phase I) to 1000 MW (Phase II).

Model
-----
This benchmark creates a simplified 22.5-degree toroidal sector (1/16 of the
full torus, corresponding to one TF coil period) of the CFETR tokamak with
the Water-Cooled Ceramic Breeder (WCCB) blanket concept. The WCCB blanket
uses Li4SiO4 ceramic pebbles as the tritium breeder, beryllium pebbles as
the neutron multiplier, RAFM steel (CLF-1) as the structural material, and
pressurized water as the coolant.

The geometry is built from concentric ZTorus surfaces (toroidal shells) cut
by two reflective planes that define the 22.5-degree sector. Each layer
represents a major functional component of the tokamak: plasma chamber,
first wall, breeding blanket, back support structure, vacuum vessel, thermal
shield, and TF coil casing.

Tritium Breeding Physics
------------------------
In a D-T fusion reactor, each fusion event produces one 14.1 MeV neutron and
one 3.5 MeV alpha particle. Since tritium has a half-life of only 12.3 years
and does not occur naturally in useful quantities, the blanket must breed
tritium from lithium to close the fuel cycle:

  Li-6 + n -> He-4 + T + 4.78 MeV   (exothermic, large thermal cross section)
  Li-7 + n -> He-4 + T + n' - 2.47 MeV  (endothermic, threshold ~2.5 MeV)

The Li4SiO4 breeder in the WCCB blanket uses lithium enriched to 60% Li-6
to enhance the breeding ratio. Beryllium pebbles serve as the neutron
multiplier via:

  Be-9 + n -> 2 He-4 + 2n - 1.57 MeV  (threshold ~1.85 MeV)

This (n,2n) reaction is essential for achieving a Tritium Breeding Ratio
(TBR) greater than unity, as it compensates for neutron losses to parasitic
absorption and leakage.

Geometry
--------
  Coordinate system: Cartesian (x, y, z) with z as the tokamak symmetry axis.
  The toroidal sector spans from phi=0 (the y=0 plane) to phi=22.5 degrees.

  Radial build (from plasma centre outward), all distances are minor radii
  measured from the plasma magnetic axis at R0=560 cm:

    Component              Inner r (cm)   Outer r (cm)   Thickness (cm)
    -----------------------------------------------------------------
    Plasma (void)              0              130             130
    First Wall (CLF-1)       130              132               2
    Breeding Blanket (hom.)  132              182              50
    Back Support Structure   182              197              15
    VV inner wall (SS316)    197              201               4
    VV fill (SS316+H2O)     201              221              20
    VV outer wall (SS316)   221              225               4
    Thermal Shield (SS304)  225              227               2
    TF Coil (SS316+Cu)      227              252              25

  Sector boundaries:
    phi = 0 degrees:    y = 0 plane (reflective)
    phi = 22.5 degrees: general plane through origin (reflective)

  The bounding box is closed by a large vacuum sphere.

Materials
---------
  CLF-1 / RAFM steel:  China Low Activation Ferritic steel, 7.8 g/cm3
  Li4SiO4 breeder:     90% TD pebbles, 60% Li-6 enrichment, 2.32 g/cm3
  Beryllium pebbles:   63% packing fraction, effective 1.17 g/cm3
  Homogenized blanket: 30% CLF-1 + 30% Li4SiO4 + 30% Be + 10% H2O
  SS316LN:             Austenitic stainless steel, 7.93 g/cm3
  SS304:               Austenitic stainless steel, 7.93 g/cm3
  Water (H2O):         Light water coolant at 300 K, 1.0 g/cm3

Reference
---------
  S. Zhu et al., "Design and R&D progress of the CFETR blanket," Plasma Sci.
  Technol. 18 (2016) 13. DOI: 10.1088/1009-0630/18/7/13

  Y. Wan et al., "Overview of the present progress and activities on the CFETR,"
  Nucl. Fusion 57 (2017) 102009.
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
    description="CFETR 22.5-degree sector -- OpenMC fixed-source model"
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
    "--generate-ww",
    action="store_true",
    help="Generate weight windows using FW-CADIS before running",
)
args = parser.parse_args()


# =============================================================================
# Key device parameters
# =============================================================================
# CFETR tokamak geometry (Zhu et al. 2016)
R0 = 560.0       # Major radius [cm]
a_plasma = 130.0  # Plasma minor radius [cm]
kappa = 1.8       # Elongation (not used in circular ZTorus approximation)
sector_angle = 22.5  # Toroidal sector angle [degrees] (360/16 TF coils)
phi_rad = math.radians(sector_angle)

# Radial build thicknesses [cm]
# Each layer is defined by its inner minor radius and thickness.
# The minor radius is measured from the plasma magnetic axis.
t_fw = 2.0         # First wall (CLF-1 RAFM steel)
t_blanket = 50.0   # Breeding blanket (homogenized WCCB)
t_bss = 15.0       # Back support structure (CLF-1 + water)
t_vv_inner = 4.0   # Vacuum vessel inner wall (SS316)
t_vv_fill = 20.0   # Vacuum vessel fill (SS316 + water)
t_vv_outer = 4.0   # Vacuum vessel outer wall (SS316)
t_shield = 2.0     # Thermal shield (SS304)
t_tf = 25.0        # TF coil casing (SS316 + Cu)

# Cumulative minor radii for each interface
r_plasma = a_plasma                                          # 130 cm
r_fw = r_plasma + t_fw                                       # 132 cm
r_blanket = r_fw + t_blanket                                 # 182 cm
r_bss = r_blanket + t_bss                                    # 197 cm
r_vv_inner = r_bss + t_vv_inner                              # 201 cm
r_vv_fill = r_vv_inner + t_vv_fill                           # 221 cm
r_vv_outer = r_vv_fill + t_vv_outer                          # 225 cm
r_shield = r_vv_outer + t_shield                             # 227 cm
r_tf = r_shield + t_tf                                       # 252 cm


# =============================================================================
# Materials
# =============================================================================
# All materials use ENDF/B-VIII.0 cross sections.

# --- CLF-1 / RAFM steel (China Low Activation Ferritic steel) ---
# CLF-1 is a reduced activation ferritic-martensitic (RAFM) steel developed
# in China as the primary structural material for fusion blankets. It is
# analogous to EUROFER-97 (EU) and F82H (Japan). The key design principle is
# to replace high-activation elements (Mo, Nb, Ni, Co) with low-activation
# alternatives (W, V, Ta) while maintaining mechanical properties.
#
# Composition (weight percent):
#   Fe: 88.5% (balance), Cr: 9.0%, W: 1.5%, V: 0.2%, Ta: 0.1%,
#   Mn: 0.5%, C: 0.1%, Si: 0.1%
# Density: 7.8 g/cm3
clf1_steel = openmc.Material(name="CLF-1 RAFM steel")
clf1_steel.add_element("Fe", 0.885, "wo")   # balance
clf1_steel.add_element("Cr", 0.090, "wo")   # 9% chromium
clf1_steel.add_element("W",  0.015, "wo")   # 1.5% tungsten (low activation)
clf1_steel.add_element("V",  0.002, "wo")   # 0.2% vanadium
clf1_steel.add_element("Ta", 0.001, "wo")   # 0.1% tantalum
clf1_steel.add_element("Mn", 0.005, "wo")   # 0.5% manganese
clf1_steel.add_nuclide("C12", 0.001 * 0.9893, "wo")  # carbon-12
clf1_steel.add_nuclide("C13", 0.001 * 0.0107, "wo")  # carbon-13
clf1_steel.add_element("Si", 0.001, "wo")   # 0.1% silicon
clf1_steel.set_density("g/cm3", 7.8)

# --- Li4SiO4 ceramic breeder pebbles ---
# Lithium orthosilicate (Li4SiO4) is the primary tritium breeder material
# in the WCCB blanket concept. The pebbles are fabricated at ~90% of
# theoretical density (2.58 g/cm3), giving an effective density of 2.32 g/cm3.
#
# The lithium is enriched to 60 at% Li-6 (up from natural 7.5%) to enhance
# the tritium breeding ratio. The Li-6(n,t) reaction has a 1/v cross section
# that reaches ~940 barns at thermal energies, making Li-6 enrichment one
# of the most effective ways to increase TBR.
#
# Stoichiometry: Li4SiO4 -> 4 Li + 1 Si + 4 O
li4sio4 = openmc.Material(name="Li4SiO4 breeder (60% Li-6)")
li4sio4.add_nuclide("Li6", 4.0 * 0.60, "ao")   # 60% of 4 Li atoms enriched
li4sio4.add_nuclide("Li7", 4.0 * 0.40, "ao")   # 40% Li-7 remainder
li4sio4.add_element("Si", 1.0, "ao")             # 1 silicon per formula unit
li4sio4.add_element("O",  4.0, "ao")             # 4 oxygen per formula unit
li4sio4.set_density("g/cm3", 2.32)

# --- Beryllium neutron multiplier pebbles ---
# Beryllium pebbles serve as the neutron multiplier in the WCCB blanket.
# The key reaction is Be-9(n,2n), which converts one high-energy neutron
# into two lower-energy neutrons (threshold ~1.85 MeV, Q = -1.57 MeV).
# This multiplication is ESSENTIAL for achieving TBR > 1.
#
# Packing fraction: ~63% (random close packing of spherical pebbles)
# Bulk Be density: 1.85 g/cm3
# Effective pebble bed density: 1.85 * 0.63 = 1.17 g/cm3
be_pebbles = openmc.Material(name="Beryllium pebbles (multiplier)")
be_pebbles.add_nuclide("Be9", 1.0, "ao")   # 100% Be-9 (monoisotopic)
be_pebbles.set_density("g/cm3", 1.17)

# --- Water coolant (H2O) ---
# Pressurized water at ~300 K serves as the coolant in the WCCB blanket.
# Water also provides some neutron moderation, which enhances the thermal
# neutron flux available for Li-6(n,t) breeding.
water = openmc.Material(name="Water coolant (H2O)")
water.add_nuclide("H1", 2.0, "ao")
water.add_element("O", 1.0, "ao")
water.set_density("g/cm3", 1.0)
# Note: S(a,b) thermal scattering omitted because mix_materials() does not
# support materials with S(a,b) tables, and at fusion blanket temperatures
# (~600 K) the thermal scattering correction is less important.

# --- Homogenized WCCB blanket ---
# The breeding blanket is homogenized from its four constituent materials
# by volume fraction:
#   30% CLF-1 RAFM steel  (structural)
#   30% Li4SiO4 pebbles   (breeder)
#   30% Be pebbles         (multiplier)
#   10% H2O                (coolant)
#
# This homogenization is a standard simplification for scoping calculations.
# A heterogeneous model with distinct pebble beds and coolant channels would
# be more accurate but is beyond the scope of this CSG benchmark.
#
# Effective density: 0.30*7.8 + 0.30*2.32 + 0.30*1.17 + 0.10*1.0 = 3.487 g/cm3
blanket_homog = openmc.Material.mix_materials(
    [clf1_steel, li4sio4, be_pebbles, water],
    [0.30, 0.30, 0.30, 0.10],
    "vo",
    name="Homogenized WCCB blanket",
)

# --- SS316LN stainless steel ---
# Austenitic stainless steel used for the vacuum vessel, similar to ITER
# specification. SS316LN has controlled nitrogen content for improved
# mechanical properties at cryogenic temperatures.
#
# Composition (weight percent):
#   Fe: 62.5% (balance), Cr: 17.5%, Ni: 12.5%, Mo: 2.5%,
#   Mn: 2.0%, Si: 1.0%
# Density: 7.93 g/cm3
ss316 = openmc.Material(name="SS316LN stainless steel")
ss316.add_element("Fe", 0.625, "wo")   # balance
ss316.add_element("Cr", 0.175, "wo")   # 17.5% chromium
ss316.add_element("Ni", 0.125, "wo")   # 12.5% nickel
ss316.add_element("Mo", 0.025, "wo")   # 2.5% molybdenum
ss316.add_element("Mn", 0.020, "wo")   # 2% manganese
ss316.add_element("Si", 0.010, "wo")   # 1% silicon
ss316.set_density("g/cm3", 7.93)

# --- Back support structure (CLF-1 steel + water, 70/30 by volume) ---
# The back support structure sits behind the breeding blanket and provides
# mechanical support. It is a mixture of CLF-1 steel and water coolant
# channels at a 70/30 volume ratio.
bss_material = openmc.Material.mix_materials(
    [clf1_steel, water],
    [0.70, 0.30],
    "vo",
    name="Back support structure (CLF-1 + H2O)",
)

# --- Vacuum vessel fill (SS316 + water, 60/40 by volume) ---
# The vacuum vessel has a double-wall structure filled with a mixture of
# SS316 steel ribs and water for neutron shielding. The 60/40 ratio
# reflects the typical steel-to-water fraction in VV shielding blocks.
vv_fill_material = openmc.Material.mix_materials(
    [ss316, water],
    [0.60, 0.40],
    "vo",
    name="VV fill (SS316 + H2O)",
)

# --- SS304 stainless steel (thermal shield) ---
# The thermal shield between the vacuum vessel and the TF coils provides
# thermal insulation and radiation shielding for the superconducting magnets.
ss304 = openmc.Material(name="SS304 (thermal shield)")
ss304.add_element("Fe", 0.690, "wo")
ss304.add_element("Cr", 0.190, "wo")
ss304.add_element("Ni", 0.095, "wo")
ss304.add_element("Mn", 0.015, "wo")
ss304.add_element("Si", 0.010, "wo")
ss304.set_density("g/cm3", 7.93)

# --- TF coil casing (SS316 + Cu composite) ---
# Simplified representation of the toroidal field coil as a composite of
# SS316 steel casing and copper winding pack. In reality the TF coil
# contains Nb3Sn superconductor, but for shielding purposes the SS316+Cu
# approximation captures the dominant neutron interactions.
copper = openmc.Material(name="Copper (TF winding)")
copper.add_element("Cu", 1.0, "ao")
copper.set_density("g/cm3", 8.96)

tf_coil_material = openmc.Material.mix_materials(
    [ss316, copper],
    [0.50, 0.50],
    "vo",
    name="TF coil (SS316 + Cu composite)",
)

# --- Collect all materials ---
materials = openmc.Materials([
    clf1_steel, li4sio4, be_pebbles, water, blanket_homog,
    ss316, bss_material, vv_fill_material, ss304, copper, tf_coil_material,
])
materials.cross_sections = "/data/endfb-viii.0-hdf5/cross_sections.xml"


# =============================================================================
# Geometry
# =============================================================================
# The CFETR sector is built from concentric ZTorus surfaces centred at the
# origin. Each ZTorus has major radius R0=560 cm and circular cross-section
# (b = c = minor_radius). The sector is bounded by two reflective planes
# at phi=0 and phi=22.5 degrees.
#
# Using circular cross-sections (b = c) is a simplification; the actual
# plasma has elongation kappa=1.8. For a more accurate model, one would
# set b = kappa * c to create an elliptical cross-section, but this
# introduces complications with the radial build that are beyond the scope
# of this scoping benchmark.

# --- Toroidal surfaces (concentric ZTorus shells) ---
# ZTorus(x0, y0, z0, a, b, c) where:
#   a = major radius (distance from z-axis to tube centre)
#   b = minor radius in z-direction (parallel to axis of revolution)
#   c = minor radius in R-direction (perpendicular to axis of revolution)

# Plasma boundary (innermost surface)
torus_plasma = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_plasma, c=r_plasma,
    name="Plasma boundary (minor radius 130 cm)",
)

# First wall outer surface
torus_fw = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_fw, c=r_fw,
    name="First wall outer surface (minor radius 132 cm)",
)

# Blanket outer surface
torus_blanket = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_blanket, c=r_blanket,
    name="Blanket outer surface (minor radius 182 cm)",
)

# Back support structure outer surface
torus_bss = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_bss, c=r_bss,
    name="BSS outer surface (minor radius 197 cm)",
)

# Vacuum vessel inner wall outer surface
torus_vv_inner = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_vv_inner, c=r_vv_inner,
    name="VV inner wall outer surface (minor radius 201 cm)",
)

# Vacuum vessel fill outer surface
torus_vv_fill = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_vv_fill, c=r_vv_fill,
    name="VV fill outer surface (minor radius 221 cm)",
)

# Vacuum vessel outer wall outer surface
torus_vv_outer = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_vv_outer, c=r_vv_outer,
    name="VV outer wall outer surface (minor radius 225 cm)",
)

# Thermal shield outer surface
torus_shield = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_shield, c=r_shield,
    name="Thermal shield outer surface (minor radius 227 cm)",
)

# TF coil outer surface (outermost component)
torus_tf = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_tf, c=r_tf,
    name="TF coil outer surface (minor radius 252 cm)",
)

# --- Sector boundary planes ---
# The 22.5-degree sector is defined by two planes through the z-axis.
# Plane at phi=0: the y=0 plane (normal in -y direction)
# Plane at phi=22.5 deg: general plane Ax + By + Cz = D with D=0
#
# For the phi=0 boundary, particles on the +y side are inside the sector.
# For the phi=22.5 boundary, the plane normal points outward (away from
# the sector interior).
plane_phi0 = openmc.YPlane(
    y0=0.0,
    boundary_type="reflective",
    name="Sector boundary phi=0 (y=0 plane)",
)

# The phi=22.5 degree plane passes through the z-axis with normal
# perpendicular to the radial direction at that angle.
# A point at angle phi from the x-axis has direction (cos(phi), sin(phi), 0).
# The plane equation: sin(phi)*x - cos(phi)*y = 0
# This gives: A=sin(22.5), B=-cos(22.5), C=0, D=0
plane_phi22 = openmc.Plane(
    a=math.sin(phi_rad),
    b=-math.cos(phi_rad),
    c=0.0,
    d=0.0,
    boundary_type="reflective",
    name="Sector boundary phi=22.5 degrees",
)

# --- Vacuum boundary sphere ---
# A large sphere enclosing the entire sector. The outermost component
# (TF coil) extends to R0 + r_tf = 560 + 252 = 812 cm from the z-axis,
# and the highest point is at z = r_tf = 252 cm. A sphere of radius
# 900 cm centred at the origin is sufficient.
boundary_sphere = openmc.Sphere(
    x0=0.0, y0=0.0, z0=0.0,
    r=900.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Helper: sector region ---
# The toroidal sector is the region between the two boundary planes.
# For phi=0 plane (YPlane at y=0): the sector is on the +y side.
# For phi=22.5 plane: the sector is on the -halfspace side (toward phi=0).
def in_sector():
    """Return the angular region defining the 22.5-degree sector."""
    return +plane_phi0 & +plane_phi22


# =============================================================================
# Cell definitions
# =============================================================================
# Build cells from innermost (plasma) to outermost (TF coil), each bounded
# by consecutive torus surfaces and the sector planes.

cells = []

# --- Plasma chamber (void) ---
# The plasma region is treated as vacuum. In a real tokamak, this contains
# the burning D-T plasma, but for the fixed-source neutronics calculation,
# the plasma is replaced by a volumetric neutron source (defined below)
# and the region itself is void.
plasma_cell = openmc.Cell(
    name="Plasma chamber (void)",
    region=-torus_plasma & in_sector(),
)
cells.append(plasma_cell)

# --- First Wall (CLF-1 RAFM steel, 2 cm) ---
# The first wall is the plasma-facing component that directly receives the
# 14.1 MeV neutron flux. It must withstand extreme heat loads and neutron
# damage. CLF-1 steel is chosen for its low activation properties.
# The first wall causes some neutron scattering and parasitic capture,
# slightly reducing the flux reaching the blanket.
fw_cell = openmc.Cell(
    name="First wall (CLF-1, 2 cm)",
    fill=clf1_steel,
    region=+torus_plasma & -torus_fw & in_sector(),
)
cells.append(fw_cell)

# --- Breeding Blanket (homogenized WCCB, 50 cm) ---
# This is the MOST IMPORTANT component for tritium self-sufficiency.
# The homogenized blanket contains:
#   - Li4SiO4 breeder pebbles (30%): provides lithium for T breeding
#   - Beryllium pebbles (30%): neutron multiplication via Be-9(n,2n)
#   - CLF-1 steel structure (30%): mechanical support
#   - Water coolant (10%): heat removal
#
# The tritium breeding reactions occur primarily in the Li4SiO4:
#   Li-6 + n -> T + He-4 + 4.78 MeV  (dominant, large thermal sigma)
#   Li-7 + n -> T + He-4 + n' - 2.47 MeV  (secondary, threshold reaction)
#
# The beryllium multiplies neutrons via Be-9(n,2n), which is essential
# because each D-T fusion produces only ONE neutron, and some are
# inevitably lost to parasitic capture and leakage.
blanket_cell = openmc.Cell(
    name="Breeding blanket (WCCB homogenized, 50 cm)",
    fill=blanket_homog,
    region=+torus_fw & -torus_blanket & in_sector(),
)
cells.append(blanket_cell)

# --- Back Support Structure (CLF-1 + H2O, 15 cm) ---
# Provides mechanical support for the blanket modules and acts as an
# additional neutron shield for the vacuum vessel. The 70/30 steel/water
# mixture captures and moderates neutrons that penetrate the blanket.
bss_cell = openmc.Cell(
    name="Back support structure (CLF-1 + H2O, 15 cm)",
    fill=bss_material,
    region=+torus_blanket & -torus_bss & in_sector(),
)
cells.append(bss_cell)

# --- Vacuum Vessel inner wall (SS316, 4 cm) ---
# The inner structural wall of the double-walled vacuum vessel.
vv_inner_cell = openmc.Cell(
    name="VV inner wall (SS316, 4 cm)",
    fill=ss316,
    region=+torus_bss & -torus_vv_inner & in_sector(),
)
cells.append(vv_inner_cell)

# --- Vacuum Vessel fill (SS316 + H2O, 20 cm) ---
# The space between the VV double walls is filled with a mixture of
# steel ribs and water, providing both structural integrity and neutron
# shielding. This layer is critical for protecting the TF coils from
# neutron damage and nuclear heating.
vv_fill_cell = openmc.Cell(
    name="VV fill (SS316 + H2O, 20 cm)",
    fill=vv_fill_material,
    region=+torus_vv_inner & -torus_vv_fill & in_sector(),
)
cells.append(vv_fill_cell)

# --- Vacuum Vessel outer wall (SS316, 4 cm) ---
vv_outer_cell = openmc.Cell(
    name="VV outer wall (SS316, 4 cm)",
    fill=ss316,
    region=+torus_vv_fill & -torus_vv_outer & in_sector(),
)
cells.append(vv_outer_cell)

# --- Thermal Shield (SS304, 2 cm) ---
# An actively cooled thermal radiation shield that prevents heat from
# the warm vacuum vessel from reaching the cryogenic TF coils.
shield_cell = openmc.Cell(
    name="Thermal shield (SS304, 2 cm)",
    fill=ss304,
    region=+torus_vv_outer & -torus_shield & in_sector(),
)
cells.append(shield_cell)

# --- TF Coil casing (SS316 + Cu, 25 cm) ---
# Simplified representation of the toroidal field coil. The actual coil
# contains Nb3Sn superconducting cable-in-conduit conductors, but for
# neutronics purposes the SS316+Cu composite is adequate.
tf_cell = openmc.Cell(
    name="TF coil (SS316 + Cu composite, 25 cm)",
    fill=tf_coil_material,
    region=+torus_shield & -torus_tf & in_sector(),
)
cells.append(tf_cell)

# --- Void region ---
# Everything inside the vacuum boundary sphere but outside the tokamak
# components and sector boundaries.
tokamak_region = -torus_tf & in_sector()
void_cell = openmc.Cell(
    name="Void (surrounding space)",
    region=-boundary_sphere & ~tokamak_region,
)
cells.append(void_cell)

# --- Build geometry ---
root_universe = openmc.Universe(name="Root universe", cells=cells)
geometry = openmc.Geometry(root_universe)


# =============================================================================
# Source definition
# =============================================================================
# The D-T fusion neutron source is modelled as a uniform volumetric source
# distributed throughout the plasma region within the 22.5-degree sector.
# The source uses cylindrical coordinates (R, phi, z):
#   R: uniform from 450 to 670 cm (R0 +/- 110 cm, approximating the
#      plasma cross-section projected onto the midplane)
#   phi: uniform from 0 to 22.5 degrees (within the sector)
#   z: uniform from -100 to +100 cm (vertical extent, less than full
#      minor radius due to flux-weighted source profile)
#
# The energy is monoenergetic at 14.1 MeV, which is the neutron energy
# from the D + T -> He-4 + n reaction at rest. In reality, the D-T neutron
# energy has a slight spread due to ion temperature (~10-20 keV), but the
# monoenergetic approximation is standard for scoping studies.

source = openmc.IndependentSource()
r_dist = openmc.stats.Uniform(a=450.0, b=670.0)
phi_dist = openmc.stats.Uniform(a=0.0, b=phi_rad)
z_dist = openmc.stats.Uniform(a=-100.0, b=100.0)
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
settings.temperature = {"default": 600.0}   # 600 K average operating temperature
settings.output = {"tallies": True}


# =============================================================================
# Tallies
# =============================================================================
# Define tallies to extract the key neutronics quantities for CFETR:
#
#   1. Tritium Breeding Ratio (TBR): total triton production in the blanket
#      per source neutron. This is THE critical figure of merit -- TBR > 1.0
#      (ideally > 1.2 including margins) is required for tritium self-sufficiency.
#
#   2. Nuclear heating: energy deposition in each component. This determines
#      the power conversion efficiency and cooling requirements.
#
#   3. Be(n,2n) reaction rate: measures neutron multiplication efficiency
#      in the beryllium pebbles within the blanket.
#
#   4. Fast neutron flux at the vacuum vessel: determines radiation damage
#      rates (dpa) and component lifetime.

tallies = openmc.Tallies()

# --- Energy filter for fast flux ---
# Fast neutrons (E > 0.1 MeV = 1e5 eV) cause displacement damage in
# structural materials. The fast flux at the VV determines the component
# lifetime and replacement schedule.
fast_energy_filter = openmc.EnergyFilter([1.0e5, 20.0e6])

# --- Material filter for blanket TBR ---
blanket_mat_filter = openmc.MaterialFilter([blanket_homog])

# --- Cell filters for each component ---
blanket_cell_filter = openmc.CellFilter([blanket_cell])
fw_cell_filter = openmc.CellFilter([fw_cell])
bss_cell_filter = openmc.CellFilter([bss_cell])
vv_inner_cell_filter = openmc.CellFilter([vv_inner_cell])
vv_fill_cell_filter = openmc.CellFilter([vv_fill_cell])
vv_outer_cell_filter = openmc.CellFilter([vv_outer_cell])

# --- Tally 1: Tritium Breeding Ratio (TBR) ---
# The (n,Xt) score tallies total triton production from ALL reactions that
# produce tritons. In the homogenized blanket, this comes from:
#   Li-6(n,t)He-4 (MT=105): dominant, thermal/epithermal energies
#   Li-7(n,n't)He-4 (MT=112): secondary, threshold at 2.47 MeV
#
# The TBR is the total (n,Xt) score summed over the entire blanket,
# normalised per source neutron. For CFETR to be tritium self-sufficient,
# TBR must exceed ~1.05 (accounting for tritium losses in processing).
# The design target is TBR > 1.2 to provide margin.
tbr_tally = openmc.Tally(name="TBR")
tbr_tally.filters = [blanket_mat_filter]
tbr_tally.scores = ["(n,Xt)"]
tallies.append(tbr_tally)

# --- Tally 2: Nuclear heating in all major components ---
# The 'heating' score gives the total nuclear heating (neutron + photon
# kerma) deposited in each component per source neutron. This is needed
# for thermal-hydraulic design and power balance calculations.
component_cells = [
    fw_cell, blanket_cell, bss_cell,
    vv_inner_cell, vv_fill_cell, vv_outer_cell,
    shield_cell, tf_cell,
]
component_names = [
    "First wall", "Blanket", "Back support",
    "VV inner", "VV fill", "VV outer",
    "Thermal shield", "TF coil",
]
heating_cell_filter = openmc.CellFilter(component_cells)
heating_tally = openmc.Tally(name="nuclear_heating")
heating_tally.filters = [heating_cell_filter]
heating_tally.scores = ["heating"]
tallies.append(heating_tally)

# --- Tally 3: Be(n,2n) reaction rate in blanket ---
# The (n,2n) score in the blanket measures the neutron multiplication
# from beryllium. This is the key reaction enabling TBR > 1.
n2n_tally = openmc.Tally(name="Be_n2n_rate")
n2n_tally.filters = [blanket_cell_filter]
n2n_tally.scores = ["(n,2n)"]
tallies.append(n2n_tally)

# --- Tally 4: Fast neutron flux at vacuum vessel ---
# Fast flux (E > 0.1 MeV) at the VV determines displacement damage rates.
# The ITER limit for VV lifetime is ~1e22 n/cm2 (E > 0.1 MeV) over the
# design life. CFETR targets a similar constraint.
#
# We tally in the VV inner wall, VV fill, and VV outer wall separately.
vv_cells_filter = openmc.CellFilter([vv_inner_cell, vv_fill_cell, vv_outer_cell])
fast_flux_tally = openmc.Tally(name="fast_flux_VV")
fast_flux_tally.filters = [vv_cells_filter, fast_energy_filter]
fast_flux_tally.scores = ["flux"]
tallies.append(fast_flux_tally)

# --- Tally 5: Total neutron flux in blanket (for spectrum analysis) ---
# Energy-resolved flux in the blanket to assess neutron spectrum softening.
spectrum_energy_filter = openmc.EnergyFilter(
    np.logspace(np.log10(1.0e-5), np.log10(15.0e6), 201)
)
spectrum_tally = openmc.Tally(name="blanket_spectrum")
spectrum_tally.filters = [blanket_cell_filter, spectrum_energy_filter]
spectrum_tally.scores = ["flux"]
tallies.append(spectrum_tally)


# =============================================================================
# Build and export model
# =============================================================================
model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies,
)
model.export_to_xml()


# =============================================================================
# Weight window generation (optional, via --generate-ww flag)
# =============================================================================
# Weight windows improve the efficiency of fixed-source Monte Carlo transport
# by biasing the particle population toward regions of interest (typically the
# blanket and vacuum vessel). The FW-CADIS (Forward-Weighted Consistent
# Adjoint Driven Importance Sampling) method uses a deterministic adjoint
# calculation to generate spatially- and energy-dependent weight windows.
#
# This implementation uses OpenMC's random ray solver to perform the
# multigroup deterministic calculation needed for FW-CADIS.

if args.generate_ww:
    print("\n--- Generating weight windows (FW-CADIS) ---")

    model_rr = copy.deepcopy(model)
    model_rr.convert_to_multigroup(
        method="material_wise", groups="CASMO-4", nparticles=3000
    )
    model_rr.convert_to_random_ray()

    # Create a regular mesh covering the bounding box of the geometry
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
    print("Weight window model exported. Run OpenMC to generate weight windows.")
    print(f"  Mesh dimensions: ({nx}, {ny}, {nz})")
    print(f"  Lower left:  {list(bbox.lower_left)}")
    print(f"  Upper right: {list(bbox.upper_right)}")


# =============================================================================
# Summary
# =============================================================================
print("=" * 72)
print("CFETR 22.5-Degree Toroidal Sector -- OpenMC Model")
print("=" * 72)
print(f"  Major radius R0:    {R0} cm ({R0/100:.1f} m)")
print(f"  Plasma minor radius: {a_plasma} cm ({a_plasma/100:.1f} m)")
print(f"  Sector angle:        {sector_angle} degrees (1/{int(360/sector_angle)} "
      f"of full torus)")
print(f"  Radial build:")
print(f"    Plasma:            0 - {r_plasma} cm (void)")
print(f"    First wall:        {r_plasma} - {r_fw} cm (CLF-1, {t_fw} cm)")
print(f"    Blanket:           {r_fw} - {r_blanket} cm (WCCB, {t_blanket} cm)")
print(f"    Back support:      {r_blanket} - {r_bss} cm (CLF-1+H2O, {t_bss} cm)")
print(f"    VV inner wall:     {r_bss} - {r_vv_inner} cm (SS316, {t_vv_inner} cm)")
print(f"    VV fill:           {r_vv_inner} - {r_vv_fill} cm "
      f"(SS316+H2O, {t_vv_fill} cm)")
print(f"    VV outer wall:     {r_vv_fill} - {r_vv_outer} cm "
      f"(SS316, {t_vv_outer} cm)")
print(f"    Thermal shield:    {r_vv_outer} - {r_shield} cm "
      f"(SS304, {t_shield} cm)")
print(f"    TF coil:           {r_shield} - {r_tf} cm (SS316+Cu, {t_tf} cm)")
print(f"  Source:              14.1 MeV D-T, cylindrical volume")
print(f"                       R: 450-670 cm, phi: 0-{sector_angle} deg, "
      f"z: -100 to +100 cm")
print(f"  Temperature:         600 K (average)")
print(f"  Particles:           {args.particles:,} per batch x {args.batches} batches")
print(f"  Cross sections:      ENDF/B-VIII.0")
print(f"  Tallies:             TBR, nuclear heating, Be(n,2n), fast flux, "
      f"blanket spectrum")
print("=" * 72)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print(f"\nRunning simulation...")
model.run()
print("Simulation complete.")
