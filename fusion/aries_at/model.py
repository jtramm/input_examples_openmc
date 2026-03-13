#!/usr/bin/env python3
"""
ARIES-AT Advanced Tokamak Fusion Power Plant Benchmark
========================================================

Device
------
ARIES-AT (Advanced Tokamak) is a US design study for an advanced commercial
fusion power plant, conducted as part of the ARIES series of conceptual design
studies. The ARIES-AT design represents the most optimistic end of the tokamak
design space, assuming advanced physics and technology extrapolations to achieve
high performance in a compact device.

The ARIES-AT design features:
  - High plasma beta (~9.2%) enabled by advanced MHD stability regimes
  - Very high bootstrap current fraction (~91%), reducing external current
    drive power requirements to ~35 MW
  - SiC/LiPb blanket concept with SiC composite as the structural material
    and lithium-lead eutectic (Pb-17Li) as both breeder and coolant
  - Compact geometry (R = 5.2 m, a = 1.3 m) with 1755 MW fusion power
  - High outlet temperature (~1100 C) enabled by SiC, allowing high thermal
    efficiency (~59%) via a Brayton cycle

This model implements a simplified CSG toroidal sector of the ARIES-AT device,
capturing the essential radial build from the plasma through the first wall,
SiC/LiPb blanket, back wall/manifold, WC shield, vacuum vessel, and TF coil.

Radial Build (from plasma axis outward)
-----------------------------------------
The radial build uses concentric torus surfaces to define the layered structure.
All dimensions are measured from the plasma centre (minor radius a = 130 cm):

  Component          Inner r (cm)   Outer r (cm)   Thickness (cm)
  ---------          ------------   ------------   --------------
  Plasma (void)          0              130.0           130.0
  First wall (SiC)     130.0            130.4             0.4
  Blanket (SiC/LiPb)  130.4            160.4            30.0
  Back wall/manifold   160.4            170.4            10.0
  Shield (WC/H2O)      170.4            195.4            25.0
  Gap                  195.4            198.0             2.6
  Vacuum vessel (F82H) 198.0            213.0            15.0
  TF coil              213.0            243.0            30.0

The model uses a 22.5-degree toroidal sector (1/16 of the full torus,
corresponding to one of the 16 TF coils) with reflective boundary conditions
on the sector faces.

Blanket Physics
---------------
The ARIES-AT blanket is a self-cooled concept where the Pb-17Li eutectic
serves as both the tritium breeder and the primary coolant. The SiC composite
provides structural support with excellent high-temperature capability.

Key nuclear interactions:
  - Li-6(n,t)He-4: Primary tritium breeding reaction (sigma_th ~ 940 b).
    The Li is enriched to 90% Li-6 to maximize TBR.
  - Li-7(n,n't)He-4: Secondary breeding, threshold ~ 2.47 MeV
  - Pb-208(n,2n): Neutron multiplication in lead. This is the ARIES-AT
    equivalent of beryllium multiplication in solid breeder blankets.
    The Pb(n,2n) reaction has a threshold of ~7 MeV and converts one
    high-energy neutron into two lower-energy neutrons, boosting TBR.
  - Si-28(n,p), Si-28(n,alpha): Parasitic reactions in SiC that reduce TBR

The combination of Li-6 enrichment and Pb neutron multiplication gives
ARIES-AT a TBR of approximately 1.1, providing adequate tritium self-
sufficiency margin.

Source Definition
-----------------
The D-T fusion source is modelled as a uniform cylindrical source within
the plasma volume of a single 22.5-degree sector. The source neutron
energy is 14.1 MeV (mono-energetic D-T).

In reality, the fusion neutron source is peaked near the plasma centre
(approximately proportional to n^2 * <sigma*v>, which peaks on-axis in
an advanced tokamak with peaked pressure profiles). The uniform source
approximation is standard for neutronics scoping studies.

Materials
---------
  - SiC composite: 3.2 g/cm3, Si:C = 50:50 at%
  - Pb-17Li: 9.4 g/cm3, Pb:Li = 83:17 at%, Li enriched to 90% Li-6
  - Homogenized blanket: 40% SiC + 55% LiPb + 5% void
  - WC/H2O shield: 80% WC + 20% H2O by volume
  - F82H ferritic steel: 7.89 g/cm3 (reduced-activation steel for VV)
  - TF coil: homogenized SS316 + Cu mixture

Reference
---------
  F. Najmabadi et al., "The ARIES-AT advanced tokamak, Advanced technology
  fusion power plant," Fusion Engineering and Design, vol. 80, pp. 3-23, 2006.

  M.S. Tillack et al., "Fusion power core engineering for the ARIES-AT
  power plant," Fusion Engineering and Design, vol. 80, pp. 71-85, 2006.

  A.R. Raffray et al., "Engineering design and analysis of the ARIES-AT
  blanket," Fusion Engineering and Design, vol. 80, pp. 79-98, 2006.

  M.A. Abdou et al., "Blanket comparison and selection study -- final
  report," ANL/FPP-84-1, Argonne National Laboratory, 1984.

  M.E. Sawan, M.A. Abdou, "Physics and technology conditions for attaining
  tritium self-sufficiency for the DT fuel cycle," Fusion Engineering and
  Design, vol. 81, pp. 1131-1144, 2006.
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
    description="ARIES-AT Advanced Tokamak -- OpenMC fixed-source model"
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
    help="Generate weight windows using FW-CADIS (requires random ray support)",
)
args = parser.parse_args()


# =============================================================================
# Device parameters
# =============================================================================
# ARIES-AT major and minor radii define the plasma torus geometry.
# The major radius R0 is the distance from the machine axis to the centre
# of the plasma cross-section. The minor radius a is the half-width of the
# plasma (simplified as circular here; the actual elongated plasma has
# kappa = 2.2, but the ZTorus in OpenMC uses semi-axes b and c for the
# poloidal cross-section).
R0 = 520.0          # cm, major radius (5.2 m)
a_plasma = 130.0    # cm, minor radius (1.3 m)

# Toroidal field coil periodicity: 16 TF coils -> 22.5-degree sector
n_tf_coils = 16
sector_angle_deg = 360.0 / n_tf_coils    # 22.5 degrees
sector_angle_rad = math.radians(sector_angle_deg)

# Fusion power: 1755 MW
# Neutron wall loading: ~3.2 MW/m2 at outboard midplane
# These are not used in the transport calculation but provide context.

# Radial build thicknesses (cm)
first_wall_thickness = 0.4       # SiC first wall
blanket_thickness = 30.0         # SiC/LiPb breeding blanket
back_wall_thickness = 10.0       # SiC + LiPb back wall / manifold
shield_thickness = 25.0          # WC + water neutron shield
gap_thickness = 2.6              # vacuum gap between shield and VV
vv_thickness = 15.0              # F82H ferritic steel vacuum vessel
tf_coil_thickness = 30.0         # superconducting TF coil winding pack

# Compute radial positions (minor radius from plasma centre)
r_plasma = a_plasma                                          # 130.0
r_fw_outer = r_plasma + first_wall_thickness                 # 130.4
r_blanket_outer = r_fw_outer + blanket_thickness             # 160.4
r_backwall_outer = r_blanket_outer + back_wall_thickness     # 170.4
r_shield_outer = r_backwall_outer + shield_thickness         # 195.4
r_gap_outer = r_shield_outer + gap_thickness                 # 198.0
r_vv_outer = r_gap_outer + vv_thickness                      # 213.0
r_tf_outer = r_vv_outer + tf_coil_thickness                  # 243.0


# =============================================================================
# Materials
# =============================================================================
# Temperature for all materials: 800 K (SiC/LiPb blanket operates at
# approximately 700-1100 K; 800 K is a representative average)
material_temperature = 800.0  # K

# --- SiC composite (first wall and blanket structural material) ---
# Silicon carbide fibre-reinforced SiC matrix composite (SiC/SiC). The
# composite has ~90% of theoretical density (3.21 g/cm3 for single crystal).
# SiC is the primary structural material in ARIES-AT, chosen for its
# exceptional high-temperature strength, low activation, and compatibility
# with Pb-17Li at high temperatures.
#
# Nuclear properties: SiC is relatively transparent to neutrons compared
# to steel, which allows more neutrons to reach the lithium breeder.
# However, Si-28(n,p) and Si-28(n,alpha) are parasitic reactions that
# consume neutrons without producing tritium.
sic = openmc.Material(name="SiC composite")
sic.add_element("Si", 50.0, "ao")   # 50 atom% silicon
sic.add_element("C", 50.0, "ao")    # 50 atom% carbon
sic.set_density("g/cm3", 3.2)
sic.temperature = material_temperature

# --- Pb-17Li eutectic (breeder and coolant) ---
# The lead-lithium eutectic alloy Pb-17Li (17 atomic percent lithium in
# lead) serves as both the tritium breeder and the primary coolant in
# ARIES-AT. It flows through channels in the SiC blanket structure.
#
# The lithium is enriched to 90% Li-6 (from natural 7.5%) to maximise
# the tritium breeding ratio. At 90% enrichment:
#   Li-6 fraction: 0.90 * 17% = 15.3 at% of total alloy
#   Li-7 fraction: 0.10 * 17% = 1.7 at% of total alloy
#   Pb fraction: 83 at% (natural isotopic composition)
#
# Density at ~700 K operating temperature: 9.4 g/cm3
# (Pb-17Li has a melting point of ~235 C and is liquid during operation)
lipb = openmc.Material(name="Pb-17Li eutectic (90% Li-6)")
lipb.add_nuclide("Li6", 15.3, "ao")    # 90% of 17 at% Li
lipb.add_nuclide("Li7", 1.7, "ao")     # 10% of 17 at% Li
lipb.add_element("Pb", 83.0, "ao")     # 83 at% natural lead
lipb.set_density("g/cm3", 9.4)
lipb.temperature = material_temperature

# --- Homogenized SiC/LiPb blanket material ---
# The blanket is a composite of SiC structural elements (flow channel
# inserts, first wall backing) and flowing Pb-17Li breeder/coolant.
# Volume fractions from the ARIES-AT design:
#   40% SiC structure
#   55% Pb-17Li breeder/coolant
#    5% void (gaps, tolerances, helium purge gas channels)
#
# Effective density: 0.40 * 3.2 + 0.55 * 9.4 + 0.05 * 0.0 = 6.45 g/cm3
#
# We create this as an openmc.Material.mix_materials() blend.
blanket_mat = openmc.Material.mix_materials(
    [sic, lipb],
    [0.40, 0.55],     # volume fractions (5% void omitted, renormalize below)
    percent_type="vo",
    name="SiC/LiPb blanket (40% SiC, 55% LiPb, 5% void)",
)
# Adjust density to account for 5% void: multiply by 0.95/0.95 is already
# handled by the volume fractions summing to 0.95. mix_materials normalizes
# internally, so we manually set the correct effective density.
blanket_mat.set_density("g/cm3", 6.45)
blanket_mat.temperature = material_temperature

# --- Back wall / manifold material ---
# The back wall and manifold region behind the blanket uses the same
# SiC/LiPb material system but with a higher SiC fraction for structural
# support of the piping manifolds. We approximate it with the same
# homogenized blanket material for simplicity.
backwall_mat = openmc.Material.mix_materials(
    [sic, lipb],
    [0.50, 0.50],
    percent_type="vo",
    name="Back wall/manifold (SiC/LiPb)",
)
backwall_mat.set_density("g/cm3", 0.50 * 3.2 + 0.50 * 9.4)
backwall_mat.temperature = material_temperature

# --- WC/H2O shield material ---
# Tungsten carbide is an extremely effective neutron shield due to
# tungsten's high atomic number (Z=74) and large inelastic scattering
# cross section. Water provides hydrogen for neutron moderation.
# The combination is standard for fusion reactor shielding:
#   80% WC by volume (density 15.63 g/cm3)
#   20% H2O by volume (density 1.0 g/cm3)
#   Effective density: 0.80 * 15.63 + 0.20 * 1.0 = 12.7 g/cm3
#
# WC composition: W:C = 50:50 atomic ratio
wc = openmc.Material(name="Tungsten carbide (WC)")
wc.add_element("W", 50.0, "ao")
wc.add_element("C", 50.0, "ao")
wc.set_density("g/cm3", 15.63)

water = openmc.Material(name="Water (H2O)")
water.add_element("H", 2.0, "ao")
water.add_element("O", 1.0, "ao")
water.set_density("g/cm3", 1.0)
# Note: S(a,b) thermal scattering omitted because mix_materials() does not
# support materials with S(a,b) tables. At 800 K operating temperature,
# the thermal scattering correction is negligible.

shield_mat = openmc.Material.mix_materials(
    [wc, water],
    [0.80, 0.20],
    percent_type="vo",
    name="WC/H2O shield (80% WC, 20% H2O)",
)
shield_mat.set_density("g/cm3", 12.7)
shield_mat.temperature = material_temperature

# --- F82H reduced-activation ferritic/martensitic steel (vacuum vessel) ---
# F82H (8Cr-2W-VTa) is a reduced-activation ferritic/martensitic steel
# developed for fusion applications. It is the reference structural
# material for the ARIES-AT vacuum vessel. The reduced-activation
# composition (W instead of Mo, V and Ta instead of Nb) minimises
# long-lived activation products, enabling simpler waste management.
#
# Composition (weight percent):
#   Fe: 89.0% (balance), Cr: 7.5%, W: 2.0%, V: 0.2%,
#   Ta: 0.04%, Mn: 0.1%, Si: 0.1%, C: 0.1%
# Density: 7.89 g/cm3
f82h = openmc.Material(name="F82H ferritic steel (vacuum vessel)")
f82h.add_element("Fe", 89.0, "wo")     # balance
f82h.add_element("Cr", 7.5, "wo")      # 7.5 wt% chromium
f82h.add_element("W", 2.0, "wo")       # 2.0 wt% tungsten
f82h.add_element("V", 0.2, "wo")       # 0.2 wt% vanadium
f82h.add_element("Ta", 0.04, "wo")     # 0.04 wt% tantalum
f82h.add_element("Mn", 0.1, "wo")      # 0.1 wt% manganese
f82h.add_element("Si", 0.16, "wo")     # 0.16 wt% silicon (adjusted for balance)
f82h.add_element("C", 0.1, "wo")       # 0.1 wt% carbon
f82h.set_density("g/cm3", 7.89)
f82h.temperature = material_temperature

# --- TF coil winding pack (homogenized) ---
# The ARIES-AT TF coils use Nb3Sn superconductor with copper stabilizer
# in a stainless steel conduit. The winding pack is homogenized as:
#   40% SS316 (structural conduit and case)
#   30% Copper (stabilizer)
#   20% Nb3Sn (superconductor)
#   10% epoxy (insulation)
# We approximate this as an SS316 + Cu mixture for simplicity, since the
# Nb3Sn and epoxy fractions are small and the TF coil is behind heavy
# shielding.
ss316 = openmc.Material(name="SS316 stainless steel")
ss316.add_element("Fe", 65.42, "wo")
ss316.add_element("Cr", 17.0, "wo")
ss316.add_element("Ni", 12.0, "wo")
ss316.add_element("Mo", 2.5, "wo")
ss316.add_element("Mn", 2.0, "wo")
ss316.add_element("Si", 1.0, "wo")
ss316.add_element("C", 0.08, "wo")
ss316.set_density("g/cm3", 7.93)

copper = openmc.Material(name="Copper (TF coil stabilizer)")
copper.add_element("Cu", 100.0, "ao")
copper.set_density("g/cm3", 8.96)

tf_coil_mat = openmc.Material.mix_materials(
    [ss316, copper],
    [0.60, 0.40],     # approximate 60% steel, 40% copper by volume
    percent_type="vo",
    name="TF coil winding pack (SS316/Cu)",
)
tf_coil_mat.set_density("g/cm3", 0.60 * 7.93 + 0.40 * 8.96)
tf_coil_mat.temperature = material_temperature

# --- Collect all materials ---
materials = openmc.Materials([
    sic, lipb, blanket_mat, backwall_mat,
    shield_mat, f82h, tf_coil_mat,
])
materials.cross_sections = "/data/endfb-viii.0-hdf5/cross_sections.xml"
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# The geometry uses concentric ZTorus surfaces to define the radial build,
# bounded by two planes forming a 22.5-degree toroidal sector with
# reflective boundary conditions.
#
# ZTorus parameters:
#   x0, y0, z0 = 0, 0, 0 (torus centred at origin)
#   a = R0 = 520 cm (major radius)
#   b = minor radius along z (vertical, = minor radius for circular cross-section)
#   c = minor radius in R-direction (horizontal, = minor radius for circular)
#
# Note: For ARIES-AT, the plasma has elongation kappa = 2.2, meaning
# b_plasma = kappa * a = 2.2 * 130 = 286 cm for the actual plasma shape.
# However, the blanket and shield components have a more circular cross-
# section in practice. We use circular (b = c = r) for all components
# to keep the CSG geometry tractable.

# --- Torus surfaces for each radial layer ---
# Plasma boundary
torus_plasma = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_plasma, c=r_plasma,
    name="Plasma boundary (a = 130 cm)",
)

# First wall outer surface
torus_fw = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_fw_outer, c=r_fw_outer,
    name=f"First wall outer (r = {r_fw_outer} cm)",
)

# Blanket outer surface
torus_blanket = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_blanket_outer, c=r_blanket_outer,
    name=f"Blanket outer (r = {r_blanket_outer} cm)",
)

# Back wall outer surface
torus_backwall = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_backwall_outer, c=r_backwall_outer,
    name=f"Back wall outer (r = {r_backwall_outer} cm)",
)

# Shield outer surface
torus_shield = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_shield_outer, c=r_shield_outer,
    name=f"Shield outer (r = {r_shield_outer} cm)",
)

# Gap outer surface (vacuum vessel inner)
torus_gap = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_gap_outer, c=r_gap_outer,
    name=f"Gap outer / VV inner (r = {r_gap_outer} cm)",
)

# Vacuum vessel outer surface
torus_vv = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_vv_outer, c=r_vv_outer,
    name=f"Vacuum vessel outer (r = {r_vv_outer} cm)",
)

# TF coil outer surface
torus_tf = openmc.ZTorus(
    x0=0.0, y0=0.0, z0=0.0,
    a=R0, b=r_tf_outer, c=r_tf_outer,
    name=f"TF coil outer (r = {r_tf_outer} cm)",
)

# --- Sector boundary planes ---
# The 22.5-degree sector is bounded by two planes passing through the z-axis.
# Plane 1: y = 0 (the xz-plane)
# Plane 2: rotated 22.5 degrees from the xz-plane
#
# A general plane through the z-axis has the form: a*x + b*y = 0
# For the y=0 plane: a=0, b=1
# For the rotated plane at angle phi: a = sin(phi), b = -cos(phi)
# (This gives a plane whose normal points into the sector.)
phi = sector_angle_rad
plane_0 = openmc.YPlane(y0=0.0, boundary_type="reflective",
                         name="Sector boundary at phi=0")
plane_22 = openmc.Plane(a=math.sin(phi), b=-math.cos(phi), c=0.0, d=0.0,
                         boundary_type="reflective",
                         name=f"Sector boundary at phi={sector_angle_deg} deg")

# --- Bounding box ---
# The outer boundary encloses the entire TF coil torus region plus margin.
# We use a large sphere centred at the origin.
r_outer_max = R0 + r_tf_outer + 50.0   # generous margin
bounding_sphere = openmc.Sphere(
    x0=0.0, y0=0.0, z0=0.0,
    r=r_outer_max,
    boundary_type="vacuum",
    name="Outer vacuum boundary",
)

# --- Helper: sector region (between the two planes) ---
# The sector lies where y > 0 AND below the rotated plane.
# For the half-space signs:
#   +plane_0 means y > 0
#   -plane_22 means sin(phi)*x - cos(phi)*y < 0
sector_region = +plane_0 & -plane_22


# =============================================================================
# Cells
# =============================================================================
cells = []

# --- Plasma (void) ---
# The plasma region is modelled as void. In reality, the D-T plasma has
# negligible density (~1e-5 g/cm3) for neutron transport purposes.
plasma_cell = openmc.Cell(
    name="Plasma (void)",
    region=-torus_plasma & sector_region,
)
cells.append(plasma_cell)

# --- First wall (SiC composite, 0.4 cm) ---
# The first wall is the plasma-facing component. It must withstand the
# neutron wall loading (~3.2 MW/m2) and surface heat flux. SiC composite
# is chosen for its high-temperature capability and low activation.
fw_cell = openmc.Cell(
    name="First wall (SiC, 0.4 cm)",
    fill=sic,
    region=+torus_plasma & -torus_fw & sector_region,
)
cells.append(fw_cell)

# --- Blanket (SiC/LiPb, 30 cm) ---
# The breeding blanket is the region where tritium is produced from
# Li-6(n,t) and Li-7(n,n't) reactions in the flowing Pb-17Li.
# Lead provides neutron multiplication via Pb(n,2n) reactions.
# This is the most neutronically important region.
blanket_cell = openmc.Cell(
    name="Blanket (SiC/LiPb, 30 cm)",
    fill=blanket_mat,
    region=+torus_fw & -torus_blanket & sector_region,
)
cells.append(blanket_cell)

# --- Back wall / manifold (SiC/LiPb, 10 cm) ---
# Structural region behind the blanket containing flow manifolds for
# the Pb-17Li coolant. Higher SiC fraction than the blanket.
backwall_cell = openmc.Cell(
    name="Back wall/manifold (SiC/LiPb, 10 cm)",
    fill=backwall_mat,
    region=+torus_blanket & -torus_backwall & sector_region,
)
cells.append(backwall_cell)

# --- Shield (WC/H2O, 25 cm) ---
# The neutron shield uses tungsten carbide mixed with water to attenuate
# the neutron flux by several orders of magnitude before it reaches the
# vacuum vessel and TF coils. WC is extremely effective due to tungsten's
# large inelastic scattering cross section and high density.
shield_cell = openmc.Cell(
    name="Shield (WC/H2O, 25 cm)",
    fill=shield_mat,
    region=+torus_backwall & -torus_shield & sector_region,
)
cells.append(shield_cell)

# --- Gap (void, 2.6 cm) ---
# Vacuum gap between shield and vacuum vessel for thermal insulation
# and assembly tolerances.
gap_cell = openmc.Cell(
    name="Gap (void, 2.6 cm)",
    region=+torus_shield & -torus_gap & sector_region,
)
cells.append(gap_cell)

# --- Vacuum vessel (F82H, 15 cm) ---
# The vacuum vessel provides the primary vacuum boundary and structural
# support. F82H reduced-activation ferritic/martensitic steel is used
# for its favourable activation characteristics.
vv_cell = openmc.Cell(
    name="Vacuum vessel (F82H, 15 cm)",
    fill=f82h,
    region=+torus_gap & -torus_vv & sector_region,
)
cells.append(vv_cell)

# --- TF coil (SS316/Cu, 30 cm) ---
# The superconducting TF coil winding pack. Neutron flux and nuclear
# heating at the TF coil are critical design parameters -- the coils
# must be protected to maintain superconducting operation and limit
# radiation damage to the insulation.
tf_cell = openmc.Cell(
    name="TF coil (SS316/Cu, 30 cm)",
    fill=tf_coil_mat,
    region=+torus_vv & -torus_tf & sector_region,
)
cells.append(tf_cell)

# --- Void outside TF coil, within sector ---
outer_void_cell = openmc.Cell(
    name="Void (outside TF coil)",
    region=+torus_tf & -bounding_sphere & sector_region,
)
cells.append(outer_void_cell)

# --- Void outside sector (rest of bounding sphere) ---
# Everything inside the bounding sphere but outside the sector
outside_sector_cell = openmc.Cell(
    name="Void (outside sector)",
    region=-bounding_sphere & ~sector_region,
)
cells.append(outside_sector_cell)

# --- Build root universe and geometry ---
root_universe = openmc.Universe(name="Root universe", cells=cells)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The D-T source is modelled as a uniform cylindrical source within the
# plasma volume of the 22.5-degree sector. The source occupies a cylindrical
# annulus centred on the z-axis:
#   R: from R0 - 0.77*a to R0 + 0.77*a (covers ~77% of the minor radius,
#      representing the region where most fusion reactions occur)
#   phi: from 0 to 22.5 degrees (matching the sector)
#   z: from -100 to +100 cm (covers most of the plasma height)
#
# Energy: mono-energetic 14.1 MeV (D-T fusion neutrons)
# Angular distribution: isotropic (fusion neutrons are emitted isotropically
# in the centre-of-mass frame, and the thermal ion velocities are negligible
# compared to the neutron velocity)

source = openmc.IndependentSource()
r_dist = openmc.stats.Uniform(a=420.0, b=620.0)
phi_dist = openmc.stats.Uniform(a=0.0, b=sector_angle_rad)
z_dist = openmc.stats.Uniform(a=-100.0, b=100.0)
source.space = openmc.stats.CylindricalIndependent(
    r=r_dist, phi=phi_dist, z=z_dist,
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
settings.photon_transport = True     # track photons for nuclear heating
settings.output = {"tallies": True}

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# Key neutronics parameters for a fusion power plant assessment:
#   1. Tritium Breeding Ratio (TBR): must exceed ~1.05 for tritium self-
#      sufficiency (accounting for losses in the tritium processing system)
#   2. Nuclear heating in each component (for thermal-hydraulic design)
#   3. Energy multiplication factor (ratio of total deposited energy to
#      source neutron energy; measures the contribution of exothermic
#      nuclear reactions like Li-6(n,t) and captures)
#   4. Neutron flux at the TF coil (for radiation damage assessment and
#      superconductor lifetime)

tallies = openmc.Tallies()

# --- Tally 1: Tritium Breeding Ratio (TBR) ---
# (n,Xt) scores total triton production from all reactions.
# In the SiC/LiPb blanket, tritium is produced by:
#   Li-6(n,t)He-4: dominant (sigma_th ~ 940 b, 1/v)
#   Li-7(n,n't)He-4: secondary (threshold 2.47 MeV)
# The back wall also contributes some tritium production.
blanket_filter = openmc.CellFilter([blanket_cell])
backwall_filter = openmc.CellFilter([backwall_cell])
breeding_filter = openmc.CellFilter([blanket_cell, backwall_cell])

tbr_tally = openmc.Tally(name="TBR")
tbr_tally.filters = [breeding_filter]
tbr_tally.scores = ["(n,Xt)"]
tallies.append(tbr_tally)

# Separate TBR tallies for blanket and back wall
tbr_blanket = openmc.Tally(name="TBR_blanket")
tbr_blanket.filters = [blanket_filter]
tbr_blanket.scores = ["(n,Xt)"]
tallies.append(tbr_blanket)

tbr_backwall = openmc.Tally(name="TBR_backwall")
tbr_backwall.filters = [backwall_filter]
tbr_backwall.scores = ["(n,Xt)"]
tallies.append(tbr_backwall)

# --- Tally 2: Nuclear heating in each component ---
# The "heating" score gives the energy deposited per source particle
# in each cell, including contributions from neutron reactions,
# gamma-ray interactions, and charged particle deposition.
all_component_cells = [
    fw_cell, blanket_cell, backwall_cell,
    shield_cell, vv_cell, tf_cell,
]
component_names = [
    "First wall", "Blanket", "Back wall",
    "Shield", "Vacuum vessel", "TF coil",
]

for cell, name in zip(all_component_cells, component_names):
    heating_tally = openmc.Tally(name=f"heating_{name}")
    heating_tally.filters = [openmc.CellFilter([cell])]
    heating_tally.scores = ["heating"]
    tallies.append(heating_tally)

# Total heating across all components
total_heating = openmc.Tally(name="heating_total")
total_heating.filters = [openmc.CellFilter(all_component_cells)]
total_heating.scores = ["heating"]
tallies.append(total_heating)

# --- Tally 3: Neutron flux at TF coil ---
# The neutron flux spectrum at the TF coil is critical for assessing:
#   - Radiation damage to the superconductor (dose to insulation)
#   - Nuclear heating in the cryogenic coil (must stay below ~10-20 kW)
#   - Nb3Sn critical current degradation from fast neutron fluence
energy_bins = np.logspace(
    np.log10(1.0e-3),     # 1 meV (thermal)
    np.log10(15.0e6),      # 15 MeV
    101,                    # 101 edges -> 100 bins
)
energy_filter = openmc.EnergyFilter(energy_bins)

tf_filter = openmc.CellFilter([tf_cell])
tf_flux = openmc.Tally(name="flux_TF_coil")
tf_flux.filters = [tf_filter, energy_filter]
tf_flux.scores = ["flux"]
tallies.append(tf_flux)

# Total flux at TF coil (energy-integrated)
tf_total_flux = openmc.Tally(name="flux_TF_coil_total")
tf_total_flux.filters = [tf_filter]
tf_total_flux.scores = ["flux"]
tallies.append(tf_total_flux)

# --- Tally 4: Neutron flux in blanket (for NWL estimate) ---
blanket_flux = openmc.Tally(name="flux_blanket")
blanket_flux.filters = [blanket_filter]
blanket_flux.scores = ["flux"]
tallies.append(blanket_flux)

# --- Tally 5: Pb(n,2n) neutron multiplication in blanket ---
# Lead neutron multiplication is the key mechanism boosting TBR in
# liquid metal blankets (analogous to Be(n,2n) in solid breeder blankets).
pb_n2n_tally = openmc.Tally(name="Pb_n2n_blanket")
pb_n2n_tally.filters = [breeding_filter]
pb_n2n_tally.scores = ["(n,2n)"]
tallies.append(pb_n2n_tally)

tallies.export_to_xml()


# =============================================================================
# Weight window generation (optional)
# =============================================================================
# The WC shield is extremely dense (12.7 g/cm3) and attenuates the neutron
# flux by many orders of magnitude. Weight windows are essential for
# getting meaningful statistics at the TF coil location.
#
# This uses the FW-CADIS (Forward-Weighted Consistent Adjoint Driven
# Importance Sampling) method via OpenMC's random ray solver to generate
# spatially dependent weight windows.

def generate_weight_windows(model):
    """Generate weight windows using FW-CADIS with random ray.

    Parameters
    ----------
    model : openmc.Model
        The OpenMC model to generate weight windows for.
    """
    model_rr = copy.deepcopy(model)
    model_rr.convert_to_multigroup(
        method="material_wise", groups="CASMO-4", nparticles=3000,
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
    print("Re-run the main model to use them.")


if args.weight_windows:
    print("Generating weight windows via FW-CADIS random ray...")
    model = openmc.Model(
        geometry=geometry, materials=materials, settings=settings, tallies=tallies,
    )
    generate_weight_windows(model)
else:
    # ==========================================================================
    # Summary
    # ==========================================================================
    print("=" * 72)
    print("ARIES-AT Advanced Tokamak -- OpenMC Model")
    print("=" * 72)
    print(f"  Major radius R0:    {R0} cm ({R0/100:.1f} m)")
    print(f"  Minor radius a:     {a_plasma} cm ({a_plasma/100:.1f} m)")
    print(f"  Sector angle:       {sector_angle_deg} deg (1/{n_tf_coils} of torus)")
    print(f"  Fusion power:       1755 MW")
    print(f"  NWL (outboard):     ~3.2 MW/m2")
    print()
    print("  Radial build (minor radius from plasma centre):")
    print(f"    Plasma (void):       0 - {r_plasma:.1f} cm")
    print(f"    First wall (SiC):    {r_plasma:.1f} - {r_fw_outer:.1f} cm "
          f"({first_wall_thickness:.1f} cm)")
    print(f"    Blanket (SiC/LiPb):  {r_fw_outer:.1f} - {r_blanket_outer:.1f} cm "
          f"({blanket_thickness:.1f} cm)")
    print(f"    Back wall:           {r_blanket_outer:.1f} - {r_backwall_outer:.1f} cm "
          f"({back_wall_thickness:.1f} cm)")
    print(f"    Shield (WC/H2O):     {r_backwall_outer:.1f} - {r_shield_outer:.1f} cm "
          f"({shield_thickness:.1f} cm)")
    print(f"    Gap (void):          {r_shield_outer:.1f} - {r_gap_outer:.1f} cm "
          f"({gap_thickness:.1f} cm)")
    print(f"    Vacuum vessel:       {r_gap_outer:.1f} - {r_vv_outer:.1f} cm "
          f"({vv_thickness:.1f} cm)")
    print(f"    TF coil:             {r_vv_outer:.1f} - {r_tf_outer:.1f} cm "
          f"({tf_coil_thickness:.1f} cm)")
    print()
    print("  Materials:")
    print(f"    SiC composite:     {sic.density:.2f} g/cm3")
    print(f"    Pb-17Li (90% Li6): {lipb.density:.2f} g/cm3")
    print(f"    Blanket (homog.):  {blanket_mat.density:.2f} g/cm3")
    print(f"    WC/H2O shield:     {shield_mat.density:.2f} g/cm3")
    print(f"    F82H steel (VV):   {f82h.density:.2f} g/cm3")
    print(f"    TF coil (homog.):  {tf_coil_mat.density:.2f} g/cm3")
    print()
    print(f"  Temperature:   {material_temperature:.0f} K")
    print(f"  Cross sections: ENDF/B-VIII.0")
    print(f"  Particles:     {args.particles:,} per batch x {args.batches} batches "
          f"= {args.particles * args.batches:,} total")
    print()
    print("  Tallies:")
    print("    - TBR: tritium production (n,Xt) in blanket + back wall")
    print("    - Nuclear heating in all components")
    print("    - Neutron flux spectrum at TF coil")
    print("    - Pb(n,2n) neutron multiplication rate")
    print()
    print("  Use --weight-windows to generate variance reduction for deep shielding.")
    print("=" * 72)
    print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
    print("Run with: openmc")
