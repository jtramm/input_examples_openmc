#!/usr/bin/env python3
"""
FNG HCPB Tritium Breeder Module Mock-up Benchmark
====================================================

Facility
--------
The Frascati Neutron Generator (FNG) is located at ENEA Frascati, near Rome,
Italy. FNG produces 14.1 MeV neutrons via the D-T (deuterium-tritium) fusion
reaction using a 300 keV deuteron beam impinging on a tritium target. The
facility has been a principal European facility for integral validation of
nuclear data relevant to ITER and future fusion power plants, with a maximum
yield of approximately 1e11 neutrons per second.

Experiment
----------
This benchmark models a mock-up of the Helium-Cooled Pebble Bed (HCPB)
breeding blanket concept from the SINBAD database. The HCPB is one of the
two candidate Test Blanket Module (TBM) designs being developed by the
European Union for testing in ITER, alongside the Helium-Cooled Lithium-Lead
(HCLL) concept. In the HCPB design, tritium is bred in lithium-containing
ceramic pebbles (Li4SiO4 or Li2TiO3) while beryllium pebbles serve as the
neutron multiplier. The FNG mock-up uses Li2CO3 powder as the breeder
surrogate and metallic beryllium blocks as the neutron multiplier.

The mock-up consists of:
  - A main AISI-303 stainless steel box (31.0 x 29.0 x 30.9 cm external)
    filled with metallic beryllium and containing two double-layer breeder
    sections of Li2CO3 powder
  - A rear AISI-316 stainless steel cassette (31.0 x 14.8 x 30.9 cm external)
    filled entirely with Li2CO3 powder

The two double-layer breeder sections within the main box each consist of
two 1.2 cm layers of Li2CO3 powder separated by a 1 mm thick stainless
steel divider wall, giving 2.5 cm per double-layer section.

Tritium Breeding Physics
------------------------
In a fusion power plant, the D-T reaction produces one neutron per fusion
event. Since tritium does not occur naturally in sufficient quantities, the
blanket must breed tritium from lithium to close the fuel cycle. The two
principal breeding reactions are:

  Li-6 + n -> He-4 + T + 4.78 MeV   (exothermic, large thermal cross section)
  Li-7 + n -> He-4 + T + n' - 2.47 MeV  (endothermic, threshold ~2.5 MeV)

The Li-6(n,t) reaction has a very large cross section at thermal energies
(~940 barns at 0.025 eV) following a 1/v law, making it the dominant
breeding reaction. However, since each D-T fusion event produces only ONE
neutron, and some neutrons are inevitably lost to parasitic absorption and
leakage, a neutron MULTIPLIER is needed to achieve a Tritium Breeding Ratio
(TBR) greater than unity.

Beryllium is the neutron multiplier of choice for the HCPB concept. The
key reaction is:

  Be-9 + n -> 2 He-4 + 2n - 1.57 MeV  (threshold ~1.85 MeV)

This (n,2n) reaction converts one high-energy neutron into two lower-energy
neutrons, effectively doubling the neutron population available for breeding.
The 14.1 MeV D-T neutrons have more than enough energy to drive this
reaction, and the resulting moderated neutrons are efficiently captured by
Li-6 to produce tritium.

The interplay between beryllium multiplication and lithium breeding is the
FUNDAMENTAL physics being validated by this benchmark. The accuracy of the
Be-9(n,2n) cross section, the Li-6(n,t) cross section, and neutron
slowing-down in beryllium directly determine the predicted TBR.

Geometry
--------
  Coordinate system:
    x-axis: horizontal transverse (width of assembly)
    y-axis: beam axis (depth into assembly, neutrons travel in +y)
    z-axis: vertical

  Assembly layout along beam axis (y):
    y = -5.3 cm     : D-T point source position
    y = 0.0 cm      : Front face of main SS box (0.5 cm AISI-303 wall)
    y = 0.5 cm      : Start of internal fill region
    y = 0.5 - 4.0 cm: Beryllium zone 1
    y = 4.0 - 4.15 cm: SS wall (0.15 cm, cassette wall of first breeder)
    y = 4.15 - 5.35 cm: Li2CO3 breeder layer 1a (1.2 cm)
    y = 5.35 - 5.45 cm: SS divider wall (0.1 cm = 1 mm)
    y = 5.45 - 6.65 cm: Li2CO3 breeder layer 1b (1.2 cm)
    y = 6.65 - 6.80 cm: SS wall (0.15 cm, cassette wall)
    y = 6.80 - 16.3 cm: Beryllium zone 2 (large central block)
    y = 16.3 - 16.45 cm: SS wall (0.15 cm, cassette wall of second breeder)
    y = 16.45 - 17.65 cm: Li2CO3 breeder layer 2a (1.2 cm)
    y = 17.65 - 17.75 cm: SS divider wall (0.1 cm = 1 mm)
    y = 17.75 - 18.95 cm: Li2CO3 breeder layer 2b (1.2 cm)
    y = 18.95 - 19.10 cm: SS wall (0.15 cm, cassette wall)
    y = 19.10 - 28.5 cm: Beryllium zone 3
    y = 28.5 cm     : Rear face of main box (0.5 cm wall -> y = 29.0 cm)
    y = 29.0 cm     : Start of rear cassette (0.5 cm AISI-316 wall)
    y = 29.5 cm     : Start of rear Li2CO3 fill
    y = 43.3 cm     : End of rear Li2CO3 fill (0.5 cm wall -> y = 43.8 cm)
    y = 43.8 cm     : Rear face of rear cassette

  Note: The main box external y-dimension is 29.0 cm (front wall at y=0,
  rear wall ending at y=29.0). The rear cassette external y-dimension is
  14.8 cm (from y=29.0 to y=43.8).

  The detector measurement positions are at y = 4.2, 10.5, 16.8, and
  23.1 cm from the front surface, which fall within:
    - y = 4.2 cm  -> first breeder double-layer region
    - y = 10.5 cm -> central beryllium region (neutron spectrum measurement)
    - y = 16.8 cm -> second breeder double-layer region
    - y = 23.1 cm -> rear beryllium region

Materials
---------
  - Beryllium metal: density 1.85 g/cm3
  - Li2CO3 powder (main box): density 1.123 g/cm3, natural Li (7.5% Li-6)
  - Li2CO3 powder (rear cassette): density 0.9413 g/cm3, natural Li
  - AISI-303 SS: density 7.954 g/cm3
      Fe ~69%, Cr ~18%, Ni ~9%, Mn ~2%, Si ~1%, C ~0.08%, S ~0.25%
  - AISI-316 SS: density 7.99 g/cm3
      Fe ~65%, Cr ~17%, Ni ~12%, Mo ~2.5%, Mn ~2%, Si ~1%, C ~0.08%

Reference
---------
  SINBAD database: "FNG HCPB Tritium Breeder Module Mock-up"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_hcpb/fnghcpb-a.htm

  P. Batistoni et al., "Neutronics experiment on a mock-up of the
  Helium-Cooled Pebble Bed breeding blanket for the European DEMO fusion
  reactor at the Frascati Neutron Generator," Nucl. Fusion 41 (2001) 1535.

  U. Fischer et al., "Neutronic analyses of the FNG HCPB mock-up
  experiment," Fusion Eng. Des. 69 (2003) 359-363.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNG HCPB Tritium Breeder Module Mock-up -- OpenMC fixed-source model"
)
parser.add_argument(
    "--particles",
    type=int,
    default=1_000_000,
    help="Number of source particles to simulate (default: 1,000,000)",
)
parser.add_argument(
    "--batches",
    type=int,
    default=10,
    help="Number of batches (default: 10)",
)
args = parser.parse_args()


# =============================================================================
# Materials
# =============================================================================
# In the HCPB blanket concept, FOUR classes of materials interact to produce
# tritium: the neutron multiplier (Be), the breeder (Li2CO3), the structural
# steel (SS), and the surrounding void. Each plays a distinct nuclear role.

# --- Beryllium (neutron multiplier) ---
# Beryllium is the KEY material in the HCPB concept. Its nuclear role is to
# multiply the neutron population via the Be-9(n,2n) reaction:
#   Be-9 + n -> 2 He-4 + 2n   (threshold 1.85 MeV, Q = -1.57 MeV)
# This reaction converts one fast neutron into two lower-energy neutrons,
# effectively doubling the number of neutrons available for tritium breeding.
# Without beryllium multiplication, TBR > 1 would be impossible with natural
# lithium and realistic blanket geometries.
#
# Additional nuclear effects in beryllium:
#   - Elastic scattering: Be-9 is a light nucleus (A=9), so elastic
#     scattering efficiently moderates (slows down) neutrons. This is
#     beneficial because the Li-6(n,t) cross section increases as 1/v
#     at low energies.
#   - (n,alpha) reaction: Be-9(n,alpha)He-6 competes with (n,2n) but has
#     a smaller cross section at 14 MeV.
#
# Density: 1.85 g/cm3 (standard for metallic beryllium blocks)
beryllium = openmc.Material(name="Beryllium (neutron multiplier)")
beryllium.add_nuclide("Be9", 1.0, "ao")  # 100% Be-9 (monoisotopic)
beryllium.set_density("g/cm3", 1.85)

# --- Li2CO3 powder (front breeder sections, main box) ---
# Lithium carbonate serves as a surrogate for the ceramic breeder pebbles
# (Li4SiO4 or Li2TiO3) used in the actual HCPB TBM design. The tritium is
# produced primarily by the Li-6(n,t)He-4 reaction:
#   Li-6 + n -> T + He-4 + 4.78 MeV   (exothermic, sigma_th ~ 940 b)
# with a secondary contribution from Li-7:
#   Li-7 + n -> T + He-4 + n' - 2.47 MeV   (endothermic, threshold ~2.5 MeV)
#
# The Li-6(n,t) reaction dominates tritium production because:
#   1. Its cross section follows a 1/v law and is enormous at thermal energies
#   2. Beryllium multiplication + moderation produces a large thermal flux
#   3. Even at natural enrichment (7.5% Li-6), the breeding is significant
#
# Composition: Li2CO3 -> 2 Li + 1 C + 3 O (per formula unit)
# Lithium isotopics: 7.5 at% Li-6, 92.5 at% Li-7 (natural abundance)
# Density: 1.123 g/cm3 (powder, not bulk crystal)
li2co3_front = openmc.Material(name="Li2CO3 powder (front breeder sections)")
li2co3_front.add_nuclide("Li6", 2.0 * 0.075, "ao")    # 7.5% of 2 Li atoms
li2co3_front.add_nuclide("Li7", 2.0 * 0.925, "ao")    # 92.5% of 2 Li atoms
li2co3_front.add_nuclide("C12", 0.9893, "ao")          # natural C (98.93% C-12)
li2co3_front.add_nuclide("C13", 0.0107, "ao")          # natural C (1.07% C-13)
li2co3_front.add_nuclide("O16", 3.0 * 0.99757, "ao")   # natural O (99.757% O-16)
li2co3_front.add_nuclide("O17", 3.0 * 0.00038, "ao")   # natural O (0.038% O-17)
li2co3_front.add_nuclide("O18", 3.0 * 0.00205, "ao")   # natural O (0.205% O-18)
li2co3_front.set_density("g/cm3", 1.123)

# --- Li2CO3 powder (rear cassette) ---
# Same chemical composition as front breeder, but lower packing density
# because the rear cassette was filled separately.
# Density: 0.9413 g/cm3 (looser packing in larger volume)
# Total mass measured: 11,690.4 +/- 0.1 g
li2co3_rear = openmc.Material(name="Li2CO3 powder (rear cassette)")
li2co3_rear.add_nuclide("Li6", 2.0 * 0.075, "ao")
li2co3_rear.add_nuclide("Li7", 2.0 * 0.925, "ao")
li2co3_rear.add_nuclide("C12", 0.9893, "ao")
li2co3_rear.add_nuclide("C13", 0.0107, "ao")
li2co3_rear.add_nuclide("O16", 3.0 * 0.99757, "ao")
li2co3_rear.add_nuclide("O17", 3.0 * 0.00038, "ao")
li2co3_rear.add_nuclide("O18", 3.0 * 0.00205, "ao")
li2co3_rear.set_density("g/cm3", 0.9413)

# --- AISI-303 stainless steel (main box walls, breeder cassette walls) ---
# Austenitic stainless steel used for the main containment box.
# AISI 303 is a free-machining variant of 304 with added sulfur.
# From a neutronics standpoint, the steel acts as a parasitic absorber --
# the Fe, Cr, and Ni all have non-negligible capture cross sections that
# compete with tritium breeding. Minimizing steel mass is important for TBR.
#
# Nominal composition (weight percent):
#   Fe: balance (~69.245%), Cr: 18%, Ni: 9%, Mn: 2%, Si: 1%,
#   C: 0.08%, S: 0.25%, P: 0.045%, N: 0.08%
# Density: 7.954 g/cm3
aisi303 = openmc.Material(name="AISI-303 stainless steel")
aisi303.add_element("Fe", 0.69245, "wo")  # balance
aisi303.add_element("Cr", 0.18, "wo")     # 18% chromium
aisi303.add_element("Ni", 0.09, "wo")     # 9% nickel
aisi303.add_element("Mn", 0.02, "wo")     # 2% manganese
aisi303.add_element("Si", 0.01, "wo")     # 1% silicon
aisi303.add_element("C", 0.0008, "wo")    # 0.08% carbon (max)
aisi303.add_element("S", 0.0025, "wo")    # 0.25% sulfur (free-machining)
aisi303.add_element("P", 0.00045, "wo")   # 0.045% phosphorus
aisi303.add_element("N", 0.0008, "wo")    # 0.08% nitrogen
aisi303.set_density("g/cm3", 7.954)

# --- AISI-316 stainless steel (rear cassette walls) ---
# AISI 316 is a molybdenum-bearing austenitic stainless steel with improved
# corrosion resistance. Used for the rear cassette structure.
# The Mo addition affects neutronics: Mo has significant resonance capture.
#
# Nominal composition (weight percent):
#   Fe: balance (~65.345%), Cr: 17%, Ni: 12%, Mo: 2.5%, Mn: 2%,
#   Si: 1%, C: 0.08%, P: 0.045%, S: 0.03%
# Density: 7.99 g/cm3
aisi316 = openmc.Material(name="AISI-316 stainless steel")
aisi316.add_element("Fe", 0.65345, "wo")  # balance
aisi316.add_element("Cr", 0.17, "wo")     # 17% chromium
aisi316.add_element("Ni", 0.12, "wo")     # 12% nickel
aisi316.add_element("Mo", 0.025, "wo")    # 2.5% molybdenum
aisi316.add_element("Mn", 0.02, "wo")     # 2% manganese
aisi316.add_element("Si", 0.01, "wo")     # 1% silicon
aisi316.add_element("C", 0.0008, "wo")    # 0.08% carbon (max)
aisi316.add_element("P", 0.00045, "wo")   # 0.045% phosphorus
aisi316.add_element("S", 0.0003, "wo")    # 0.03% sulfur
aisi316.set_density("g/cm3", 7.99)

# Collect all materials and export to XML
materials = openmc.Materials([
    beryllium, li2co3_front, li2co3_rear, aisi303, aisi316
])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# Coordinate system:
#   y-axis: beam axis (depth into assembly, neutrons travel in +y direction)
#   x-axis: horizontal transverse (width of assembly)
#   z-axis: vertical
#   Origin: centre of the front face of the main SS box
#
# The D-T source is at y = -5.3 cm (upstream of front face).
# The assembly extends from y = 0 to y = 43.8 cm (main box + rear cassette).
#
# Cross-sectional dimensions:
#   Main box: 31.0 cm (x) x 30.9 cm (z) -> x: -15.5 to +15.5, z: -15.45 to +15.45
#   Rear cassette: same x and z dimensions (31.0 x 30.9 cm)

# --- Key dimensions ---
# Main box (AISI-303)
main_box_x = 31.0     # cm, external width
main_box_y = 29.0     # cm, external depth along beam
main_box_z = 30.9     # cm, external height
wall_303 = 0.5        # cm, wall thickness of main box

# Rear cassette (AISI-316)
rear_box_x = 31.0     # cm, external width
rear_box_y = 14.8     # cm, external depth along beam
rear_box_z = 30.9     # cm, external height
wall_316 = 0.5        # cm, wall thickness of rear cassette

# Breeder layer parameters
breeder_thickness = 1.2    # cm, thickness of each Li2CO3 layer
divider_thickness = 0.1    # cm (1 mm), SS divider between paired layers
cassette_wall = 0.15       # cm, thin SS wall enclosing each breeder cassette

# Double-layer total thickness: wall + layer + divider + layer + wall
double_layer_total = (cassette_wall + breeder_thickness + divider_thickness
                      + breeder_thickness + cassette_wall)
# = 0.15 + 1.2 + 0.1 + 1.2 + 0.15 = 2.80 cm

# Source position
source_y = -5.3  # cm, D-T point source upstream of front face

# --- Layer positions along beam axis (y) ---
# Front wall of main box
y_main_front = 0.0                          # front face of main box
y_main_front_inner = y_main_front + wall_303  # 0.5 cm

# First double breeder layer section
# Positioned so that the measurement depth y=4.2 cm falls within it.
# The first breeder layer centre should be near y=4.2 cm.
# We place the first cassette starting at y=3.5 cm from the front face:
y_breeder1_start = 3.50                     # start of first breeder cassette wall
y_breeder1_wall1_end = y_breeder1_start + cassette_wall     # 3.65
y_breeder1a_end = y_breeder1_wall1_end + breeder_thickness  # 4.85
y_breeder1_div_end = y_breeder1a_end + divider_thickness    # 4.95
y_breeder1b_end = y_breeder1_div_end + breeder_thickness    # 6.15
y_breeder1_end = y_breeder1b_end + cassette_wall            # 6.30

# Second double breeder layer section
# Positioned so that measurement depth y=16.8 cm falls within it.
# We place the second cassette starting at y=15.90 cm:
y_breeder2_start = 15.90
y_breeder2_wall1_end = y_breeder2_start + cassette_wall     # 16.05
y_breeder2a_end = y_breeder2_wall1_end + breeder_thickness  # 17.25
y_breeder2_div_end = y_breeder2a_end + divider_thickness    # 17.35
y_breeder2b_end = y_breeder2_div_end + breeder_thickness    # 18.55
y_breeder2_end = y_breeder2b_end + cassette_wall            # 18.70

# Rear wall of main box
y_main_rear_inner = y_main_front + main_box_y - wall_303  # 28.5
y_main_rear = y_main_front + main_box_y                    # 29.0

# Rear cassette
y_rear_front = y_main_rear                                 # 29.0
y_rear_front_inner = y_rear_front + wall_316               # 29.5
y_rear_rear_inner = y_rear_front + rear_box_y - wall_316   # 43.3
y_rear_rear = y_rear_front + rear_box_y                    # 43.8

# --- Cross-sectional half-dimensions ---
half_x = main_box_x / 2.0   # 15.5 cm
half_z = main_box_z / 2.0   # 15.45 cm

# =============================================================================
# Surface definitions
# =============================================================================

# --- Y-planes (along beam axis) ---
# Main box boundaries
py_main_front = openmc.YPlane(y0=y_main_front, name="Main box front face (y=0)")
py_main_front_inner = openmc.YPlane(y0=y_main_front_inner,
                                     name="Main box front wall inner surface")
py_main_rear_inner = openmc.YPlane(y0=y_main_rear_inner,
                                    name="Main box rear wall inner surface")
py_main_rear = openmc.YPlane(y0=y_main_rear, name="Main box rear face")

# First double breeder layer boundaries
py_b1_start = openmc.YPlane(y0=y_breeder1_start,
                             name="Breeder 1 cassette front wall start")
py_b1_w1_end = openmc.YPlane(y0=y_breeder1_wall1_end,
                              name="Breeder 1a front (after front cassette wall)")
py_b1a_end = openmc.YPlane(y0=y_breeder1a_end,
                            name="Breeder layer 1a rear")
py_b1_div_end = openmc.YPlane(y0=y_breeder1_div_end,
                               name="Breeder 1 divider rear")
py_b1b_end = openmc.YPlane(y0=y_breeder1b_end,
                            name="Breeder layer 1b rear")
py_b1_end = openmc.YPlane(y0=y_breeder1_end,
                           name="Breeder 1 cassette rear wall end")

# Second double breeder layer boundaries
py_b2_start = openmc.YPlane(y0=y_breeder2_start,
                             name="Breeder 2 cassette front wall start")
py_b2_w1_end = openmc.YPlane(y0=y_breeder2_wall1_end,
                              name="Breeder 2a front (after front cassette wall)")
py_b2a_end = openmc.YPlane(y0=y_breeder2a_end,
                            name="Breeder layer 2a rear")
py_b2_div_end = openmc.YPlane(y0=y_breeder2_div_end,
                               name="Breeder 2 divider rear")
py_b2b_end = openmc.YPlane(y0=y_breeder2b_end,
                            name="Breeder layer 2b rear")
py_b2_end = openmc.YPlane(y0=y_breeder2_end,
                           name="Breeder 2 cassette rear wall end")

# Rear cassette boundaries
py_rear_front = openmc.YPlane(y0=y_rear_front,
                               name="Rear cassette front face")
py_rear_front_inner = openmc.YPlane(y0=y_rear_front_inner,
                                     name="Rear cassette front wall inner")
py_rear_rear_inner = openmc.YPlane(y0=y_rear_rear_inner,
                                    name="Rear cassette rear wall inner")
py_rear_rear = openmc.YPlane(y0=y_rear_rear,
                              name="Rear cassette rear face")

# --- X-planes (horizontal transverse) ---
px_min = openmc.XPlane(x0=-half_x, name="Assembly left face (x-)")
px_max = openmc.XPlane(x0=+half_x, name="Assembly right face (x+)")

# --- Z-planes (vertical) ---
pz_min = openmc.ZPlane(z0=-half_z, name="Assembly bottom face (z-)")
pz_max = openmc.ZPlane(z0=+half_z, name="Assembly top face (z+)")

# --- Vacuum boundary sphere ---
# Centred near the midpoint of the full assembly (main box + rear cassette).
# Total assembly length: 43.8 cm. Midpoint: ~22 cm.
# Sphere radius 60 cm is sufficient to enclose everything with margin.
boundary_sphere = openmc.Sphere(
    x0=0.0,
    y0=y_rear_rear / 2.0,  # ~21.9 cm, midpoint of full assembly
    z0=0.0,
    r=60.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Helper: cross-sectional region (within the box footprint) ---
# All cells in the assembly share the same x-z cross section.
def in_box_xz():
    """Return the x-z region common to all cells in the assembly."""
    return +px_min & -px_max & +pz_min & -pz_max


# =============================================================================
# Cell definitions
# =============================================================================
# Build cells from front to back along the beam axis.
# Each cell is commented with its physical role in the HCPB mock-up.

cells = []

# --- Main box front wall (AISI-303, 0.5 cm) ---
# The SS front wall is the first material the 14.1 MeV source neutrons
# encounter. It causes some scattering and a small amount of parasitic
# capture, slightly reducing the neutron flux entering the beryllium.
front_wall_cell = openmc.Cell(
    name="Main box front wall (AISI-303, 0.5 cm)",
    fill=aisi303,
    region=+py_main_front & -py_main_front_inner & in_box_xz(),
)
cells.append(front_wall_cell)

# --- Beryllium zone 1 (between front wall and first breeder section) ---
# This is the first beryllium region the neutrons pass through. Here, the
# 14.1 MeV neutrons begin to multiply via Be-9(n,2n) and slow down via
# elastic scattering. The thickness of this zone determines how many
# multiplied, moderated neutrons reach the first breeder layer.
be_zone1_cell = openmc.Cell(
    name="Beryllium zone 1 (front wall to first breeder)",
    fill=beryllium,
    region=+py_main_front_inner & -py_b1_start & in_box_xz(),
)
cells.append(be_zone1_cell)

# --- First breeder cassette front wall (AISI-303, 0.15 cm) ---
breeder1_front_wall_cell = openmc.Cell(
    name="Breeder cassette 1 front wall (AISI-303)",
    fill=aisi303,
    region=+py_b1_start & -py_b1_w1_end & in_box_xz(),
)
cells.append(breeder1_front_wall_cell)

# --- First breeder layer 1a (Li2CO3, 1.2 cm) ---
# This is the first tritium-producing layer. Neutrons arriving here have
# been partially moderated by the beryllium and some have been multiplied.
# The Li-6(n,t) reaction rate depends on the neutron spectrum -- lower
# energy neutrons are more efficiently captured.
breeder1a_cell = openmc.Cell(
    name="Breeder layer 1a (Li2CO3, 1.2 cm, first of double layer)",
    fill=li2co3_front,
    region=+py_b1_w1_end & -py_b1a_end & in_box_xz(),
)
cells.append(breeder1a_cell)

# --- First breeder divider wall (AISI-303, 1 mm) ---
# Thin SS wall separating the two layers in the double-layer configuration.
# This adds a small parasitic absorption but is structurally necessary.
breeder1_divider_cell = openmc.Cell(
    name="Breeder 1 divider wall (AISI-303, 1 mm)",
    fill=aisi303,
    region=+py_b1a_end & -py_b1_div_end & in_box_xz(),
)
cells.append(breeder1_divider_cell)

# --- First breeder layer 1b (Li2CO3, 1.2 cm) ---
# Second layer of the first double-layer section. Being deeper in the
# assembly, the neutron spectrum here is softer (more thermal) than in
# layer 1a, leading to a higher Li-6(n,t) reaction rate per unit flux.
breeder1b_cell = openmc.Cell(
    name="Breeder layer 1b (Li2CO3, 1.2 cm, second of double layer)",
    fill=li2co3_front,
    region=+py_b1_div_end & -py_b1b_end & in_box_xz(),
)
cells.append(breeder1b_cell)

# --- First breeder cassette rear wall (AISI-303, 0.15 cm) ---
breeder1_rear_wall_cell = openmc.Cell(
    name="Breeder cassette 1 rear wall (AISI-303)",
    fill=aisi303,
    region=+py_b1b_end & -py_b1_end & in_box_xz(),
)
cells.append(breeder1_rear_wall_cell)

# --- Beryllium zone 2 (between first and second breeder sections) ---
# This is the largest beryllium region, providing maximum neutron
# multiplication. Neutrons that pass through the first breeder section
# (without being captured by lithium) continue to multiply and moderate
# in this thick beryllium block before reaching the second breeder.
# The cumulative multiplication effect is crucial for achieving TBR > 1.
be_zone2_cell = openmc.Cell(
    name="Beryllium zone 2 (between breeder sections, main multiplier)",
    fill=beryllium,
    region=+py_b1_end & -py_b2_start & in_box_xz(),
)
cells.append(be_zone2_cell)

# --- Second breeder cassette front wall (AISI-303, 0.15 cm) ---
breeder2_front_wall_cell = openmc.Cell(
    name="Breeder cassette 2 front wall (AISI-303)",
    fill=aisi303,
    region=+py_b2_start & -py_b2_w1_end & in_box_xz(),
)
cells.append(breeder2_front_wall_cell)

# --- Second breeder layer 2a (Li2CO3, 1.2 cm) ---
# The second breeder section sees a different neutron spectrum than the
# first: lower flux but softer spectrum (more multiplied, moderated
# neutrons). The tritium production rate here provides crucial data on
# how Be multiplication compensates for flux attenuation with depth.
breeder2a_cell = openmc.Cell(
    name="Breeder layer 2a (Li2CO3, 1.2 cm, first of second double layer)",
    fill=li2co3_front,
    region=+py_b2_w1_end & -py_b2a_end & in_box_xz(),
)
cells.append(breeder2a_cell)

# --- Second breeder divider wall (AISI-303, 1 mm) ---
breeder2_divider_cell = openmc.Cell(
    name="Breeder 2 divider wall (AISI-303, 1 mm)",
    fill=aisi303,
    region=+py_b2a_end & -py_b2_div_end & in_box_xz(),
)
cells.append(breeder2_divider_cell)

# --- Second breeder layer 2b (Li2CO3, 1.2 cm) ---
breeder2b_cell = openmc.Cell(
    name="Breeder layer 2b (Li2CO3, 1.2 cm, second of second double layer)",
    fill=li2co3_front,
    region=+py_b2_div_end & -py_b2b_end & in_box_xz(),
)
cells.append(breeder2b_cell)

# --- Second breeder cassette rear wall (AISI-303, 0.15 cm) ---
breeder2_rear_wall_cell = openmc.Cell(
    name="Breeder cassette 2 rear wall (AISI-303)",
    fill=aisi303,
    region=+py_b2b_end & -py_b2_end & in_box_xz(),
)
cells.append(breeder2_rear_wall_cell)

# --- Beryllium zone 3 (after second breeder to rear of main box) ---
# Final beryllium region in the main box. This zone provides additional
# neutron reflection (backscatter) that can enhance breeding in the
# second breeder section. It also moderates neutrons that will enter
# the rear Li2CO3 cassette.
be_zone3_cell = openmc.Cell(
    name="Beryllium zone 3 (after second breeder to main box rear wall)",
    fill=beryllium,
    region=+py_b2_end & -py_main_rear_inner & in_box_xz(),
)
cells.append(be_zone3_cell)

# --- Main box rear wall (AISI-303, 0.5 cm) ---
rear_wall_cell = openmc.Cell(
    name="Main box rear wall (AISI-303, 0.5 cm)",
    fill=aisi303,
    region=+py_main_rear_inner & -py_main_rear & in_box_xz(),
)
cells.append(rear_wall_cell)

# --- Rear cassette front wall (AISI-316, 0.5 cm) ---
# The rear cassette uses AISI-316 instead of 303, which has a slightly
# different composition (Mo addition). This affects neutron transport
# through the steel slightly due to Mo resonances.
rear_cass_front_wall_cell = openmc.Cell(
    name="Rear cassette front wall (AISI-316, 0.5 cm)",
    fill=aisi316,
    region=+py_rear_front & -py_rear_front_inner & in_box_xz(),
)
cells.append(rear_cass_front_wall_cell)

# --- Rear cassette Li2CO3 fill ---
# The rear cassette is entirely filled with Li2CO3 powder at a lower
# density (0.9413 g/cm3) than the front breeder layers (1.123 g/cm3).
# This large volume of breeder material provides additional tritium
# production from deeply penetrating neutrons. The lower density results
# in a different self-shielding and moderation environment.
# Total mass: 11,690.4 g, measured to 0.1 g precision.
rear_cass_fill_cell = openmc.Cell(
    name="Rear cassette Li2CO3 fill (density 0.9413 g/cm3)",
    fill=li2co3_rear,
    region=+py_rear_front_inner & -py_rear_rear_inner & in_box_xz(),
)
cells.append(rear_cass_fill_cell)

# --- Rear cassette rear wall (AISI-316, 0.5 cm) ---
rear_cass_rear_wall_cell = openmc.Cell(
    name="Rear cassette rear wall (AISI-316, 0.5 cm)",
    fill=aisi316,
    region=+py_rear_rear_inner & -py_rear_rear & in_box_xz(),
)
cells.append(rear_cass_rear_wall_cell)

# --- Void region ---
# Everything inside the vacuum boundary sphere but outside the assembly.
# This includes the space between the source and the front face, and
# all surrounding void/air.
assembly_region = +py_main_front & -py_rear_rear & in_box_xz()
void_cell = openmc.Cell(
    name="Void (air/vacuum surrounding assembly)",
    region=-boundary_sphere & ~assembly_region,
)
cells.append(void_cell)

# --- Build the root universe and geometry ---
root_universe = openmc.Universe(name="Root universe", cells=cells)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The FNG D-T source is approximated as an isotropic point source of 14.1 MeV
# neutrons located 5.3 cm upstream of the assembly front face.
#
# In our coordinate system, the source is at (0, -5.3, 0).
#
# The real FNG source has an angular dependence of both intensity and mean
# energy due to D-T kinematics with 300 keV deuterons (the neutron energy
# varies from ~13.5 to ~14.8 MeV depending on angle). However, the isotropic
# 14.1 MeV approximation is the standard approach used in the SINBAD reference
# calculations and published analyses.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(0.0, source_y, 0.0))   # 5.3 cm from front face
source.angle = openmc.stats.Isotropic()                         # isotropic emission
source.energy = openmc.stats.Discrete([14.1e6], [1.0])          # 14.1 MeV mono-energetic
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"           # external source, not eigenvalue
settings.source = source                      # D-T point source defined above
settings.particles = args.particles           # particles per batch (from CLI)
settings.batches = args.batches               # number of batches (from CLI)
settings.photon_transport = False             # neutron-only (tritium production focus)
settings.output = {"tallies": True}          # write tally results to statepoint

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We define tallies to extract:
#   1. Tritium production rate in each breeder layer
#      - Li-6(n,t)He-4: the dominant breeding reaction (MT=205 = (n,Xt))
#      - We use the (n,Xt) score which tallies tritium production from
#        all reactions producing tritons (both Li-6(n,t) and Li-7(n,n't))
#   2. Neutron flux spectrum in each breeder layer
#   3. Cell flux at the 4 experimental measurement depths
#
# The (n,Xt) score in OpenMC tallies the total number of tritons produced
# per source particle in a given cell. This is the most direct measure of
# tritium breeding and corresponds to what pellet detectors measure in
# the experiment (total tritium produced in each Li2CO3 layer).

tallies = openmc.Tallies()

# --- Energy bin structure for spectral tallies ---
# Logarithmic bins from 1e-5 eV (thermal) to 15 MeV (above D-T peak).
# We extend to very low energies because the Li-6(n,t) cross section is
# largest at thermal energies, and beryllium moderation produces a
# significant thermal population.
energy_bins = np.logspace(
    np.log10(1.0e-5),   # 0.01 meV (deep thermal)
    np.log10(15.0e6),    # 15 MeV (above D-T source energy)
    201,                  # 201 edges -> 200 bins (~5 bins per decade)
)
energy_filter = openmc.EnergyFilter(energy_bins)

# --- Breeder layer cells for tritium production tallies ---
# These are the four Li2CO3 layers in the main box plus the rear cassette fill.
breeder_cells = [
    breeder1a_cell,     # Layer 1a (first double-layer, front)
    breeder1b_cell,     # Layer 1b (first double-layer, rear)
    breeder2a_cell,     # Layer 2a (second double-layer, front)
    breeder2b_cell,     # Layer 2b (second double-layer, rear)
    rear_cass_fill_cell,  # Rear cassette (large Li2CO3 volume)
]

breeder_names = [
    "breeder_1a", "breeder_1b",
    "breeder_2a", "breeder_2b",
    "rear_cassette",
]

# --- Tally: Tritium production rate in each breeder layer ---
# The (n,Xt) score counts the number of tritons produced per source neutron
# in each cell. This is the KEY result of this benchmark -- it directly
# measures the tritium breeding performance of the HCPB mock-up.
#
# Note on (n,Xt): This score tallies triton production from ALL reactions
# that produce at least one triton. In Li2CO3, the dominant contributors are:
#   - Li-6(n,t)He-4  (MT=105): dominant at thermal/epithermal energies
#   - Li-7(n,n't)He-4 (MT=112): secondary, requires E > 2.47 MeV
# The experiment measured total tritium production, so (n,Xt) is the
# appropriate score for comparison.
for i, (bcell, bname) in enumerate(zip(breeder_cells, breeder_names)):
    # Tritium production (energy-integrated)
    t_tally = openmc.Tally(name=f"tritium_production_{bname}")
    t_tally.filters = [openmc.CellFilter([bcell])]
    t_tally.scores = ["(n,Xt)"]  # total triton production
    tallies.append(t_tally)

    # Neutron flux spectrum in the breeder layer
    # This helps diagnose whether the spectrum is sufficiently thermalised
    # for efficient Li-6(n,t) breeding.
    spec_tally = openmc.Tally(name=f"flux_spectrum_{bname}")
    spec_tally.filters = [openmc.CellFilter([bcell]), energy_filter]
    spec_tally.scores = ["flux"]
    tallies.append(spec_tally)

# --- Tally: Neutron flux at experimental measurement depths ---
# The experiment had activation foils at 4 depths. We create thin tally
# cells at these depths to capture the neutron flux and activation rates.
# Measurement depths: y = 4.2, 10.5, 16.8, 23.1 cm from front surface.
measurement_depths = [4.2, 10.5, 16.8, 23.1]  # cm from front face

# --- Tally: Total tritium production (all breeder layers combined) ---
# This gives the overall TBR figure for the mock-up.
all_breeder_filter = openmc.CellFilter(breeder_cells)
total_tritium_tally = openmc.Tally(name="total_tritium_production")
total_tritium_tally.filters = [all_breeder_filter]
total_tritium_tally.scores = ["(n,Xt)"]
tallies.append(total_tritium_tally)

# --- Tally: Neutron flux in beryllium zones ---
# Useful for understanding the neutron multiplication and moderation.
# The Be-9(n,2n) reaction rate indicates the multiplication efficiency.
be_cells = [be_zone1_cell, be_zone2_cell, be_zone3_cell]
be_names = ["be_zone1", "be_zone2", "be_zone3"]

for bcell, bname in zip(be_cells, be_names):
    be_flux_tally = openmc.Tally(name=f"flux_{bname}")
    be_flux_tally.filters = [openmc.CellFilter([bcell])]
    be_flux_tally.scores = ["flux"]
    tallies.append(be_flux_tally)

    # Be-9(n,2n) reaction rate -- directly measures neutron multiplication
    be_n2n_tally = openmc.Tally(name=f"n2n_rate_{bname}")
    be_n2n_tally.filters = [openmc.CellFilter([bcell])]
    be_n2n_tally.scores = ["(n,2n)"]
    tallies.append(be_n2n_tally)

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print("=" * 72)
print("FNG HCPB Tritium Breeder Module Mock-up -- OpenMC Model")
print("=" * 72)
print(f"  Main box:       {main_box_x} x {main_box_y} x {main_box_z} cm "
      f"(AISI-303, {wall_303} cm walls)")
print(f"  Rear cassette:  {rear_box_x} x {rear_box_y} x {rear_box_z} cm "
      f"(AISI-316, {wall_316} cm walls)")
print(f"  Beryllium:      density {beryllium.density:.2f} g/cm3, "
      f"neutron multiplier via Be-9(n,2n)")
print(f"  Breeder (front): Li2CO3 powder, density 1.123 g/cm3, "
      f"4 layers x {breeder_thickness} cm")
print(f"  Breeder (rear):  Li2CO3 powder, density 0.9413 g/cm3, "
      f"fills rear cassette")
print(f"  Lithium:        Natural enrichment (7.5% Li-6, 92.5% Li-7)")
print(f"  Source:         14.1 MeV isotropic point at y = {source_y} cm")
print(f"  Particles:      {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins:    {len(energy_bins)-1} logarithmic bins, "
      f"{energy_bins[0]:.1e} - {energy_bins[-1]/1e6:.1f} MeV")
print(f"  Tallies:        Tritium production (5 layers + total), "
      f"flux spectra (5), Be(n,2n) rates (3)")
print()
print("  Key physics:")
print("    - Be-9(n,2n): neutron multiplication (threshold 1.85 MeV)")
print("    - Li-6(n,t):  tritium breeding (sigma_th ~ 940 b, 1/v law)")
print("    - Li-7(n,n't): secondary breeding (threshold 2.47 MeV)")
print("=" * 72)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
