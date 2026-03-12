#!/usr/bin/env python3
"""
FNS Dogleg Duct Streaming Experiment
======================================

Facility
--------
The Fusion Neutronics Source (FNS) facility is located at JAEA (formerly JAERI)
in Tokai-mura, Japan. FNS produces 14 MeV neutrons via the D-T reaction using
a 350 keV deuteron beam impinging on a rotating tritium target.

Experiment
----------
This benchmark models the "FNS Dogleg Duct Streaming Experiment" from the
SINBAD database. A large iron slab assembly (1700 x 1400 x 1800 mm) contains
a doubly-bent (dogleg) duct with a 300 x 300 mm square cross-section. The duct
consists of three legs connected at right angles:

  1st horizontal leg:  1150 mm along the beam direction (x-axis)
  Vertical leg:         600 mm perpendicular (z-axis)
  2nd horizontal leg:   650 mm parallel to the first but offset vertically

D-T neutrons enter the duct through its entrance on the front face of the iron
assembly. As neutrons travel through the dogleg duct, they must scatter off the
duct walls to navigate around the bends. This geometry tests the ability of
transport codes to model neutron streaming through voids with right-angle bends,
which is a challenging problem because:

  - Line-of-sight streaming dominates the first leg
  - Wall-scattering (albedo) dominates transport around the bends
  - Deep penetration through iron competes with streaming through the duct
  - The spectrum progressively softens as high-energy neutrons scatter at bends
  - Statistical convergence is difficult because few particles navigate both bends

These streaming and deep-penetration phenomena are directly relevant to fusion
reactor design, where ducts and penetrations (for heating systems, diagnostics,
and maintenance ports) must be routed through thick biological shields.

Geometry
--------
  - Iron assembly: 1700 mm (x) x 1400 mm (y) x 1800 mm (z)
  - Duct cross-section: 300 x 300 mm
  - 1st horizontal leg: extends 1150 mm along x-axis
  - Vertical leg: extends 600 mm along z-axis (connecting the two horizontal legs)
  - 2nd horizontal leg: extends 650 mm along x-axis (offset from 1st leg)
  - Right-angle connections between all legs
  - D-T point source aligned with the 1st horizontal duct leg entrance

Coordinate System
-----------------
  - x-axis: beam direction (horizontal, into the iron assembly)
  - y-axis: horizontal transverse
  - z-axis: vertical
  - Origin: centre of the assembly front face
  - The 1st duct leg is centred at y=0, z=0 on the front face
  - The vertical leg rises from z=0 to z=+60 cm
  - The 2nd duct leg is offset to z=+60 cm

Material
--------
  - Iron: natural Fe, density 7.874 g/cm3

Detectors
---------
  Four measurement positions within and beyond the duct:
    Position 1: Middle of 1st horizontal leg (x ~ 57.5 cm from entrance)
    Position 2: Junction of 1st leg and vertical leg
    Position 3: Junction of vertical leg and 2nd leg (top bend)
    Position 4: Middle of 2nd horizontal leg

  Measured quantities:
    - NE213 scintillator neutron spectra (> 2 MeV)
    - Activation dosimetry: Nb-93(n,2n), In-115(n,n'), Au-197(n,g)

Reference
---------
  SINBAD database: "FNS Dogleg Duct Streaming Experiment"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_duct/fnsstr-a.htm
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNS Dogleg Duct Streaming Experiment -- OpenMC fixed-source model"
)
parser.add_argument(
    "--particles",
    type=int,
    default=2_000_000,
    help="Number of source particles to simulate (default: 2,000,000)",
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

# --- Iron (natural) ---
# The FNS assembly was constructed from iron blocks. Natural iron is
# predominantly Fe-56 (91.7%) with minor contributions from Fe-54, Fe-57,
# and Fe-58. The density of pure iron is 7.874 g/cm3.
iron = openmc.Material(name="Iron (natural)")
iron.add_element("Fe", 1.0)          # 100% natural iron
iron.set_density("g/cm3", 7.874)      # standard density of pure iron

# Collect all materials and export to XML
materials = openmc.Materials([iron])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# The geometry consists of a large rectangular iron block with a dogleg duct
# carved through it using CSG boolean operations.
#
# The duct has three legs forming a Z-shape (or dogleg) when viewed from the
# side (XZ plane):
#
#   Side view (XZ plane, y=0):
#
#   z=75 cm  +---------+------------------------------------------+
#            |  IRON    |           IRON                            |
#   z=60 cm  |          +-----+                                     |
#            |    2nd leg duct |                                     |
#   z=30 cm  |          +-----+-----+                               |
#            |  IRON    | vert leg  |          IRON                 |
#   z=0  cm  +-----+   +-----+-----+                               |
#            | duct 1st leg   |                                     |
#   z=-15 cm +-----+   +-----+                                     |
#            |  IRON    |           IRON                            |
#   z=-75 cm +---------+------------------------------------------+
#            x=0       x=115  x=115+30                             x=170
#
# Note: dimensions in the diagram are approximate / schematic.

# --- Assembly dimensions (in cm) ---
# The iron assembly is 1700 x 1400 x 1800 mm = 170 x 140 x 180 cm.
assembly_x = 170.0   # cm, along beam direction
assembly_y = 140.0   # cm, horizontal transverse
assembly_z = 180.0   # cm, vertical

# --- Duct dimensions (in cm) ---
duct_width = 30.0     # cm (300 mm square cross-section)

# Duct leg lengths (along their respective axes):
leg1_length = 115.0   # cm (1150 mm) -- 1st horizontal leg along x
vert_length = 60.0    # cm (600 mm)  -- vertical leg along z
leg2_length = 65.0    # cm (650 mm)  -- 2nd horizontal leg along x

# --- Duct placement ---
# The 1st horizontal leg enters the front face of the assembly, centred at
# y=0, z=0. It extends from x=0 to x=leg1_length along the beam axis.
# The duct cross-section spans y=[-15, +15] and z=[-15, +15].
#
# The vertical leg connects the end of the 1st leg to the start of the 2nd
# leg. It rises from z=+15 (top of 1st leg) to z=+15+vert_length = z=75.
# In x, it occupies the same range as the junction: x in [leg1_length - duct_width, leg1_length + duct_width].
# Actually, the vertical leg shares one wall with the end of leg 1 and one
# wall with the start of leg 2. We place the vertical leg so it connects
# the two horizontal legs properly:
#   x-range: [leg1_length - duct_width, leg1_length]  (i.e., [85, 115])
#            Actually we need the vertical leg to span from the END of leg1
#            to the START of leg2. Both horizontal legs share the same x-range
#            for their connection to the vertical leg.

# Let's define the duct geometry precisely:
#
# 1st horizontal leg:
#   x: [0, leg1_length]   = [0, 115]
#   y: [-15, +15]
#   z: [-15, +15]
#
# Vertical leg (connects top of 1st leg to bottom of 2nd leg):
#   x: [leg1_length - duct_width, leg1_length]  = [85, 115]
#   y: [-15, +15]
#   z: [+15, +15 + vert_length] = [15, 75]
#   (This shares the top wall of leg1 and the bottom wall of leg2)
#
# 2nd horizontal leg:
#   x: [leg1_length - duct_width, leg1_length - duct_width + leg2_length] = [85, 150]
#   y: [-15, +15]
#   z: [+15 + vert_length - duct_width, +15 + vert_length] = [45, 75]
#   Wait -- the 2nd leg needs to be at the TOP of the vertical leg.
#   The vertical leg top is at z = 15 + 60 = 75.
#   So the 2nd leg spans z: [75 - 30, 75] = [45, 75].
#
# Let me reconsider. The vertical leg connects the two horizontal legs.
# The geometry should be:
#   - 1st leg centre at z=0, extends from z=-15 to z=+15
#   - Vertical leg goes UP from z=+15 to z=+15+60 = z=+75
#   - 2nd leg centre at z=+60 (offset by vert_length from 1st leg centre)
#     so it spans z=[60-15, 60+15] = [45, 75]
#
# The vertical leg in x should be at the junction, sharing walls:
#   x: [leg1_length - duct_width, leg1_length] = [85, 115]
#
# The 2nd horizontal leg extends AWAY from the vertical leg:
#   The 2nd leg could go either direction; physically it goes in the SAME
#   direction as the 1st leg (further into the assembly) but OFFSET.
#   Since the vertical leg is at x=[85,115], the 2nd leg starts at the
#   vertical leg and extends back toward the front face:
#   x: [leg1_length - duct_width - leg2_length + duct_width, leg1_length]
#   Hmm, let's think about this differently.
#
# The vertical leg is at the END of the 1st leg, so the 2nd leg must start
# from the vertical leg and go in some direction. Looking at the physical
# setup: the duct exits the rear of the assembly (or close to it). The 2nd
# leg extends from the top of the vertical leg BACK toward the entrance
# (or further into the assembly).
#
# For a classic dogleg: 1st leg goes forward, vertical leg goes up, 2nd leg
# goes forward (same direction). This means:
#   2nd leg x: [leg1_length - duct_width, leg1_length - duct_width + leg2_length]
#            = [85, 150]
#   Actually the 2nd leg should start where the vertical leg is and extend further.
#   The vertical leg occupies x=[85, 115]. The 2nd leg shares the vertical leg's
#   x-range and extends further:
#   x: [85, 85 + leg2_length] = [85, 150]
#
# Final duct coordinates (all in cm):

# Half-width of duct cross-section
hw = duct_width / 2.0  # 15 cm

# 1st horizontal leg: enters front face at (0, 0, 0), goes along +x
leg1_x_min = 0.0
leg1_x_max = leg1_length          # 115 cm
leg1_y_min = -hw                   # -15 cm
leg1_y_max = +hw                   # +15 cm
leg1_z_min = -hw                   # -15 cm
leg1_z_max = +hw                   # +15 cm

# Vertical leg: rises from top of leg1 to bottom of leg2
# In x, it aligns with the end of leg1 (sharing the last duct_width of leg1)
vert_x_min = leg1_length - duct_width  # 85 cm
vert_x_max = leg1_length               # 115 cm
vert_y_min = -hw                        # -15 cm
vert_y_max = +hw                        # +15 cm
vert_z_min = leg1_z_max                 # +15 cm (top of leg1)
vert_z_max = leg1_z_max + vert_length   # +75 cm

# 2nd horizontal leg: offset vertically, extends further into assembly
# Its z-range is at the top of the vertical leg
leg2_z_min = vert_z_max - duct_width   # 45 cm
leg2_z_max = vert_z_max                # 75 cm
leg2_y_min = -hw                        # -15 cm
leg2_y_max = +hw                        # +15 cm
# In x, it starts where the vertical leg is and extends further
leg2_x_min = vert_x_min                # 85 cm
leg2_x_max = vert_x_min + leg2_length  # 150 cm

print("Duct geometry (cm):")
print(f"  1st leg:  x=[{leg1_x_min}, {leg1_x_max}], "
      f"y=[{leg1_y_min}, {leg1_y_max}], z=[{leg1_z_min}, {leg1_z_max}]")
print(f"  Vert leg: x=[{vert_x_min}, {vert_x_max}], "
      f"y=[{vert_y_min}, {vert_y_max}], z=[{vert_z_min}, {vert_z_max}]")
print(f"  2nd leg:  x=[{leg2_x_min}, {leg2_x_max}], "
      f"y=[{leg2_y_min}, {leg2_y_max}], z=[{leg2_z_min}, {leg2_z_max}]")


# --- Define surfaces for the iron assembly ---
# The assembly is centred at y=0, and positioned so that the front face
# is at x=0. In z, the duct 1st leg is centred at z=0, so the assembly
# extends from z = -(assembly_z/2) to z = +(assembly_z/2). But we should
# verify the duct fits inside. The duct extends from z=-15 to z=+75, a
# total span of 90 cm. The assembly is 180 cm tall, so centring on the
# midpoint of the duct z-range: mid_z = (-15+75)/2 = 30. Let's just
# position the assembly so it fully contains the duct.
#
# We'll centre the assembly in y at y=0, and in z such that the duct is
# well within the iron. Place bottom of assembly at z = -75 and top at
# z = +105. This gives 180 cm total height and contains the duct
# (z = -15 to +75) with generous margins.

assembly_z_min = -75.0    # cm
assembly_z_max = assembly_z_min + assembly_z  # -75 + 180 = 105 cm
assembly_y_min = -assembly_y / 2.0  # -70 cm
assembly_y_max = +assembly_y / 2.0  # +70 cm
assembly_x_min = 0.0                 # front face at origin
assembly_x_max = assembly_x          # 170 cm

# Iron assembly bounding surfaces
sx_min = openmc.XPlane(x0=assembly_x_min, name="Assembly front face (x-min)")
sx_max = openmc.XPlane(x0=assembly_x_max, name="Assembly rear face (x-max)")
sy_min = openmc.YPlane(y0=assembly_y_min, name="Assembly left face (y-min)")
sy_max = openmc.YPlane(y0=assembly_y_max, name="Assembly right face (y-max)")
sz_min = openmc.ZPlane(z0=assembly_z_min, name="Assembly bottom face (z-min)")
sz_max = openmc.ZPlane(z0=assembly_z_max, name="Assembly top face (z-max)")

# The full assembly region (before duct is carved out)
assembly_region = +sx_min & -sx_max & +sy_min & -sy_max & +sz_min & -sz_max


# --- Define surfaces for the duct legs ---
# Each duct leg is a rectangular parallelepiped (box) defined by 6 planes.
# We'll create the duct as the UNION of 3 box regions, then subtract it
# from the iron assembly to create the void channel.

# 1st horizontal leg surfaces
leg1_xmin_s = openmc.XPlane(x0=leg1_x_min, name="Leg1 x-min")
leg1_xmax_s = openmc.XPlane(x0=leg1_x_max, name="Leg1 x-max")
leg1_ymin_s = openmc.YPlane(y0=leg1_y_min, name="Leg1 y-min")
leg1_ymax_s = openmc.YPlane(y0=leg1_y_max, name="Leg1 y-max")
leg1_zmin_s = openmc.ZPlane(z0=leg1_z_min, name="Leg1 z-min")
leg1_zmax_s = openmc.ZPlane(z0=leg1_z_max, name="Leg1 z-max")

leg1_region = (
    +leg1_xmin_s & -leg1_xmax_s &
    +leg1_ymin_s & -leg1_ymax_s &
    +leg1_zmin_s & -leg1_zmax_s
)

# Vertical leg surfaces
vert_xmin_s = openmc.XPlane(x0=vert_x_min, name="Vert leg x-min")
vert_xmax_s = openmc.XPlane(x0=vert_x_max, name="Vert leg x-max")
vert_ymin_s = openmc.YPlane(y0=vert_y_min, name="Vert leg y-min")
vert_ymax_s = openmc.YPlane(y0=vert_y_max, name="Vert leg y-max")
vert_zmin_s = openmc.ZPlane(z0=vert_z_min, name="Vert leg z-min")
vert_zmax_s = openmc.ZPlane(z0=vert_z_max, name="Vert leg z-max")

vert_region = (
    +vert_xmin_s & -vert_xmax_s &
    +vert_ymin_s & -vert_ymax_s &
    +vert_zmin_s & -vert_zmax_s
)

# 2nd horizontal leg surfaces
leg2_xmin_s = openmc.XPlane(x0=leg2_x_min, name="Leg2 x-min")
leg2_xmax_s = openmc.XPlane(x0=leg2_x_max, name="Leg2 x-max")
leg2_ymin_s = openmc.YPlane(y0=leg2_y_min, name="Leg2 y-min")
leg2_ymax_s = openmc.YPlane(y0=leg2_y_max, name="Leg2 y-max")
leg2_zmin_s = openmc.ZPlane(z0=leg2_z_min, name="Leg2 z-min")
leg2_zmax_s = openmc.ZPlane(z0=leg2_z_max, name="Leg2 z-max")

leg2_region = (
    +leg2_xmin_s & -leg2_xmax_s &
    +leg2_ymin_s & -leg2_ymax_s &
    +leg2_zmin_s & -leg2_zmax_s
)

# --- Combined duct region (union of all three legs) ---
# The union of the three rectangular boxes forms the complete dogleg duct.
# Note: where legs overlap (at the junctions), the union naturally handles
# the shared volume -- no double-counting occurs because it's a single void.
duct_region = leg1_region | vert_region | leg2_region


# --- Detector cells ---
# We place thin detector volumes (5 cm thick, matching duct cross-section)
# at 4 positions within the duct to tally neutron flux.
#
# Detector positions (chosen to capture streaming physics at key locations):
#   Position 1: Middle of 1st horizontal leg
#               x = leg1_length / 2 = 57.5 cm, z = 0
#   Position 2: Bottom bend (junction of 1st leg and vertical leg)
#               x = leg1_length - duct_width/2 = 100 cm, z = duct_width/2 = 15 cm
#               (inside the vertical leg, just above the 1st leg)
#   Position 3: Top bend (junction of vertical leg and 2nd leg)
#               x = vert_x_min + duct_width/2 = 100 cm, z = vert_z_max - duct_width/2 = 60 cm
#               (inside the vertical leg, just below the 2nd leg)
#   Position 4: Middle of 2nd horizontal leg
#               x = leg2_x_min + leg2_length/2 = 117.5 cm, z = (leg2_z_min+leg2_z_max)/2 = 60 cm

# Detector thickness along the leg direction
det_thickness = 5.0  # cm

# Each detector is a thin slab within the duct, oriented perpendicular to
# the local flow direction.

detector_cells = []
detector_regions = []
detector_descriptions = []

# --- Detector 1: Middle of 1st horizontal leg ---
# Thin slab perpendicular to x-axis at x = 57.5 cm
det1_x_centre = leg1_length / 2.0  # 57.5 cm
det1_xmin = openmc.XPlane(x0=det1_x_centre - det_thickness / 2.0,
                           name="Det1 x-min")
det1_xmax = openmc.XPlane(x0=det1_x_centre + det_thickness / 2.0,
                           name="Det1 x-max")
det1_region = +det1_xmin & -det1_xmax & +leg1_ymin_s & -leg1_ymax_s & +leg1_zmin_s & -leg1_zmax_s
detector_regions.append(det1_region)
detector_descriptions.append(
    f"Detector 1: middle of 1st leg, x={det1_x_centre} cm, z=0"
)

# --- Detector 2: Bottom bend (just inside vertical leg above 1st leg) ---
# Thin slab perpendicular to z-axis at z = 15 + 2.5 = 17.5 cm
det2_z_centre = vert_z_min + det_thickness / 2.0  # 17.5 cm
det2_zmin = openmc.ZPlane(z0=det2_z_centre - det_thickness / 2.0,
                           name="Det2 z-min")
det2_zmax = openmc.ZPlane(z0=det2_z_centre + det_thickness / 2.0,
                           name="Det2 z-max")
det2_region = +vert_xmin_s & -vert_xmax_s & +vert_ymin_s & -vert_ymax_s & +det2_zmin & -det2_zmax
detector_regions.append(det2_region)
detector_descriptions.append(
    f"Detector 2: bottom bend, x=[{vert_x_min},{vert_x_max}] cm, z={det2_z_centre} cm"
)

# --- Detector 3: Top bend (just inside vertical leg below 2nd leg) ---
# Thin slab perpendicular to z-axis at z = 75 - 2.5 = 72.5 cm
det3_z_centre = vert_z_max - det_thickness / 2.0  # 72.5 cm
det3_zmin = openmc.ZPlane(z0=det3_z_centre - det_thickness / 2.0,
                           name="Det3 z-min")
det3_zmax = openmc.ZPlane(z0=det3_z_centre + det_thickness / 2.0,
                           name="Det3 z-max")
det3_region = +vert_xmin_s & -vert_xmax_s & +vert_ymin_s & -vert_ymax_s & +det3_zmin & -det3_zmax
detector_regions.append(det3_region)
detector_descriptions.append(
    f"Detector 3: top bend, x=[{vert_x_min},{vert_x_max}] cm, z={det3_z_centre} cm"
)

# --- Detector 4: Middle of 2nd horizontal leg ---
# Thin slab perpendicular to x-axis at x = 85 + 65/2 = 117.5 cm
det4_x_centre = leg2_x_min + leg2_length / 2.0  # 117.5 cm
det4_xmin = openmc.XPlane(x0=det4_x_centre - det_thickness / 2.0,
                           name="Det4 x-min")
det4_xmax = openmc.XPlane(x0=det4_x_centre + det_thickness / 2.0,
                           name="Det4 x-max")
det4_region = +det4_xmin & -det4_xmax & +leg2_ymin_s & -leg2_ymax_s & +leg2_zmin_s & -leg2_zmax_s
detector_regions.append(det4_region)
detector_descriptions.append(
    f"Detector 4: middle of 2nd leg, x={det4_x_centre} cm, z={(leg2_z_min+leg2_z_max)/2} cm"
)

# Create detector cells (void -- the duct is empty)
for i, (det_reg, det_desc) in enumerate(zip(detector_regions, detector_descriptions)):
    det_cell = openmc.Cell(
        name=det_desc,
        fill=None,       # void (air) -- the duct is an empty channel
        region=det_reg,
    )
    detector_cells.append(det_cell)


# --- Duct void cell (excluding detector volumes) ---
# The duct itself is void (air). We subtract the detector volumes from the
# duct region so there are no overlapping cells.
duct_void_region = duct_region
for det_reg in detector_regions:
    duct_void_region = duct_void_region & ~det_reg

duct_cell = openmc.Cell(
    name="Duct void (air channel through iron)",
    fill=None,         # void -- neutrons stream freely through the duct
    region=duct_void_region,
)


# --- Bulk iron cell (assembly minus duct) ---
# The iron assembly is the rectangular block with the duct carved out.
iron_cell = openmc.Cell(
    name="Iron assembly (bulk shielding)",
    fill=iron,
    region=assembly_region & ~duct_region,
)


# --- Vacuum boundary ---
# A large sphere enclosing the entire geometry. Everything outside the
# assembly but inside this sphere is void. The boundary is vacuum (particles
# that reach it are killed).
boundary_sphere = openmc.Sphere(
    x0=assembly_x / 2.0,   # centre at assembly midpoint
    y0=0.0,
    z0=(assembly_z_min + assembly_z_max) / 2.0,
    r=200.0,                 # large enough to enclose everything + source
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# Void outside the assembly but inside the boundary sphere
exterior_void = openmc.Cell(
    name="Exterior void (outside assembly)",
    region=-boundary_sphere & ~assembly_region,
)


# --- Build geometry ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[iron_cell, duct_cell, *detector_cells, exterior_void],
)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The FNS D-T source is modelled as an isotropic point source of 14.1 MeV
# neutrons. The source is aligned with the 1st horizontal duct leg entrance,
# positioned just outside the front face of the iron assembly.
#
# Source position: x = -10 cm (100 mm in front of the assembly), y = 0, z = 0
# This places the source on the beam axis, aligned with the centre of the
# 1st duct leg entrance.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(-10.0, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV D-T neutrons
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = source
settings.particles = args.particles
settings.batches = args.batches
settings.photon_transport = False          # neutron-only transport
settings.output = {"tallies": True}

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We tally neutron flux in each of the 4 detector cells, both as energy
# spectra and as energy-integrated totals. This allows us to:
#   - Compare spectral shapes at different duct positions
#   - Quantify the streaming attenuation through the dogleg bends
#   - Observe spectral softening as neutrons scatter around the bends

# --- Energy bin structure ---
# Logarithmic bins from 10 keV to 15 MeV (100 bins).
# This covers the NE213 measurement range (> 2 MeV) and extends to lower
# energies where scattered neutrons accumulate after wall reflections.
energy_bins = np.logspace(
    np.log10(1.0e4),    # 10 keV in eV
    np.log10(15.0e6),   # 15 MeV in eV
    101,                  # 101 edges -> 100 bins
)

energy_filter = openmc.EnergyFilter(energy_bins)

tallies = openmc.Tallies()

# --- Spectral flux tallies at each detector position ---
for i, det_cell in enumerate(detector_cells):
    cell_filter = openmc.CellFilter([det_cell])
    tally = openmc.Tally(name=f"flux_detector_{i+1}")
    tally.filters = [cell_filter, energy_filter]
    tally.scores = ["flux"]
    tallies.append(tally)

# --- Energy-integrated flux tallies at each detector position ---
# These give a single number for quick attenuation comparisons.
for i, det_cell in enumerate(detector_cells):
    cell_filter = openmc.CellFilter([det_cell])
    tally = openmc.Tally(name=f"total_flux_detector_{i+1}")
    tally.filters = [cell_filter]
    tally.scores = ["flux"]
    tallies.append(tally)

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print()
print("=" * 70)
print("FNS Dogleg Duct Streaming Experiment -- OpenMC Model")
print("=" * 70)
print(f"  Assembly:    {assembly_x} x {assembly_y} x {assembly_z} cm iron block")
print(f"  Duct:        {duct_width} x {duct_width} cm cross-section, 3-leg dogleg")
print(f"    1st leg:   {leg1_length} cm along x (beam direction)")
print(f"    Vert leg:  {vert_length} cm along z (perpendicular)")
print(f"    2nd leg:   {leg2_length} cm along x (offset vertically)")
print(f"  Material:    Natural iron, rho = 7.874 g/cm3")
print(f"  Source:      14.1 MeV isotropic point at (-10, 0, 0) cm")
print(f"  Detectors:   {len(detector_cells)} positions in duct")
for desc in detector_descriptions:
    print(f"    - {desc}")
print(f"  Particles:   {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins: {len(energy_bins)-1} logarithmic bins, "
      f"{energy_bins[0]/1e6:.4f} - {energy_bins[-1]/1e6:.1f} MeV")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
