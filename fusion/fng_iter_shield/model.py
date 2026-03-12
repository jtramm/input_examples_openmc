#!/usr/bin/env python3
"""
FNG-ITER Bulk Shield Mock-up Benchmark: Inboard Shield Simulation
===================================================================

Facility
--------
The Frascati Neutron Generator (FNG) is located at ENEA Frascati, near Rome,
Italy. FNG produces 14.1 MeV neutrons via the D-T (deuterium-tritium) fusion
reaction using a 300 keV deuteron beam impinging on a tritium target. The
facility was purpose-built for fusion neutronics experiments, with a maximum
yield of approximately 1e11 neutrons per second. It has been a principal
European facility for integral validation of nuclear cross-section data
relevant to ITER and future fusion power plants.

Experiment
----------
This benchmark models the "FNG-ITER Bulk Shield Mock-up" from the SINBAD
database. A ~94 cm thick assembly simulating the ITER inboard shielding was
irradiated by 14.1 MeV D-T neutrons. The mock-up reproduces the layered
structure of the ITER inboard shielding, comprising (from front to back):
  - First wall (copper)
  - Blanket region (stainless steel SS316 and Perspex/PMMA sandwich)
  - Vacuum vessel (stainless steel SS316 and Perspex/PMMA sandwich)
  - Toroidal field (TF) coil simulation (alternating copper and SS316 plates)

Activation foils were placed at 14 depth positions to measure neutron reaction
rate profiles through the mock-up. Four reaction types were measured:
  - Nb-93(n,2n)Nb-92m    (threshold ~9 MeV, sensitive to 14 MeV peak)
  - Al-27(n,alpha)Na-24   (threshold ~3.3 MeV, fast neutron indicator)
  - Ni-58(n,p)Co-58       (threshold ~1 MeV, intermediate-fast neutrons)
  - Au-197(n,gamma)Au-198 (thermal/epithermal capture, low-energy sensitivity)

TLD-300 (CaF2:Tm) thermoluminescent dosimeters were placed at 17 positions
to measure nuclear heating profiles.

ITER Shielding Physics Background
----------------------------------
The ITER inboard shield must attenuate the 14 MeV fusion neutron flux by
roughly 10 orders of magnitude to protect the superconducting toroidal field
coils. The design relies on a multi-material layered approach:

  * Copper (first wall): Excellent thermal conductivity removes surface
    heat loads. Cu has moderate neutron scattering cross sections and
    contributes some moderation via inelastic scattering. The (n,2n)
    reaction in Cu is significant above ~10 MeV, effectively multiplying
    the neutron population while degrading the mean energy.

  * Stainless steel SS316: The primary structural and shielding material.
    Iron-based alloys are effective neutron attenuators due to:
      - High inelastic scattering cross section (Fe-56) above ~0.85 MeV
      - Deep minima in the elastic scattering cross section near 24 keV
        ("iron window") that require careful treatment in shielding codes
      - Significant (n,gamma) capture cross section at thermal energies
      - Nickel and chromium provide additional capture and scattering
    The Mn, Mo, and Si minor constituents also contribute capture.

  * Perspex (PMMA, C5O2H8): Acts as a water-equivalent moderator. The
    hydrogen in Perspex is critical for thermalizing fast neutrons via
    elastic scattering (the most efficient moderator nuclide). The
    moderated neutrons are then captured by the surrounding steel.
    This moderation-capture cycle is the essential mechanism for
    deep-penetration shielding effectiveness.

  * Alternating Cu/SS section: Simulates the TF coil winding pack, which
    in ITER consists of Nb3Sn superconductor embedded in a copper
    stabilizer matrix, enclosed in a stainless steel conduit (the
    "cable-in-conduit conductor" or CICC design).

Deep-penetration shielding calculations are numerically challenging because
the neutron flux attenuates by many orders of magnitude. Statistical
uncertainties in Monte Carlo codes grow exponentially with shield thickness.
Variance reduction techniques (weight windows, source biasing) are typically
needed for production calculations but are not applied in this basic model.

Layer Configuration
-------------------
The mock-up consists of the following layers (from front face to rear),
based on the published experimental configuration from Batistoni et al.:

  Layer  Material    Thickness (cm)  Cumulative (cm)  ITER Component
  -----  --------    --------------  ---------------  --------------
  1      Cu          0.53            0.53             First wall
  2      SS316       5.36            5.89             First wall back plate
  3      Perspex     5.70            11.59            Blanket breeder zone
  4      SS316       5.56            17.15            Blanket back wall
  5      Perspex     6.80            23.95            Blanket coolant/gap
  6      SS316       6.85            30.80            Vacuum vessel inner shell
  7      Perspex     4.25            35.05            VV shielding fill
  8      SS316       6.80            41.85            VV mid-plate
  9      Perspex     5.00            46.85            VV shielding fill
  10     SS316       6.95            53.80            VV outer shell
  11     Cu          6.75            60.55            TF coil case
  12     SS316       6.85            67.40            TF coil structure
  13     Cu          7.00            74.40            TF coil conductor
  14     SS316       6.70            81.10            TF coil structure
  15     Cu          6.65            87.75            TF coil conductor
  16     SS316       4.40            92.15            TF coil structure
  17     Cu          2.15            94.30            TF coil rear

Source-to-mockup distance: 5.3 cm
Total mock-up thickness: ~94.3 cm

Materials
---------
  - Copper: density 8.96 g/cm^3, natural isotopic composition (69.17% Cu-63,
    30.83% Cu-65). OFHC (Oxygen-Free High-Conductivity) grade.

  - SS316 stainless steel: density 7.93 g/cm^3, composition by weight:
    Fe ~65.4%, Cr 17.0%, Ni 12.0%, Mo 2.5%, Mn 2.0%, Si 1.0%, C 0.08%
    Standard austenitic stainless steel used extensively in nuclear
    applications for its corrosion resistance and weldability.

  - Perspex (PMMA, polymethyl methacrylate): density 1.19 g/cm^3,
    chemical formula (C5O2H8)n, used as a water-equivalent moderator.
    Perspex provides hydrogen for neutron moderation without the
    containment issues of liquid water.

Detector Positions
------------------
  14 activation foil positions at depths from front face (cm):
    3.43, 10.32, 17.15, 23.95, 30.80, 41.85, 46.85, 53.80,
    60.55, 67.40, 74.40, 81.10, 87.75, 92.15

  These positions correspond to interfaces between major layers,
  allowing measurement of the neutron field after each shielding
  component.

Reference
---------
  SINBAD database, NEA-1553: "FNG/ITER Bulk Shield"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/

  P. Batistoni et al., "Neutronics Shield Experiment for ITER at the
  Frascati Neutron Generator FNG," Fusion Eng. Des., vol. 47, pp. 25-60,
  1999.

  M. Angelone, P. Batistoni, M. Pillon, "Neutron Streaming Experiment
  and Analysis of the FNG-ITER Bulk Shield Mock-up," Fusion Eng. Des.,
  2001.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNG-ITER Bulk Shield Mock-up -- OpenMC fixed-source model"
)
parser.add_argument(
    "--particles",
    type=int,
    default=2_000_000,
    # Default is 2M particles. Deep penetration through ~94 cm of mixed
    # shielding attenuates the neutron flux by many orders of magnitude.
    # The rear detectors will have very poor statistics without variance
    # reduction or very high particle counts. For production runs, consider
    # 10M-100M particles or use weight windows.
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

# --- Copper (OFHC, Oxygen-Free High-Conductivity) ---
# Copper is used in the first wall and TF coil simulation sections.
# In ITER, CuCrZr and CuAl25 alloys serve as heat sinks bonded to
# beryllium armour tiles. For this benchmark, pure copper is used.
# Natural isotopic composition: 69.17% Cu-63, 30.83% Cu-65
# Density: 8.96 g/cm^3 (standard for OFHC copper at room temperature)
copper = openmc.Material(name="Copper")
copper.add_nuclide("Cu63", 0.6917, "ao")   # 69.17 atom% Cu-63
copper.add_nuclide("Cu65", 0.3083, "ao")   # 30.83 atom% Cu-65
copper.set_density("g/cm3", 8.96)

# --- Stainless steel SS316 ---
# SS316 is the workhorse structural material for ITER's vacuum vessel,
# blanket support structures, and port plugs. Its composition includes:
#   - Iron (balance, ~65.4 wt%): primary shielding element; Fe-56 has
#     strong inelastic scattering above 0.85 MeV and resonance capture
#   - Chromium (~17 wt%): corrosion resistance; Cr-52 contributes
#     inelastic scattering and has capture resonances
#   - Nickel (~12 wt%): austenite stabilizer; Ni-58 has important
#     (n,p) and (n,alpha) reactions above ~1 MeV
#   - Molybdenum (~2.5 wt%): strengthening at high temperature;
#     Mo isotopes have large capture cross sections
#   - Manganese (~2 wt%): deoxidizer; Mn-55 has significant capture
#   - Silicon (~1 wt%): deoxidizer; Si-28 has small cross sections
#   - Carbon (~0.08 wt%): strengthening; kept low to avoid sensitization
# Density: 7.93 g/cm^3 (typical for SS316 at room temperature)
ss316 = openmc.Material(name="SS316")
ss316.add_element("Fe", 0.6542, "wo")  # iron, balance (~65.42 wt%)
ss316.add_element("Cr", 0.1700, "wo")  # chromium, 17.0 wt%
ss316.add_element("Ni", 0.1200, "wo")  # nickel, 12.0 wt%
ss316.add_element("Mo", 0.0250, "wo")  # molybdenum, 2.5 wt%
ss316.add_element("Mn", 0.0200, "wo")  # manganese, 2.0 wt%
ss316.add_element("Si", 0.0100, "wo")  # silicon, 1.0 wt%
ss316.add_element("C",  0.0008, "wo")  # carbon, 0.08 wt%
ss316.set_density("g/cm3", 7.93)

# --- Perspex (PMMA, polymethyl methacrylate) ---
# Chemical formula: (C5O2H8)n
# Perspex serves as a water-equivalent moderator in the mock-up. The
# hydrogen content is the key neutronics feature: H-1 is the most
# efficient moderator nuclide because its mass is nearly equal to the
# neutron mass, allowing maximum energy transfer per collision.
# A 14 MeV neutron requires only ~25 elastic collisions with hydrogen
# to thermalize (compared to ~110 for carbon, ~500 for iron).
# Once thermalized, neutrons are captured by surrounding steel.
# This moderation-capture cycle is the fundamental mechanism of the
# ITER inboard shield.
#
# Molecular formula C5O2H8 gives the following atom fractions:
#   H: 8/15 = 0.5333
#   C: 5/15 = 0.3333
#   O: 2/15 = 0.1333
# Density: 1.19 g/cm^3
perspex = openmc.Material(name="Perspex (PMMA)")
perspex.add_element("H", 8, "ao")   # 8 hydrogen atoms per monomer unit
perspex.add_element("C", 5, "ao")   # 5 carbon atoms per monomer unit
perspex.add_element("O", 2, "ao")   # 2 oxygen atoms per monomer unit
perspex.set_density("g/cm3", 1.19)

# Collect all materials into a Materials object and export
materials = openmc.Materials([copper, ss316, perspex])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# Coordinate system:
#   - x-axis: beam axis (source propagation direction, neutrons travel in +x)
#   - y-axis: horizontal transverse axis
#   - z-axis: vertical transverse axis
#   - Origin: centre of the front face of the mock-up assembly
#   - Source position: x = -5.3 cm (53 mm upstream of front face)
#   - Assembly front face: x = 0 cm
#   - Assembly rear face:  x ~ 94.3 cm
#
# The mock-up is modelled as a one-dimensional slab geometry along the beam
# axis. The actual experimental assembly was approximately 60 x 60 cm in
# cross-section, but for this benchmark the lateral dimensions are large
# enough relative to the beam spot that a 1D treatment is standard.

# --- Assembly cross-section dimensions ---
assembly_width = 60.0     # cm (y-extent: -30 to +30 cm)
assembly_height = 60.0    # cm (z-extent: -30 to +30 cm)

# --- Source position ---
source_x = -5.3  # cm (53 mm upstream of front face)

# =============================================================================
# Layer configuration
# =============================================================================
# The mock-up consists of 17 layers of three materials, reproducing the
# layered structure of the ITER inboard shield. The layer sequence is:
#
# FIRST WALL section (layers 1-2):
#   The ITER first wall faces the plasma and must withstand intense
#   neutron and heat loads. It consists of beryllium armour on a copper
#   alloy heat sink backed by a stainless steel support structure.
#
# BLANKET section (layers 3-6):
#   The breeding blanket surrounds the plasma and serves three functions:
#   (1) breed tritium fuel via Li-6(n,t)He-4, (2) convert neutron kinetic
#   energy to heat for power generation, (3) provide neutron shielding.
#   In ITER's test blanket modules, the breeder zone contains lithium
#   ceramics or Pb-Li eutectic. Perspex simulates the hydrogenous
#   breeder/coolant material.
#
# VACUUM VESSEL section (layers 7-10):
#   The ITER vacuum vessel is a massive double-walled structure of SS316
#   filled with water and steel shielding blocks. Its primary shielding
#   function is to attenuate neutron flux to protect the TF coils.
#   The SS316/Perspex sandwich in the mock-up simulates the steel/water
#   combination. This section dominates the overall shielding performance.
#
# TOROIDAL FIELD COIL section (layers 11-17):
#   The ITER TF coils use Nb3Sn superconductor (critical temperature ~18 K)
#   in a "cable-in-conduit conductor" (CICC) configuration: superconducting
#   cable cooled by forced-flow supercritical helium, enclosed in a
#   stainless steel conduit, with copper stabilizer. The alternating
#   Cu/SS316 plates simulate this composite structure.
#   Nuclear heating in the TF coils must be kept below ~10-20 kW to
#   maintain cryogenic temperatures. This is the primary design driver
#   for the inboard shield thickness.

# Each layer is defined as (material_object, thickness_cm, description)
layers = [
    # ---- FIRST WALL ----
    # Layer 1: Thin copper plate simulating the first-wall heat sink.
    # In ITER, this is a CuCrZr alloy bonded to beryllium plasma-facing tiles.
    (copper, 0.53, "First wall Cu heat sink"),

    # Layer 2: SS316 plate simulating the first-wall structural back-plate.
    # This stainless steel plate provides structural support for the first wall
    # and begins the neutron attenuation.
    (ss316,  5.36, "First wall SS316 back plate"),

    # ---- BLANKET ----
    # Layer 3: Perspex simulating the blanket breeder zone.
    # In ITER, this region contains lithium ceramic pebbles (Li4SiO4 or Li2TiO3)
    # for tritium breeding, or Pb-Li eutectic in liquid metal blanket concepts.
    # Perspex approximates the hydrogenous moderating properties.
    (perspex, 5.70, "Blanket breeder zone (Perspex)"),

    # Layer 4: SS316 plate simulating the blanket back wall / manifold plate.
    # This structural plate separates the breeder zone from the coolant manifold.
    (ss316,  5.56, "Blanket back wall SS316"),

    # Layer 5: Perspex simulating the blanket-to-VV gap / coolant space.
    # In ITER, this region contains water coolant piping and inter-component gaps.
    (perspex, 6.80, "Blanket-VV gap (Perspex)"),

    # ---- VACUUM VESSEL ----
    # Layer 6: SS316 plate - vacuum vessel inner shell.
    # The ITER VV inner shell is ~60 mm thick SS316, forming the first boundary
    # of the double-walled vessel. This is a primary structural element.
    (ss316,  6.85, "VV inner shell SS316"),

    # Layer 7: Perspex - VV shielding fill (simulates water + steel filling).
    # The space between the VV inner and outer shells is filled with water
    # and SS316 shielding blocks for neutron attenuation.
    (perspex, 4.25, "VV shielding fill (Perspex)"),

    # Layer 8: SS316 plate - VV mid-plate / internal stiffener.
    # Internal ribs and plates provide structural stiffness to the VV
    # while contributing to shielding.
    (ss316,  6.80, "VV mid-plate SS316"),

    # Layer 9: Perspex - VV shielding fill (second moderator layer).
    # Additional water/moderator layer for thermalizing remaining fast neutrons.
    (perspex, 5.00, "VV shielding fill (Perspex)"),

    # Layer 10: SS316 plate - vacuum vessel outer shell.
    # The VV outer shell completes the double-walled vessel structure.
    (ss316,  6.95, "VV outer shell SS316"),

    # ---- TOROIDAL FIELD COIL SIMULATION ----
    # Layers 11-17: Alternating Cu and SS316 plates simulating the TF coil
    # winding pack and case structure. In ITER, the TF coil consists of:
    #   - Nb3Sn superconductor (cable-in-conduit conductor)
    #   - Copper stabilizer (absorbs energy during quench events)
    #   - SS316 conduit and structural case
    #   - Insulation (glass-epoxy) between turns

    # Layer 11: Copper - TF coil case / conductor stabilizer
    (copper, 6.75, "TF coil Cu case"),

    # Layer 12: SS316 - TF coil structural element
    (ss316,  6.85, "TF coil SS316 structure"),

    # Layer 13: Copper - TF coil conductor stabilizer
    (copper, 7.00, "TF coil Cu conductor"),

    # Layer 14: SS316 - TF coil structural element
    (ss316,  6.70, "TF coil SS316 structure"),

    # Layer 15: Copper - TF coil conductor stabilizer
    (copper, 6.65, "TF coil Cu conductor"),

    # Layer 16: SS316 - TF coil rear structural plate
    (ss316,  4.40, "TF coil SS316 rear plate"),

    # Layer 17: Copper - TF coil rear plate
    (copper, 2.15, "TF coil Cu rear plate"),
]

# Compute cumulative layer boundaries (front face at x=0)
layer_boundaries = [0.0]  # front face at x = 0
for mat, thickness, desc in layers:
    layer_boundaries.append(layer_boundaries[-1] + thickness)

# Total mock-up thickness
total_thickness = layer_boundaries[-1]

# --- Bounding surfaces ---
# Lateral boundaries of the assembly (60 x 60 cm cross-section)
y_min = openmc.YPlane(y0=-assembly_width / 2.0, name="Assembly left face (y-)")
y_max = openmc.YPlane(y0=+assembly_width / 2.0, name="Assembly right face (y+)")
z_min = openmc.ZPlane(z0=-assembly_height / 2.0, name="Assembly bottom face (z-)")
z_max = openmc.ZPlane(z0=+assembly_height / 2.0, name="Assembly top face (z+)")

# --- Vacuum boundary sphere ---
# The source is at x = -5.3, assembly extends to x ~ 94.3.
# A sphere of radius 60 cm centred near the assembly midpoint encloses
# the full geometry with margin. We use a slightly larger radius to be safe.
boundary_sphere = openmc.Sphere(
    x0=total_thickness / 2.0,   # centre at midpoint of assembly along x
    y0=0.0,
    z0=0.0,
    r=65.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Layer boundary x-planes ---
# Create an x-plane at each layer interface, including front and rear faces.
layer_planes = []
for i, x_pos in enumerate(layer_boundaries):
    plane = openmc.XPlane(
        x0=x_pos,
        name=f"Layer boundary at x = {x_pos:.2f} cm",
    )
    layer_planes.append(plane)

# --- Detector cells (activation foils at 14 positions) ---
# Activation foils are placed at specific depths within the mock-up.
# Each foil is a thin disc at the interface between layers.
# The positions correspond to boundaries between major shield components,
# allowing measurement of the progressive attenuation through each section.
#
# Positions are measured from the front face of the assembly (cm):
detector_depths_cm = [
    3.43,   # Position 1:  inside first wall SS316 layer
    10.32,  # Position 2:  inside blanket Perspex layer
    17.15,  # Position 3:  at blanket back wall / Perspex interface
    23.95,  # Position 4:  at blanket-VV gap / SS316 interface
    30.80,  # Position 5:  at VV inner shell / Perspex interface
    41.85,  # Position 6:  at VV mid-plate / Perspex interface
    46.85,  # Position 7:  at VV Perspex / SS316 interface
    53.80,  # Position 8:  at VV outer shell / Cu interface
    60.55,  # Position 9:  at TF coil Cu / SS316 interface
    67.40,  # Position 10: at TF coil SS316 / Cu interface
    74.40,  # Position 11: at TF coil Cu / SS316 interface
    81.10,  # Position 12: at TF coil SS316 / Cu interface
    87.75,  # Position 13: at TF coil Cu / SS316 interface
    92.15,  # Position 14: at TF coil SS316 / Cu interface (near rear)
]

# Detector foils are modelled as thin slabs (2 mm thick, 2 cm x 2 cm)
# centred on the beam axis. They are filled with the surrounding material
# since the actual foils are negligibly thin and do not perturb transport.
detector_half_thickness = 0.1   # cm (total thickness 2 mm)
detector_half_width = 1.0       # cm (2 cm x 2 cm foil area)

# Determine which material each detector sits in, based on its depth
def get_material_at_depth(depth):
    """Return the material object for the layer containing the given depth.

    Each layer spans from layer_boundaries[i] to layer_boundaries[i+1].
    The detector at a given depth is assigned the material of the
    enclosing layer.
    """
    for i, (mat, thickness, desc) in enumerate(layers):
        if depth <= layer_boundaries[i + 1] + 1e-10:
            return mat
    # If beyond last layer, return last layer's material
    return layers[-1][0]

detector_cells = []       # will hold the Cell objects for tallying
detector_regions = []     # will hold the regions to subtract from bulk layers

for i, depth in enumerate(detector_depths_cm):
    # Each detector foil is bounded by two x-planes and y/z limits.
    # The foil extends from (depth - 0.1) to (depth + 0.1) cm in x,
    # and +/- 1 cm in y and z.
    det_x_lo = openmc.XPlane(
        x0=depth - detector_half_thickness,
        name=f"Detector {i+1} front (depth {depth:.2f} cm)",
    )
    det_x_hi = openmc.XPlane(
        x0=depth + detector_half_thickness,
        name=f"Detector {i+1} back (depth {depth:.2f} cm)",
    )
    det_y_lo = openmc.YPlane(y0=-detector_half_width, name=f"Det {i+1} y-")
    det_y_hi = openmc.YPlane(y0=+detector_half_width, name=f"Det {i+1} y+")
    det_z_lo = openmc.ZPlane(z0=-detector_half_width, name=f"Det {i+1} z-")
    det_z_hi = openmc.ZPlane(z0=+detector_half_width, name=f"Det {i+1} z+")

    # Region: inside the small rectangular foil volume
    det_region = (
        +det_x_lo & -det_x_hi
        & +det_y_lo & -det_y_hi
        & +det_z_lo & -det_z_hi
    )

    # Fill with surrounding layer material (foil is negligibly thin)
    det_material = get_material_at_depth(depth)
    det_cell = openmc.Cell(
        name=f"Detector {i+1} at depth {depth:.2f} cm",
        fill=det_material,
        region=det_region,
    )
    detector_cells.append(det_cell)
    detector_regions.append(det_region)

# --- Bulk layer cells ---
# Create individual cells for each of the 17 layers.
# Each layer extends from layer_planes[i] to layer_planes[i+1] in x,
# and across the full 60 x 60 cm cross-section in y and z.
# Detector regions that overlap with each layer are subtracted to avoid
# geometry overlaps.
layer_cells = []
for i, (mat, thickness, desc) in enumerate(layers):
    # Base region for this layer
    layer_region = (
        +layer_planes[i] & -layer_planes[i + 1]
        & +y_min & -y_max
        & +z_min & -z_max
    )

    # Subtract overlapping detector regions.
    # A detector at depth d extends from d-0.1 to d+0.1 cm.
    # It overlaps this layer if the intervals intersect.
    layer_x_lo = layer_boundaries[i]
    layer_x_hi = layer_boundaries[i + 1]
    for j, det_region in enumerate(detector_regions):
        det_depth = detector_depths_cm[j]
        det_x_lo = det_depth - detector_half_thickness
        det_x_hi = det_depth + detector_half_thickness
        if det_x_lo < layer_x_hi and det_x_hi > layer_x_lo:
            layer_region = layer_region & ~det_region

    layer_cell = openmc.Cell(
        name=f"Layer {i+1}: {desc} (x={layer_boundaries[i]:.2f}-{layer_boundaries[i+1]:.2f} cm)",
        fill=mat,
        region=layer_region,
    )
    layer_cells.append(layer_cell)

# --- Void region (air/vacuum surrounding the assembly) ---
# Everything inside the boundary sphere but outside the assembly.
front_face = layer_planes[0]
rear_face = layer_planes[-1]
assembly_box_region = +front_face & -rear_face & +y_min & -y_max & +z_min & -z_max

# Void: inside boundary sphere but outside the assembly
void_region = -boundary_sphere & ~assembly_box_region

# Subtract any detector regions that might extend outside the assembly
# (possible edge effects near front/rear faces)
for det_region in detector_regions:
    void_region = void_region & ~det_region

void_cell = openmc.Cell(
    name="Void (air/vacuum surrounding assembly)",
    region=void_region,
)

# --- Build the universe and geometry ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[*layer_cells, *detector_cells, void_cell],
)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The FNG D-T source is approximated as an isotropic point source of 14.1 MeV
# neutrons located 53 mm (5.3 cm) upstream of the assembly front face.
# In our coordinate system, the source is at x = -5.3 cm, y = 0, z = 0.
#
# The real FNG source has angular and energy dependence from the D-T reaction
# kinematics with 300 keV deuterons (the neutron energy varies from ~13.5 MeV
# at 180 degrees to ~14.8 MeV at 0 degrees in the lab frame). However, the
# isotropic 14.1 MeV approximation is standard for this benchmark as
# documented in the SINBAD database.
#
# Source intensity is measured by the associated-particle method using a
# silicon surface barrier detector to count the alpha particles from the
# D + T -> n + He-4 reaction. This provides an absolute calibration of the
# neutron yield with approximately +/- 2% accuracy.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(source_x, 0.0, 0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV mono-energetic
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"        # external D-T neutron source
settings.source = source                  # point source defined above
settings.particles = args.particles       # particles per batch (default 2M)
settings.batches = args.batches           # number of batches (default 10)
settings.photon_transport = False         # neutron-only transport for this run

# Deep penetration note:
# This is a challenging shielding problem. The neutron flux attenuates by
# roughly 8-10 orders of magnitude across the ~94 cm mock-up. Without
# variance reduction, the statistical uncertainty at rear detectors will
# be very large. For production calculations, weight windows should be
# generated (e.g., via the MAGIC method or ADVANTG) and applied.
# The default 2M particles gives reasonable statistics for the first
# ~60 cm but poor statistics at the rear.

settings.output = {"tallies": True}
settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We define tallies to extract:
#   1. Neutron flux spectrum at each of the 14 detector positions
#   2. Total (energy-integrated) flux at each detector position
#
# The spectral tallies allow comparison of the energy-dependent neutron
# field at each shield depth. Key spectral features to look for:
#   - 14.1 MeV source peak (attenuated with depth)
#   - Inelastic scattering continuum (1-10 MeV)
#   - Iron resonance structure (especially the 24 keV "window")
#   - Thermal/epithermal peak (from hydrogen moderation in Perspex)
#
# The activation foil reactions measured in the experiment have different
# energy sensitivities, probing different parts of the spectrum:
#   - Nb-93(n,2n): threshold ~9 MeV -- tracks the 14 MeV peak
#   - Al-27(n,alpha): threshold ~3.3 MeV -- fast neutrons
#   - Ni-58(n,p): threshold ~1 MeV -- intermediate-to-fast
#   - Au-197(n,gamma): no threshold -- thermal/epithermal capture

# --- Energy bin structure ---
# 100 logarithmic bins from 1e-3 eV (cold neutrons) to 15 MeV.
# Wide range needed because Perspex moderator produces thermal neutrons.
energy_bins = np.logspace(
    np.log10(1.0e-3),   # 1 meV (thermal neutrons)
    np.log10(15.0e6),   # 15 MeV (above D-T source energy)
    201,                 # 201 edges -> 200 bins (~5 bins per decade)
)

energy_filter = openmc.EnergyFilter(energy_bins)
tallies = openmc.Tallies()

# --- Pre-create cell filters to avoid duplicate warnings ---
cell_filters = []
for det_cell in detector_cells:
    cell_filters.append(openmc.CellFilter([det_cell]))

# --- Tallies 1-14: Cell flux spectra at each detector position ---
# Track-length estimator for scalar neutron flux in each thin detector cell.
# Units: [neutrons-cm / source-particle]
for i, cf in enumerate(cell_filters):
    tally = openmc.Tally(name=f"flux_detector_{i+1}")
    tally.filters = [cf, energy_filter]
    tally.scores = ["flux"]
    tallies.append(tally)

# --- Tallies 15-28: Total flux at each detector position ---
# Energy-integrated flux for quick attenuation profile.
for i, cf in enumerate(cell_filters):
    tally = openmc.Tally(name=f"total_flux_detector_{i+1}")
    tally.filters = [cf]
    tally.scores = ["flux"]
    tallies.append(tally)

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print("=" * 72)
print("FNG-ITER Bulk Shield Mock-up -- OpenMC Model")
print("=" * 72)
print(f"  Mock-up:     {len(layers)} layers, {total_thickness:.2f} cm total thickness")
print(f"  Cross-sect:  {assembly_width} x {assembly_height} cm")
print(f"  Materials:   Copper ({copper.density:.2f} g/cc), "
      f"SS316 ({ss316.density:.2f} g/cc), "
      f"Perspex ({perspex.density:.2f} g/cc)")
print(f"  Source:      14.1 MeV isotropic point at x = {source_x} cm")
print(f"  Detectors:   {len(detector_cells)} positions at depths:")
for i, d in enumerate(detector_depths_cm):
    mat_name = get_material_at_depth(d).name
    print(f"    {i+1:2d}. {d:6.2f} cm  ({mat_name})")
print(f"  Particles:   {args.particles:,} per batch x {args.batches} batches "
      f"= {args.particles * args.batches:,} total")
print(f"  Energy bins: {len(energy_bins)-1} logarithmic bins, "
      f"{energy_bins[0]:.1e} - {energy_bins[-1]:.1e} eV")
print(f"  Tallies:     Flux spectra ({len(detector_cells)}), "
      f"total flux ({len(detector_cells)})")
print()
print("Layer structure:")
print(f"  {'#':>3s}  {'Material':12s}  {'Thick (cm)':>10s}  {'Cumul (cm)':>10s}  Description")
print(f"  {'---':>3s}  {'--------':12s}  {'----------':>10s}  {'----------':>10s}  -----------")
for i, (mat, thick, desc) in enumerate(layers):
    print(f"  {i+1:3d}  {mat.name:12s}  {thick:10.2f}  {layer_boundaries[i+1]:10.2f}  {desc}")
print("=" * 72)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
