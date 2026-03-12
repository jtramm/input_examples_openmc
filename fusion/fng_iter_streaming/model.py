#!/usr/bin/env python3
"""
FNG-ITER Streaming Experiment: Neutron Streaming Through an ITER-Like Duct
===========================================================================

Facility
--------
The Frascati Neutron Generator (FNG) is located at ENEA Frascati, near Rome,
Italy. FNG produces 14.1 MeV neutrons via the D-T (deuterium-tritium) fusion
reaction using a 300 keV deuteron beam impinging on a tritium target. The
facility was purpose-built for fusion neutronics experiments, with a maximum
yield of approximately 1e11 neutrons per second.

Experiment
----------
This benchmark models the "FNG-ITER Streaming Experiment" from the SINBAD
database. It is the most complex fusion shielding benchmark performed at FNG,
designed to test computational methods for predicting neutron streaming through
narrow channels in ITER-like shielding assemblies.

The assembly is a 100 x 100 cm block, 94.26 cm deep, composed of alternating
layers of perspex (PMMA, a water-equivalent moderator), stainless steel
AISI-316, and copper. A narrow cylindrical channel (28 mm inner diameter,
39.07 cm long, 1 mm SS316 wall) penetrates the front portion of the assembly,
simulating a diagnostic or instrumentation duct through the ITER blanket/shield.
At the end of the channel, a rectangular detector cavity (52 x 148 x 48 mm)
houses activation foils. Behind the cavity, a "magnetic coil simulator" of
alternating copper and SS316 plates represents the superconducting magnet
region. The rear ~30 cm is polyethylene shielding.

Two source configurations were measured:
  - On-axis:  D-T source centred on the channel axis (streaming geometry)
  - Off-axis: D-T source shifted 5.3 cm laterally (oblique streaming)

Streaming Physics
-----------------
This experiment is fundamentally a STREAMING problem. In fusion reactor
shielding, narrow penetrations (diagnostic ports, cooling pipes, heating
ducts) can allow neutrons to travel deep into the shield without
significant interaction, bypassing tens of centimetres of shielding material.
The 28 mm channel here has an aspect ratio of ~14:1 (length:diameter),
meaning neutrons travelling within a narrow forward cone can traverse the
entire channel in a single flight.

Key physics:
  - Direct streaming: 14 MeV neutrons travel the full channel length
    without scattering, producing the dominant contribution at channel exit
  - Wall scattering: neutrons interact with the 1 mm SS316 channel wall,
    producing a secondary scattered component
  - Albedo: neutrons scatter back into the channel from surrounding material
  - Deep penetration: beyond the cavity, flux attenuation through bulk
    shielding material follows near-exponential behaviour

Monte Carlo simulation of streaming problems is challenging because:
  - The channel subtends a tiny solid angle from the source (~2.2e-5 sr)
  - Most source particles miss the channel entirely
  - Variance reduction (weight windows) is typically essential for
    acceptable statistics at depth, but this model uses analog transport

Geometry
--------
The assembly is oriented along the z-axis (beam direction):
  - z-axis: beam/channel axis (source at -z, assembly extends in +z)
  - x-axis: horizontal transverse
  - y-axis: vertical transverse
  - Origin: centre of the front face of the assembly
  - Source: z = -5.3 cm (on-axis) or offset by 5.3 cm in x (off-axis)

Layer structure (from front to back, all thicknesses in cm):
  z =   0.00 -   1.00 : Copper front plate (first wall simulator, 1 cm)
  z =   1.00 -   7.00 : Perspex (PMMA) block
  z =   7.00 -   7.30 : SS316 plate (3 mm)
  z =   7.30 -  13.30 : Perspex block
  z =  13.30 -  13.60 : SS316 plate (3 mm)
  z =  13.60 -  19.60 : Perspex block
  z =  19.60 -  19.90 : SS316 plate (3 mm)
  z =  19.90 -  25.90 : Perspex block
  z =  25.90 -  26.20 : SS316 plate (3 mm)
  z =  26.20 -  32.20 : Perspex block
  z =  32.20 -  32.50 : SS316 plate (3 mm)
  z =  32.50 -  38.50 : Perspex block
  z =  38.50 -  39.07 : SS316 plate (5.7 mm, channel end plate)
  z =  39.07 -  39.33 : Air gap (2.6 mm)
  z =  39.33 -  44.53 : Detector cavity region (52 mm in z, SS316 box)
  z =  44.53 -  46.35 : SS316 plate behind cavity (18.2 mm)
  z =  46.35 -  63.15 : Cu/SS316 coil simulator (~16.8 cm)
  z =  63.15 -  94.26 : Polyethylene shield (~31 cm)

Channel: 28 mm ID, 1 mm SS316 wall, extends from z = 0 to z = 39.07 cm

Detector cavity: inner dimensions 148 mm (x) x 48 mm (y) x 52 mm (z)
  centred on channel axis, from z = 39.33 to z = 44.53 cm

Materials
---------
  - SS316:  Fe/Cr/Ni/Mo/Mn stainless steel, 7.93 g/cm3
  - Copper: OFHC copper, 8.96 g/cm3
  - Perspex (PMMA): C5O2H8, 1.19 g/cm3
  - Polyethylene: CH2, 0.95 g/cm3

References
----------
  SINBAD database, NEA-1553: "FNG/ITER Streaming Experiment"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/fng_str/fngstr-a.htm

  P. Batistoni et al., "Neutron Streaming Experiment at FNG for ITER
  Shielding," Fusion Eng. Des., 1998.

  R. Villari et al., "Analysis of the FNG Streaming Experiment,"
  Fusion Eng. Des., 2000.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNG-ITER Streaming Experiment -- OpenMC fixed-source model"
)
parser.add_argument(
    "--particles",
    type=int,
    default=5_000_000,
    help="Number of source particles to simulate (default: 5,000,000)",
)
parser.add_argument(
    "--batches",
    type=int,
    default=10,
    help="Number of batches (default: 10)",
)
parser.add_argument(
    "--config",
    choices=["on-axis", "off-axis"],
    default="on-axis",
    help="Source configuration: on-axis (aligned with channel) or off-axis "
         "(shifted 5.3 cm laterally). Default: on-axis.",
)
args = parser.parse_args()


# =============================================================================
# Materials
# =============================================================================

# --- Stainless Steel AISI-316 ---
# Standard austenitic stainless steel used extensively in fusion reactor
# structural components: blanket modules, vacuum vessel, shield blocks.
# Composition (weight percent, typical for 316):
#   Fe: 65.395%, Cr: 17.0%, Ni: 12.0%, Mo: 2.5%, Mn: 2.0%,
#   Si: 1.0%, C: 0.08%, P: 0.025%
# Density: 7.93 g/cm3
ss316 = openmc.Material(name="SS316")
ss316.add_element("Fe", 0.65395, "wo")  # balance iron
ss316.add_element("Cr", 0.17, "wo")     # 17% chromium
ss316.add_element("Ni", 0.12, "wo")     # 12% nickel
ss316.add_element("Mo", 0.025, "wo")    # 2.5% molybdenum
ss316.add_element("Mn", 0.02, "wo")     # 2% manganese
ss316.add_element("Si", 0.01, "wo")     # 1% silicon
ss316.add_element("C", 0.0008, "wo")    # 0.08% carbon
ss316.add_element("P", 0.00025, "wo")   # 0.025% phosphorus
ss316.set_density("g/cm3", 7.93)

# --- OFHC Copper ---
# Oxygen-Free High-Conductivity copper, used for first-wall heat sinks
# and magnet stabiliser in ITER. >= 99.95% purity.
# Natural isotopic composition: 69.17% Cu-63, 30.83% Cu-65
# Density: 8.96 g/cm3
copper = openmc.Material(name="Copper")
copper.add_nuclide("Cu63", 0.6917, "ao")  # 69.17 atom% Cu-63
copper.add_nuclide("Cu65", 0.3083, "ao")  # 30.83 atom% Cu-65
copper.set_density("g/cm3", 8.96)

# --- Perspex (PMMA, polymethyl methacrylate) ---
# Chemical formula: (C5O2H8)n -- used as a water-equivalent moderator
# in the FNG assembly. It provides hydrogen-rich moderation similar to
# water without the containment issues.
# Density: 1.19 g/cm3
perspex = openmc.Material(name="Perspex (PMMA)")
perspex.add_element("H", 8, "ao")   # 8 hydrogen atoms per monomer unit
perspex.add_element("C", 5, "ao")   # 5 carbon atoms per monomer unit
perspex.add_element("O", 2, "ao")   # 2 oxygen atoms per monomer unit
perspex.set_density("g/cm3", 1.19)

# --- Polyethylene ---
# Chemical formula: (CH2)n -- used as the rear neutron shield.
# High hydrogen content makes it an excellent neutron moderator and
# thermaliser. Used behind the coil simulator to reduce backscatter.
# Density: 0.95 g/cm3
polyethylene = openmc.Material(name="Polyethylene")
polyethylene.add_element("H", 2, "ao")  # 2 hydrogen per CH2 unit
polyethylene.add_element("C", 1, "ao")  # 1 carbon per CH2 unit
polyethylene.set_density("g/cm3", 0.95)

# --- Air (for gap regions) ---
# Dry air approximation for the small gap between channel end and cavity.
air = openmc.Material(name="Air")
air.add_element("N", 0.784, "ao")   # 78.4% nitrogen
air.add_element("O", 0.216, "ao")   # 21.6% oxygen
air.set_density("g/cm3", 0.001205)  # standard air density

# Collect all materials
materials = openmc.Materials([ss316, copper, perspex, polyethylene, air])
materials.export_to_xml()


# =============================================================================
# Geometry -- layer structure along z-axis
# =============================================================================
# Coordinate system:
#   z-axis: beam axis (neutrons travel in +z direction)
#   x-axis: horizontal transverse
#   y-axis: vertical transverse
#   Origin: centre of front face of assembly
#   Assembly front face: z = 0, rear face: z = 94.26 cm
#   Source: z = -5.3 cm (on-axis)

# --- Assembly outer dimensions ---
assembly_half_width = 50.0   # cm (100 x 100 cm cross-section)
assembly_depth = 94.26       # cm total thickness

# --- Define the layer stack ---
# Each entry: (z_start, z_end, material, description)
# This layer structure is derived from the SINBAD specification.
# The front section consists of perspex blocks separated by thin SS316 plates,
# with a 1 cm copper first-wall simulator at the very front.
layers = [
    # Front copper plate (first wall simulator)
    (0.00, 1.00, copper, "Cu front plate (first wall)"),
    # Alternating perspex/SS316 layers (blanket/shield simulator)
    (1.00, 7.00, perspex, "Perspex block 1"),
    (7.00, 7.30, ss316, "SS316 plate 1 (3 mm)"),
    (7.30, 13.30, perspex, "Perspex block 2"),
    (13.30, 13.60, ss316, "SS316 plate 2 (3 mm)"),
    (13.60, 19.60, perspex, "Perspex block 3"),
    (19.60, 19.90, ss316, "SS316 plate 3 (3 mm)"),
    (19.90, 25.90, perspex, "Perspex block 4"),
    (25.90, 26.20, ss316, "SS316 plate 4 (3 mm)"),
    (26.20, 32.20, perspex, "Perspex block 5"),
    (32.20, 32.50, ss316, "SS316 plate 5 (3 mm)"),
    (32.50, 38.50, perspex, "Perspex block 6"),
    (38.50, 39.07, ss316, "SS316 channel end plate (5.7 mm)"),
    # Gap between channel exit and detector cavity
    (39.07, 39.33, air, "Air gap before cavity"),
    # Detector cavity region -- modelled as air-filled box inside SS316
    # (cavity walls are SS316; the inner cavity volume is air)
    (39.33, 44.53, ss316, "SS316 cavity block (contains detector cavity)"),
    # SS316 plate behind cavity
    (44.53, 46.35, ss316, "SS316 plate behind cavity"),
    # Cu/SS316 coil simulator -- alternating Cu and SS316 plates
    # Simplified as alternating 2 cm layers over ~16.8 cm
    (46.35, 48.35, copper, "Cu coil plate 1"),
    (48.35, 50.35, ss316, "SS316 coil plate 1"),
    (50.35, 52.35, copper, "Cu coil plate 2"),
    (52.35, 54.35, ss316, "SS316 coil plate 2"),
    (54.35, 56.35, copper, "Cu coil plate 3"),
    (56.35, 58.35, ss316, "SS316 coil plate 3"),
    (58.35, 60.35, copper, "Cu coil plate 4"),
    (60.35, 63.15, ss316, "SS316 coil plate 4"),
    # Polyethylene rear shield
    (63.15, 94.26, polyethylene, "Polyethylene rear shield"),
]

# --- Streaming channel parameters ---
channel_inner_radius = 1.4     # cm (28 mm diameter / 2)
channel_outer_radius = 1.5     # cm (inner + 1 mm SS316 wall)
channel_z_start = 0.0          # cm (starts at front face)
channel_z_end = 39.07          # cm (ends at channel end plate)

# --- Detector cavity parameters ---
# Inner dimensions: 148 mm (x) x 48 mm (y) x 52 mm (z)
# Centred on channel axis (x=0, y=0)
cavity_half_x = 7.4            # cm (148 mm / 2)
cavity_half_y = 2.4            # cm (48 mm / 2)
cavity_z_start = 39.33         # cm
cavity_z_end = 44.53           # cm (39.33 + 5.2 = 44.53, 52 mm deep)

# --- Create surfaces ---
# Assembly outer boundary planes
z_front = openmc.ZPlane(z0=0.0, name="Assembly front face")
z_rear = openmc.ZPlane(z0=assembly_depth, name="Assembly rear face")
x_min = openmc.XPlane(x0=-assembly_half_width, name="Assembly left face (x-)")
x_max = openmc.XPlane(x0=+assembly_half_width, name="Assembly right face (x+)")
y_min = openmc.YPlane(y0=-assembly_half_width, name="Assembly bottom face (y-)")
y_max = openmc.YPlane(y0=+assembly_half_width, name="Assembly top face (y+)")

# Vacuum boundary sphere -- must enclose assembly and source
boundary_sphere = openmc.Sphere(
    x0=0.0, y0=0.0,
    z0=assembly_depth / 2.0,  # centre at assembly midpoint
    r=100.0,                  # generous radius to enclose everything
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# Channel cylinder (along z-axis, centred at x=0, y=0)
channel_inner = openmc.ZCylinder(
    x0=0.0, y0=0.0, r=channel_inner_radius,
    name="Channel inner surface (r=1.4 cm)",
)
channel_outer = openmc.ZCylinder(
    x0=0.0, y0=0.0, r=channel_outer_radius,
    name="Channel outer surface (r=1.5 cm)",
)

# Channel z-extent planes
z_channel_start = z_front  # channel starts at front face
z_channel_end = openmc.ZPlane(z0=channel_z_end, name="Channel end (z=39.07 cm)")

# Detector cavity surfaces
cav_x_lo = openmc.XPlane(x0=-cavity_half_x, name="Cavity x-min")
cav_x_hi = openmc.XPlane(x0=+cavity_half_x, name="Cavity x-max")
cav_y_lo = openmc.YPlane(y0=-cavity_half_y, name="Cavity y-min")
cav_y_hi = openmc.YPlane(y0=+cavity_half_y, name="Cavity y-max")
cav_z_lo = openmc.ZPlane(z0=cavity_z_start, name="Cavity z-min (39.33 cm)")
cav_z_hi = openmc.ZPlane(z0=cavity_z_end, name="Cavity z-max (44.53 cm)")

# --- Define regions ---
# Channel void (air-filled cylindrical bore)
channel_void_region = -channel_inner & +z_channel_start & -z_channel_end

# Channel wall (SS316 annulus)
channel_wall_region = (+channel_inner & -channel_outer
                       & +z_channel_start & -z_channel_end)

# Detector cavity (air-filled rectangular box)
cavity_region = (+cav_x_lo & -cav_x_hi & +cav_y_lo & -cav_y_hi
                 & +cav_z_lo & -cav_z_hi)

# --- Create z-planes for each layer boundary ---
z_planes = {}  # cache z-planes by position to avoid duplicates
for z_val in [0.0, assembly_depth]:
    z_planes[z_val] = z_front if z_val == 0.0 else z_rear

for z_start, z_end, mat, desc in layers:
    for z_val in [z_start, z_end]:
        if z_val not in z_planes:
            z_planes[z_val] = openmc.ZPlane(z0=z_val, name=f"z = {z_val:.2f} cm")

# Also add channel end plane to the cache
z_planes[channel_z_end] = z_channel_end

# --- Build cells for each layer ---
# For layers in the channel region (z < 39.07), we must subtract
# the channel bore and wall from the bulk material.
# For the cavity block layer, we subtract the cavity volume.
layer_cells = []
for i, (z_start, z_end, mat, desc) in enumerate(layers):
    zlo = z_planes[z_start]
    zhi = z_planes[z_end]

    # Base region: slab between z_start and z_end, within assembly cross-section
    base_region = (+zlo & -zhi & +x_min & -x_max & +y_min & -y_max)

    # Subtract channel if this layer overlaps the channel z-range
    if z_start < channel_z_end and z_end > 0.0:
        # Subtract outer channel cylinder (includes wall + bore)
        base_region = base_region & +channel_outer

    # Subtract cavity if this layer contains it
    if z_start <= cavity_z_start and z_end >= cavity_z_end:
        base_region = base_region & ~cavity_region
    elif (z_start < cavity_z_end and z_end > cavity_z_start
          and not (z_start >= cavity_z_start and z_end <= cavity_z_end)):
        # Partial overlap with cavity -- subtract the overlapping part
        base_region = base_region & ~cavity_region

    cell = openmc.Cell(name=f"Layer {i+1}: {desc}", fill=mat, region=base_region)
    layer_cells.append(cell)

# --- Channel void cell (air inside the bore) ---
channel_void_cell = openmc.Cell(
    name="Channel void (air-filled bore, 28 mm ID)",
    fill=air,
    region=channel_void_region,
)

# --- Channel wall cell (1 mm SS316 annulus) ---
channel_wall_cell = openmc.Cell(
    name="Channel wall (1 mm SS316 annulus)",
    fill=ss316,
    region=channel_wall_region,
)

# --- Detector cavity cell (air-filled rectangular box) ---
cavity_cell = openmc.Cell(
    name="Detector cavity (148x48x52 mm, air-filled)",
    fill=air,
    region=cavity_region,
)

# --- Void cell (everything inside boundary sphere but outside assembly) ---
assembly_box_region = (+z_front & -z_rear
                       & +x_min & -x_max & +y_min & -y_max)
void_region = -boundary_sphere & ~assembly_box_region
void_cell = openmc.Cell(name="Void (surrounding air/vacuum)", region=void_region)

# --- Build root universe ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[*layer_cells, channel_void_cell, channel_wall_cell,
           cavity_cell, void_cell],
)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# 14.1 MeV D-T neutron source.
# On-axis: source at (0, 0, -5.3) -- aligned with channel axis
# Off-axis: source at (5.3, 0, -5.3) -- shifted laterally by 5.3 cm in x
#
# The real FNG source has angular and energy dependence from D-T kinematics
# with 300 keV deuterons. The isotropic 14.1 MeV approximation is the
# standard SINBAD benchmark specification.

if args.config == "on-axis":
    source_pos = (0.0, 0.0, -5.3)
    print("Source configuration: ON-AXIS (aligned with channel)")
else:
    source_pos = (5.3, 0.0, -5.3)
    print("Source configuration: OFF-AXIS (shifted 5.3 cm in x)")

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=source_pos)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = source
settings.particles = args.particles
settings.batches = args.batches
settings.photon_transport = False    # neutron-only for this benchmark
settings.output = {"tallies": True}

settings.export_to_xml()


# =============================================================================
# Tallies -- detector cells along channel, in cavity, and behind
# =============================================================================
# Detector positions from the SINBAD specification:
#   Channel positions:     0.25, 12.95, 25.95, 38.65 cm from front face
#   Cavity positions:      39.12 to 43.82 cm (11 foils, ~0.47 cm spacing)
#   Behind-cavity:         46.35, 53.30, 60.05, 66.90, 73.90, 80.60,
#                          87.25, 91.65 cm from front face
#
# We model detector foils as thin cylindrical discs (for channel positions)
# or thin rectangular slabs (for cavity and behind-cavity positions).
# Each detector is a small volume centred on the beam axis.

# Channel detector positions (inside the 28 mm bore)
channel_det_z = [0.25, 12.95, 25.95, 38.65]  # cm from front face

# Cavity detector positions (11 foils evenly spaced)
cavity_det_z = np.linspace(39.12, 43.82, 11).tolist()

# Behind-cavity detector positions
behind_det_z = [46.35, 53.30, 60.05, 66.90, 73.90, 80.60, 87.25, 91.65]

# All detector positions in order
all_det_z = channel_det_z + cavity_det_z + behind_det_z
all_det_labels = (
    [f"channel_{z:.2f}cm" for z in channel_det_z]
    + [f"cavity_{z:.2f}cm" for z in cavity_det_z]
    + [f"behind_{z:.2f}cm" for z in behind_det_z]
)

# Detector cell dimensions
det_half_thickness = 0.05  # cm (1 mm total thickness for most foils)
# Channel detectors: small disc inside the bore (radius 1.0 cm < bore 1.4 cm)
chan_det_radius = 1.0  # cm
# Cavity/behind detectors: small rectangular slab (2 x 2 cm)
det_half_width = 1.0   # cm

# Create detector cells and tallies
detector_cells = []
tallies = openmc.Tallies()

# Energy bins for spectral tallies (100 logarithmic bins, 1 keV to 15 MeV)
energy_bins = np.logspace(np.log10(1.0e3), np.log10(15.0e6), 101)
energy_filter = openmc.EnergyFilter(energy_bins)

for i, (z_pos, label) in enumerate(zip(all_det_z, all_det_labels)):
    # Create bounding surfaces for this detector
    det_z_lo = openmc.ZPlane(z0=z_pos - det_half_thickness,
                             name=f"Det {i+1} z-lo ({label})")
    det_z_hi = openmc.ZPlane(z0=z_pos + det_half_thickness,
                             name=f"Det {i+1} z-hi ({label})")

    if z_pos <= channel_z_end:
        # Channel detector: cylindrical disc inside bore
        det_cyl = openmc.ZCylinder(x0=0.0, y0=0.0, r=chan_det_radius,
                                   name=f"Det {i+1} cylinder ({label})")
        det_region = -det_cyl & +det_z_lo & -det_z_hi
    else:
        # Cavity or behind-cavity detector: rectangular slab
        det_x_lo = openmc.XPlane(x0=-det_half_width,
                                 name=f"Det {i+1} x-lo ({label})")
        det_x_hi = openmc.XPlane(x0=+det_half_width,
                                 name=f"Det {i+1} x-hi ({label})")
        det_y_lo = openmc.YPlane(y0=-det_half_width,
                                 name=f"Det {i+1} y-lo ({label})")
        det_y_hi = openmc.YPlane(y0=+det_half_width,
                                 name=f"Det {i+1} y-hi ({label})")
        det_region = (+det_x_lo & -det_x_hi & +det_y_lo & -det_y_hi
                      & +det_z_lo & -det_z_hi)

    # Detector cell filled with air (activation foils are too thin to
    # perturb transport; the foil material is not modelled explicitly)
    det_cell = openmc.Cell(
        name=f"Detector {i+1}: {label}",
        fill=air,
        region=det_region,
    )
    detector_cells.append(det_cell)

    # --- Flux spectrum tally ---
    cf = openmc.CellFilter([det_cell])
    tally_spectrum = openmc.Tally(name=f"flux_spectrum_{label}")
    tally_spectrum.filters = [cf, energy_filter]
    tally_spectrum.scores = ["flux"]
    tallies.append(tally_spectrum)

    # --- Total flux tally ---
    tally_total = openmc.Tally(name=f"total_flux_{label}")
    tally_total.filters = [cf]
    tally_total.scores = ["flux"]
    tallies.append(tally_total)

# Add detector cells to geometry
# We need to rebuild the root universe with detector cells included.
# Detector cells are tiny volumes embedded in the channel void, cavity,
# or bulk material. They must be added to the universe.
# Note: these small detector volumes will overlap with existing cells.
# In OpenMC, the FIRST cell that contains a point wins, so we place
# detector cells before the other cells in the universe cell list.
root_universe = openmc.Universe(
    name="Root universe",
    cells=[*detector_cells, *layer_cells, channel_void_cell,
           channel_wall_cell, cavity_cell, void_cell],
)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print("=" * 70)
print("FNG-ITER Streaming Experiment -- OpenMC Model")
print("=" * 70)
print(f"  Assembly:    100 x 100 x {assembly_depth} cm")
print(f"  Channel:     {2*channel_inner_radius*10:.0f} mm ID, "
      f"{channel_z_end:.2f} cm long, 1 mm SS316 wall")
print(f"  Cavity:      {2*cavity_half_x*10:.0f} x {2*cavity_half_y*10:.0f} x "
      f"{(cavity_z_end - cavity_z_start)*10:.0f} mm")
print(f"  Source:      14.1 MeV isotropic at {source_pos}")
print(f"  Config:      {args.config}")
print(f"  Detectors:   {len(detector_cells)} positions "
      f"({len(channel_det_z)} channel + {len(cavity_det_z)} cavity + "
      f"{len(behind_det_z)} behind)")
print(f"  Particles:   {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins: {len(energy_bins)-1} logarithmic bins")
print(f"  NOTE: This is a streaming problem. Statistics will be poor at depth")
print(f"        without variance reduction (weight windows).")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
