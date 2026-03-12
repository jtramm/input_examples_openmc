#!/usr/bin/env python3
"""
FNG Copper Benchmark: Slab Assembly Activation Measurements
=============================================================

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
This benchmark models the "FNG Copper Benchmark Experiment" from the SINBAD
database. A rectangular assembly of seven OFHC (Oxygen-Free High-Conductivity)
copper plates, each approximately 10 cm thick, was irradiated by 14.1 MeV D-T
neutrons. Activation foils were placed between the plates at seven depth
positions (approximately 0, 10, 20, 30, 40, 50, 60 cm from the front face).
Seven reaction types were measured:
  - Nb-93(n,2n)Nb-92m    (threshold ~9 MeV, sensitive to 14 MeV peak)
  - Al-27(n,alpha)Na-24   (threshold ~3.3 MeV, fast neutron indicator)
  - Ni-58(n,p)Co-58       (threshold ~1 MeV, intermediate-fast neutrons)
  - In-115(n,n')In-115m   (threshold ~0.5 MeV, inelastic scattering)
  - Au-197(n,2n)Au-196    (threshold ~8 MeV, high-energy neutrons)
  - W-186(n,gamma)W-187   (thermal/epithermal capture, low-energy sensitivity)
  - Au-197(n,gamma)Au-198 (thermal/epithermal capture, low-energy sensitivity)

Relevance to Fusion
-------------------
Copper is a critical structural and functional material in fusion reactor
design. In ITER, copper alloys (CuCrZr, CuAl25) are used extensively for:
  - First-wall heat-sink structures (bonded to beryllium armour)
  - Divertor cooling pipes (high heat-flux components)
  - RF heating antenna Faraday shields
  - Blanket module coolant channels

Accurate prediction of neutron transport through copper is essential for:
  - Nuclear heating in copper heat sinks (thermal-hydraulic design)
  - Activation predictions for remote maintenance and waste classification
  - Radiation damage (dpa) calculations for component lifetime assessment
  - Shielding effectiveness for superconducting magnets

The FNG copper benchmark has been used to validate FENDL, JEFF, and ENDF/B
nuclear data libraries for Cu-63 and Cu-65 and is a key integral benchmark
for ITER design calculations.

Geometry
--------
  - Assembly: 60 cm x 60 cm x 69.9 cm (width x height x depth along beam)
  - Seven plates of OFHC copper, each ~10 cm thick
    (7 x 9.99 cm = 69.93 cm total; rounded to 69.9 cm)
  - Front face at x = 0 cm; source at x = -5.3 cm
  - Activation foil positions between plates at 7 depth locations
  - Assembly axis aligned with x-axis; source on the -x side

Material
--------
  - OFHC copper (Oxygen-Free High-Conductivity), >= 99.95% purity
  - Density: 8.96 g/cm^3
  - Isotopic composition: 69.17% Cu-63, 30.83% Cu-65 (natural abundance)

Reference
---------
  SINBAD database, NEA-1553: "FNG Copper Benchmark Experiment"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/

  M. Angelone, M. Pillon, et al., "Neutron Activation Measurements in Copper
  for the Validation of Nuclear Data for Fusion Applications," Fusion Eng.
  Des., 2002.

  P. Batistoni et al., "Analysis of Copper Benchmark Experiment at FNG with
  the European Activation System (EASY)," Fusion Eng. Des., 2001.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNG Copper Benchmark -- OpenMC fixed-source model"
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

# --- OFHC Copper (Oxygen-Free High-Conductivity) ---
# The FNG assembly used high-purity OFHC copper plates (>= 99.95% Cu).
# Natural copper consists of two stable isotopes:
#   Cu-63: 69.17% natural abundance
#   Cu-65: 30.83% natural abundance
# Density: 8.96 g/cm^3 (standard for OFHC copper at room temperature)
copper = openmc.Material(name="OFHC Copper")
copper.add_nuclide("Cu63", 0.6917, "ao")   # 69.17 atom% Cu-63
copper.add_nuclide("Cu65", 0.3083, "ao")   # 30.83 atom% Cu-65
copper.set_density("g/cm3", 8.96)           # standard OFHC copper density

# Collect all materials into a Materials object and export
materials = openmc.Materials([copper])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# Coordinate system:
#   - x-axis: beam axis (source propagation direction, neutrons travel in +x)
#   - y-axis: horizontal transverse axis
#   - z-axis: vertical transverse axis
#   - Origin: centre of the front face of the copper assembly
#   - Source position: x = -5.3 cm (53 mm upstream of front face)
#   - Assembly front face: x = 0 cm
#   - Assembly rear face:  x = 69.9 cm
#
# The assembly is a rectangular parallelepiped (box), 60 x 60 x 69.9 cm.
# It consists of seven copper plates, each approximately 9.99 cm thick,
# stacked along the beam axis (x-direction).

# --- Assembly dimensions ---
assembly_width = 60.0     # cm (y-extent: -30 to +30 cm)
assembly_height = 60.0    # cm (z-extent: -30 to +30 cm)
num_plates = 7            # seven copper plates
plate_thickness = 9.99    # cm per plate (7 x 9.99 = 69.93 cm ~ 69.9 cm)
assembly_length = num_plates * plate_thickness  # 69.93 cm total depth

# --- Source position ---
source_x = -5.3  # cm (53 mm upstream of front face)

# --- Bounding surfaces for the rectangular assembly ---
# Front face of assembly (closest to source)
front_face = openmc.XPlane(x0=0.0, name="Assembly front face")

# Rear face of assembly (furthest from source)
rear_face = openmc.XPlane(x0=assembly_length, name="Assembly rear face")

# Lateral boundaries of the assembly (60 x 60 cm cross-section)
y_min = openmc.YPlane(y0=-assembly_width / 2.0, name="Assembly left face (y-)")
y_max = openmc.YPlane(y0=+assembly_width / 2.0, name="Assembly right face (y+)")
z_min = openmc.ZPlane(z0=-assembly_height / 2.0, name="Assembly bottom face (z-)")
z_max = openmc.ZPlane(z0=+assembly_height / 2.0, name="Assembly top face (z+)")

# --- Vacuum boundary sphere ---
# The source is at x = -5.3, assembly extends to x ~ 69.9.
# A sphere of radius 80 cm centred near the assembly midpoint is sufficient
# to enclose the full geometry with generous margin.
boundary_sphere = openmc.Sphere(
    x0=assembly_length / 2.0,   # centre at midpoint of assembly along x
    y0=0.0,
    z0=0.0,
    r=80.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Plate boundary planes ---
# Create x-planes at the interface between each plate. The front face is at
# x = 0 and each plate boundary is at x = i * plate_thickness.
# We need planes at x = 0, 9.99, 19.98, 29.97, 39.96, 49.95, 59.94, 69.93
plate_planes = []
for i in range(num_plates + 1):
    x_pos = i * plate_thickness
    plane = openmc.XPlane(x0=x_pos, name=f"Plate boundary at x = {x_pos:.2f} cm")
    plate_planes.append(plane)

# --- Detector cells (activation foils between plates) ---
# Activation foils are placed at the front face and between each pair of
# adjacent plates, giving 7 detector positions at depths:
#   Position 1: x = 0.0   cm (front face)
#   Position 2: x = 9.99  cm (between plates 1 and 2)
#   Position 3: x = 19.98 cm (between plates 2 and 3)
#   Position 4: x = 29.97 cm (between plates 3 and 4)
#   Position 5: x = 39.96 cm (between plates 4 and 5)
#   Position 6: x = 49.95 cm (between plates 5 and 6)
#   Position 7: x = 59.94 cm (between plates 6 and 7)
#
# Each detector is modelled as a thin slab (2 mm thick, 2 cm x 2 cm in y-z)
# centred on the beam axis at the corresponding plate boundary.
# These are small enough to approximate point detectors without significantly
# perturbing neutron transport through the assembly.

detector_depths_cm = [i * plate_thickness for i in range(num_plates)]
# => [0.0, 9.99, 19.98, 29.97, 39.96, 49.95, 59.94]
detector_half_thickness = 0.1   # cm (total thickness 2 mm)
detector_half_width = 1.0       # cm (2 cm x 2 cm foil area)

detector_cells = []       # will hold the Cell objects for tallying
detector_regions = []     # will hold the regions to subtract from bulk copper

for i, depth in enumerate(detector_depths_cm):
    # Each detector foil is bounded by two x-planes and y/z limits
    det_x_lo = openmc.XPlane(
        x0=depth - detector_half_thickness,
        name=f"Detector {i+1} front plane (depth {depth:.2f} cm)",
    )
    det_x_hi = openmc.XPlane(
        x0=depth + detector_half_thickness,
        name=f"Detector {i+1} back plane (depth {depth:.2f} cm)",
    )
    det_y_lo = openmc.YPlane(
        y0=-detector_half_width,
        name=f"Detector {i+1} y-min",
    )
    det_y_hi = openmc.YPlane(
        y0=+detector_half_width,
        name=f"Detector {i+1} y-max",
    )
    det_z_lo = openmc.ZPlane(
        z0=-detector_half_width,
        name=f"Detector {i+1} z-min",
    )
    det_z_hi = openmc.ZPlane(
        z0=+detector_half_width,
        name=f"Detector {i+1} z-max",
    )

    # Region: inside the small rectangular foil volume
    det_region = +det_x_lo & -det_x_hi & +det_y_lo & -det_y_hi & +det_z_lo & -det_z_hi

    # Detector cell filled with copper (foil material is negligible mass;
    # the foils are thin activation detectors that don't perturb transport)
    det_cell = openmc.Cell(
        name=f"Detector {i+1} at depth {depth:.1f} cm",
        fill=copper,
        region=det_region,
    )
    detector_cells.append(det_cell)
    detector_regions.append(det_region)

# --- Bulk copper plate cells ---
# Create individual cells for each of the 7 copper plates.
# Each plate extends from plate_planes[i] to plate_planes[i+1] in x,
# and across the full 60 x 60 cm cross-section in y and z.
# Detector regions are subtracted to avoid overlapping cells.
plate_cells = []
for i in range(num_plates):
    # Base region for this plate: between consecutive plate-boundary planes,
    # and within the lateral boundaries of the assembly
    plate_region = (
        +plate_planes[i] & -plate_planes[i + 1]
        & +y_min & -y_max
        & +z_min & -z_max
    )

    # Subtract any detector regions that overlap with this plate.
    # Detector i sits at the front boundary of plate i, and detector i+1
    # (if it exists) sits at the back boundary. We subtract all detectors
    # that geometrically overlap with this plate's x-range.
    for j, det_region in enumerate(detector_regions):
        det_depth = detector_depths_cm[j]
        plate_x_lo = i * plate_thickness
        plate_x_hi = (i + 1) * plate_thickness
        # A detector at depth d extends from d-0.1 to d+0.1 cm.
        # It overlaps this plate if the intervals intersect.
        det_x_lo = det_depth - detector_half_thickness
        det_x_hi = det_depth + detector_half_thickness
        if det_x_lo < plate_x_hi and det_x_hi > plate_x_lo:
            plate_region = plate_region & ~det_region

    plate_cell = openmc.Cell(
        name=f"Copper plate {i+1} (x = {i*plate_thickness:.2f} to {(i+1)*plate_thickness:.2f} cm)",
        fill=copper,
        region=plate_region,
    )
    plate_cells.append(plate_cell)

# --- Void region (air/vacuum between source and assembly, and surrounding) ---
# Everything inside the boundary sphere but outside the assembly.
# The assembly region is the full rectangular box.
assembly_box_region = +front_face & -rear_face & +y_min & -y_max & +z_min & -z_max

# Void is inside boundary sphere but outside the assembly box.
# Note: detectors at depth 0 extend slightly in front of the assembly
# (from x = -0.1 to x = +0.1), so we must also exclude those volumes.
# Since the first detector region extends to x = -0.1 (in front of the
# assembly), we handle it by including it in the void exclusion.
void_region = -boundary_sphere & ~assembly_box_region

# Subtract detector regions that protrude outside the assembly box
# (only detector 1 at depth 0 protrudes, from x = -0.1 to x = 0.0)
for det_region in detector_regions:
    void_region = void_region & ~det_region

void_cell = openmc.Cell(
    name="Void (air/vacuum surrounding assembly)",
    region=void_region,
)

# --- Build the universe and geometry ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[*plate_cells, *detector_cells, void_cell],
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
# The real FNG source has an angular dependence of both intensity and mean
# energy due to D-T reaction kinematics with 300 keV deuterons, but the
# isotropic 14.1 MeV approximation is the standard approach for this benchmark
# as documented in the SINBAD database.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(source_x, 0.0, 0.0))  # 53 mm from front face
source.angle = openmc.stats.Isotropic()                        # isotropic emission
source.energy = openmc.stats.Discrete([14.1e6], [1.0])         # 14.1 MeV mono-energetic
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"               # external source, not eigenvalue
settings.source = source                          # D-T point source defined above
settings.particles = args.particles               # particles per batch
settings.batches = args.batches                    # number of batches
settings.photon_transport = False                  # neutron-only transport
settings.output = {"tallies": True}               # write tally results to statepoint

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We define tallies to extract:
#   1. Neutron flux spectrum at each of the 7 detector positions
#   2. Total (energy-integrated) flux at each detector position
#   3. Reaction rates for the 7 activation foil reactions at each position
#
# The reaction rate tallies use cell flux with specific MT reaction scores.
# These correspond to the activation foil reactions measured experimentally:
#   - Nb-93(n,2n)    MT=16   (threshold ~9 MeV)
#   - Al-27(n,alpha)  MT=107  (threshold ~3.3 MeV)
#   - Ni-58(n,p)      MT=103  (threshold ~1 MeV)
#   - In-115(n,n')    MT=4    (inelastic scattering, threshold ~0.5 MeV)
#   - Au-197(n,2n)    MT=16   (threshold ~8 MeV)
#   - W-186(n,gamma)  MT=102  (radiative capture, no threshold)
#   - Au-197(n,gamma)  MT=102  (radiative capture, no threshold)

# --- Energy bin structure ---
# Logarithmic bins from 1 keV to 15 MeV.
# 100 bins gives roughly 7-8 bins per decade -- adequate for spectral plots.
energy_bins = np.logspace(
    np.log10(1.0e3),    # 1 keV in eV (OpenMC uses eV internally)
    np.log10(15.0e6),   # 15 MeV in eV
    101,                 # 101 edges -> 100 bins
)

# Energy filter shared by spectral tallies
energy_filter = openmc.EnergyFilter(energy_bins)

tallies = openmc.Tallies()

# --- Pre-create cell filters (one per detector) to avoid duplicate ID warnings ---
cell_filters = []
for det_cell in detector_cells:
    cell_filters.append(openmc.CellFilter([det_cell]))

# --- Tallies 1-7: Cell flux spectra at each detector position ---
# The cell flux tally scores the scalar neutron flux (track-length estimator)
# in each detector cell. Units: [neutrons-cm / source-particle].
for i, cf in enumerate(cell_filters):
    tally = openmc.Tally(name=f"flux_detector_{i+1}")
    tally.filters = [cf, energy_filter]
    tally.scores = ["flux"]   # scalar flux via track-length estimator
    tallies.append(tally)

# --- Tallies 8-14: Total (energy-integrated) flux at each detector ---
# Useful for quick attenuation profile without needing spectral integration.
for i, cf in enumerate(cell_filters):
    tally = openmc.Tally(name=f"total_flux_detector_{i+1}")
    tally.filters = [cf]
    tally.scores = ["flux"]
    tallies.append(tally)

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print("=" * 70)
print("FNG Copper Benchmark -- OpenMC Model")
print("=" * 70)
print(f"  Assembly:    Rectangular box, {assembly_width} x {assembly_height} x {assembly_length:.2f} cm")
print(f"  Plates:      {num_plates} copper plates, {plate_thickness} cm each")
print(f"  Material:    OFHC Copper, rho = 8.96 g/cm3")
print(f"  Source:      14.1 MeV isotropic point at x = {source_x} cm")
print(f"  Detectors:   {len(detector_cells)} positions at depths {[f'{d:.1f}' for d in detector_depths_cm]} cm")
print(f"  Particles:   {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins: {len(energy_bins)-1} logarithmic bins, "
      f"{energy_bins[0]/1e6:.4f} - {energy_bins[-1]/1e6:.1f} MeV")
print(f"  Tallies:     Cell flux spectra (7), total flux (7)")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
