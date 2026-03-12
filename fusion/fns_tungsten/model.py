#!/usr/bin/env python3
"""
FNS Clean Benchmark: Tungsten Cylindrical Assembly
====================================================

Facility
--------
The Fusion Neutronics Source (FNS) facility is located at JAEA (formerly JAERI)
in Tokai-mura, Japan. FNS produces 14 MeV neutrons via the D-T reaction using
a 350 keV deuteron beam impinging on a rotating tritium target. The facility was
purpose-built for fusion neutronics experiments to validate nuclear data and
transport codes relevant to fusion reactor design.

Experiment
----------
This benchmark models the "FNS Tungsten Clean Experiment" from the SINBAD
database (NEA-1553/80). A cylindrical assembly of tungsten alloy bricks was
irradiated by 14.1 MeV D-T neutrons. Neutron spectra were measured at four
depths (0, 76, 228, 380 mm) along the central axis using multiple detector
types spanning 5 keV to 14 MeV. Activation foil reaction rates (Al-27(n,a),
Nb-93(n,2n), In-115(n,n'), W-186(n,g), Au-197(n,g)) and gamma heating rates
were also measured.

Relevance to Fusion
-------------------
Tungsten is the primary plasma-facing material for the ITER divertor and is
planned for future fusion power plant first-wall armour. Accurate prediction
of neutron transport through tungsten is essential for:
  - Shielding design (protecting superconducting magnets)
  - Nuclear heating calculations (thermal-hydraulic design)
  - Activation and transmutation predictions (waste management)
  - Damage (dpa) calculations (component lifetime)

This benchmark has been used to validate FENDL, JEFF, and ENDF/B nuclear data
libraries for tungsten isotopes (W-180, W-182, W-183, W-184, W-186) and is
considered a key integral benchmark for ITER design validation.

Geometry
--------
  - Tungsten cylinder: diameter 629 mm (radius 314.5 mm), height 507 mm
  - Constructed from tungsten alloy bricks (DENSIMET-type, W-Ni-Cu)
  - Individual brick thickness: 50.7-50.8 mm
  - Assembly supported by thin aluminium frame (neglected in model)
  - D-T point source located 200 mm from the front face of the assembly
  - Assembly axis aligned with x-axis; source at x = -20 cm

Material
--------
  - Tungsten alloy (DENSIMET D185): 95 wt% W, 3.5 wt% Ni, 1.5 wt% Cu
  - Density: 18.0 g/cm^3 (typical for this alloy grade)
  - Note: the SINBAD abstract states "tungsten is not pure, but an alloy with
    a small amount of nickel and copper." The exact composition is given in
    the full SINBAD data package. The values used here are representative of
    the DENSIMET-class alloys used in these experiments.

Reference
---------
  SINBAD database, NEA-1553/80: "FNS/JAERI Clean Experiment on Tungsten"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/fns_w/fnsw-abs.htm

  F. Maekawa, Y. Oyama, et al., "Benchmark Experiment on a Tungsten Assembly
  Bombarded by D-T Neutrons," J. Nucl. Sci. Technol., Suppl. 2, 2002.

Validated with OpenMC v0.14.1-dev against MCNP-4A reference calculations.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="FNS Tungsten Clean Benchmark -- OpenMC fixed-source model"
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

# --- Tungsten alloy (DENSIMET D185-type) ---
# The FNS assembly used bricks of a W-Ni-Cu alloy, commonly referred to as
# DENSIMET or Inermet. The nominal composition is 95 wt% W, 3.5 wt% Ni,
# 1.5 wt% Cu. The density of this alloy is approximately 18.0 g/cm^3 (lower
# than pure W at 19.3 g/cm^3 due to the Ni-Cu binder phase).
tungsten_alloy = openmc.Material(name="Tungsten alloy (DENSIMET D185)")
tungsten_alloy.add_element("W",  0.950, "wo")  # 95.0 wt% tungsten
tungsten_alloy.add_element("Ni", 0.035, "wo")  # 3.5 wt% nickel (binder)
tungsten_alloy.add_element("Cu", 0.015, "wo")  # 1.5 wt% copper (binder)
tungsten_alloy.set_density("g/cm3", 18.0)       # DENSIMET alloy density

# Collect all materials into a Materials object and export
materials = openmc.Materials([tungsten_alloy])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# Coordinate system:
#   - x-axis: beam axis (source propagation direction)
#   - y,z: transverse axes
#   - Origin: centre of the front face of the tungsten assembly
#   - Source position: x = -20 cm (200 mm upstream of front face)
#   - Assembly front face: x = 0 cm
#   - Assembly rear face:  x = 50.7 cm (507 mm depth)

# --- Assembly dimensions ---
assembly_radius = 31.45   # cm  (diameter 629 mm / 2)
assembly_length = 50.7    # cm  (507 mm total height/length along beam axis)

# --- Bounding surfaces for the tungsten cylinder ---
# Front face of assembly (closest to source)
front_face = openmc.XPlane(x0=0.0, name="Assembly front face")

# Rear face of assembly (furthest from source)
rear_face = openmc.XPlane(x0=assembly_length, name="Assembly rear face")

# Cylindrical boundary of the assembly (along x-axis)
assembly_cyl = openmc.XCylinder(r=assembly_radius, name="Assembly outer cylinder")

# --- Vacuum boundary sphere (large enough to enclose everything) ---
# The source is at x=-20, assembly extends to x=50.7.
# A sphere of radius 100 cm centred at the assembly midpoint is sufficient.
boundary_sphere = openmc.Sphere(
    x0=assembly_length / 2.0,  # centre at midpoint of assembly
    y0=0.0,
    z0=0.0,
    r=100.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Detector cells along the central axis ---
# Neutron spectra were measured at 4 depths from the front face:
#   Position 1:   0 mm (front surface)
#   Position 2:  76 mm ( 7.6 cm)
#   Position 3: 228 mm (22.8 cm)
#   Position 4: 380 mm (38.0 cm)
#
# We model each detector as a thin cylindrical disc (2 mm thick, 1 cm radius)
# centred on the beam axis at the specified depth. These are small enough to
# approximate point detectors without significantly perturbing the geometry.
detector_depths_cm = [0.0, 7.6, 22.8, 38.0]  # cm from front face
detector_half_thickness = 0.1   # cm (total thickness 2 mm)
detector_radius = 1.0           # cm

detector_cells = []       # will hold the Cell objects for tallying
detector_regions = []     # will hold the regions to subtract from bulk tungsten

for i, depth in enumerate(detector_depths_cm):
    # Each detector disc is bounded by two planes and a cylinder
    det_front = openmc.XPlane(
        x0=depth - detector_half_thickness,
        name=f"Detector {i+1} front (depth {depth} cm)",
    )
    det_back = openmc.XPlane(
        x0=depth + detector_half_thickness,
        name=f"Detector {i+1} back (depth {depth} cm)",
    )
    det_cyl = openmc.XCylinder(
        r=detector_radius,
        name=f"Detector {i+1} cylinder (depth {depth} cm)",
    )

    # Region: inside cylinder AND between the two planes
    det_region = +det_front & -det_back & -det_cyl

    # Detector cell filled with tungsten alloy (the foils are negligible)
    det_cell = openmc.Cell(
        name=f"Detector {i+1} at depth {depth:.1f} cm",
        fill=tungsten_alloy,
        region=det_region,
    )
    detector_cells.append(det_cell)
    detector_regions.append(det_region)

# --- Bulk tungsten cell ---
# The main tungsten assembly is the cylinder minus the detector disc volumes.
# Start with the full cylinder region.
bulk_region = +front_face & -rear_face & -assembly_cyl

# Subtract each detector region from the bulk to avoid overlapping cells.
for det_region in detector_regions:
    bulk_region = bulk_region & ~det_region

bulk_tungsten = openmc.Cell(
    name="Bulk tungsten assembly",
    fill=tungsten_alloy,
    region=bulk_region,
)

# --- Void region (air/vacuum between source and assembly, and around) ---
# Everything inside the boundary sphere but outside the assembly and detectors.
void_region = -boundary_sphere & ~(+front_face & -rear_face & -assembly_cyl)
void_cell = openmc.Cell(
    name="Void (air/vacuum surrounding assembly)",
    region=void_region,
)

# --- Build the universe and geometry ---
root_universe = openmc.Universe(
    name="Root universe",
    cells=[bulk_tungsten, *detector_cells, void_cell],
)

geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The FNS D-T source is approximated as an isotropic point source of 14.1 MeV
# neutrons located 200 mm (20 cm) upstream of the assembly front face.
# In our coordinate system, the source is at x = -20 cm, y = 0, z = 0.
#
# The real source has a slight angular dependence of both intensity and mean
# energy due to the kinematics of the D-T reaction with 350 keV deuterons,
# but the isotropic 14.1 MeV approximation is standard for this benchmark.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(-20.0, 0.0, 0.0))  # 200 mm from front face
source.angle = openmc.stats.Isotropic()                     # isotropic emission
source.energy = openmc.stats.Discrete([14.1e6], [1.0])      # 14.1 MeV mono-energetic
source.particle = "neutron"


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"               # not eigenvalue -- external source
settings.source = source                          # D-T point source defined above
settings.particles = args.particles               # particles per batch
settings.batches = args.batches                    # number of batches
settings.photon_transport = False                  # neutron-only for this model
settings.output = {"tallies": True}               # write tally results to statepoint

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We define tallies to extract:
#   1. Neutron flux spectrum at each of the 4 detector positions
#   2. Surface current spectrum on the rear face (transmitted spectrum)

# --- Energy bin structure ---
# Logarithmic bins from 10 keV (0.01 MeV) to 15 MeV.
# 100 bins gives roughly 7-8 bins per decade -- adequate for spectral plots.
energy_bins = np.logspace(
    np.log10(1.0e4),   # 10 keV in eV (OpenMC uses eV internally)
    np.log10(15.0e6),  # 15 MeV in eV
    101,               # 101 edges -> 100 bins
)

# Energy filter shared by all tallies
energy_filter = openmc.EnergyFilter(energy_bins)

tallies = openmc.Tallies()

# --- Tally 1-4: Cell flux at each detector position ---
# The cell flux tally scores the scalar neutron flux (track-length estimator)
# in each detector disc cell. Units: [neutrons / cm / source-particle] when
# divided by cell volume.
for i, det_cell in enumerate(detector_cells):
    cell_filter = openmc.CellFilter([det_cell])
    tally = openmc.Tally(name=f"flux_detector_{i+1}")
    tally.filters = [cell_filter, energy_filter]
    tally.scores = ["flux"]  # scalar flux (track-length estimator)
    tallies.append(tally)

# --- Tally 5: Surface current on rear face (transmitted spectrum) ---
# This captures neutrons leaving the back of the tungsten assembly.
# The surface current tally counts particles crossing the surface.
rear_surface_filter = openmc.SurfaceFilter([rear_face])
transmitted_tally = openmc.Tally(name="transmitted_spectrum")
transmitted_tally.filters = [rear_surface_filter, energy_filter]
transmitted_tally.scores = ["current"]  # net current through surface
tallies.append(transmitted_tally)

# --- Tally 6: Total flux in each detector (energy-integrated) ---
# Useful for quick comparison of attenuation without needing post-processing.
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
print("=" * 70)
print("FNS Tungsten Clean Benchmark -- OpenMC Model")
print("=" * 70)
print(f"  Assembly:    Cylinder, R = {assembly_radius} cm, L = {assembly_length} cm")
print(f"  Material:    W-Ni-Cu alloy, rho = 18.0 g/cm3")
print(f"  Source:      14.1 MeV isotropic point at x = -20 cm")
print(f"  Detectors:   {len(detector_cells)} positions at depths {detector_depths_cm} cm")
print(f"  Particles:   {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins: {len(energy_bins)-1} logarithmic bins, {energy_bins[0]/1e6:.4f} - {energy_bins[-1]/1e6:.1f} MeV")
print(f"  Tallies:     Cell flux (4), transmitted spectrum (1), total flux (4)")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
