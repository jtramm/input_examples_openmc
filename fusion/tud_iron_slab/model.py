#!/usr/bin/env python3
"""
TUD Iron Slab Benchmark
========================

Facility
--------
The experiment was performed at the Technische Universitaet Dresden (TUD),
Germany. TUD operates a pulsed 14 MeV D-T neutron generator used for fusion
neutronics shielding experiments. The facility is designed for spectral
measurements of neutron and photon fields behind shield assemblies to validate
nuclear data and transport codes for fusion applications.

Experiment
----------
This benchmark models the "TUD Iron Slab" experiment from the SINBAD database.
A large iron slab assembly (100 x 100 cm frontal area, 30 cm thick) built from
20 x 10 x 5 cm iron blocks was irradiated with 14 MeV D-T neutrons. Spectral
neutron and photon flux measurements were performed behind the slab using an
NE213 liquid scintillation detector. Three assembly configurations were tested:

  - A0: Solid slab (no gap) -- baseline attenuation measurement
  - A1: 5 cm vertical gap located 10 cm from the slab centre
  - A2: 5 cm vertical gap located 20 cm from the slab centre

The gap configurations (A1, A2) simulate streaming paths through iron shield
walls, relevant to penetrations in fusion reactor shielding (beam ports,
diagnostic channels, coolant pipes). By varying the gap offset from the beam
axis, the experiment quantifies the angular dependence of radiation streaming.

Relevance to Fusion
-------------------
Iron and steel are the primary structural materials for fusion reactor
components (vacuum vessel, blanket support structures, biological shield).
Accurate prediction of neutron deep-penetration through iron is essential for:
  - Biological shield design (dose rate behind thick steel walls)
  - Activation calculations for maintenance planning
  - Nuclear heating in structural components
  - Streaming through gaps and penetrations (a major design challenge)

This benchmark has been used to validate FENDL-1, EFF-2, and other nuclear
data libraries for iron, and is particularly valuable for testing streaming
corrections in deterministic transport codes.

Geometry
--------
  - Iron slab: 100 cm x 100 cm frontal area, 30 cm thick
  - Built from 20 x 10 x 5 cm iron blocks (density ~7.874 g/cm3)
  - Three configurations:
      A0: Solid slab (no gap)
      A1: 5 cm wide vertical gap, offset 10 cm from centre
      A2: 5 cm wide vertical gap, offset 20 cm from centre
  - Coordinate system:
      x-axis: beam axis (source-to-detector direction)
      y-axis: horizontal transverse (gap offset direction)
      z-axis: vertical
      Origin: centre of the front face of the iron slab
  - D-T point source at x = -19 cm (19 cm from front face)
  - Slab front face at x = 0 cm
  - Slab rear face at x = 30 cm
  - Detector at x = 330 cm (300 cm behind rear face)
  - Total source-to-detector distance: 349 cm
  - D-beam angle: 74 degrees relative to source-slab axis

Source
------
  - 14 MeV D-T pulsed neutron generator
  - Approximated as isotropic point source at 14.1 MeV
  - Source time profile: exp[-(t / 1.4 ns)^2] (Gaussian pulse)

Material
--------
  - Natural iron (Fe), density 7.874 g/cm3
  - Chemical composition: essentially pure iron (trace impurities neglected)

Reference
---------
  SINBAD database: "TUD Iron Slab Experiment"
  https://www.oecd-nea.org/science/wprs/shielding/sinbad/tud_fe/tufe-abs.htm

  Reference calculations: MCNP-4A with FENDL-1 and EFF-2 nuclear data.
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser(
    description="TUD Iron Slab Benchmark -- OpenMC fixed-source model"
)
parser.add_argument(
    "--config",
    type=str,
    choices=["A0", "A1", "A2"],
    default="A0",
    help="Slab configuration: A0 (solid), A1 (gap at 10 cm), A2 (gap at 20 cm). Default: A0",
)
parser.add_argument(
    "--particles",
    type=int,
    default=1_000_000,
    help="Number of source particles per batch (default: 1,000,000)",
)
parser.add_argument(
    "--batches",
    type=int,
    default=10,
    help="Number of batches (default: 10)",
)
args = parser.parse_args()


# =============================================================================
# Configuration parameters
# =============================================================================
# Gap offset from the slab centre (y-axis) and gap width for each configuration.
# A0 has no gap. A1 and A2 have a 5 cm wide vertical gap at different offsets.
gap_configs = {
    "A0": {"has_gap": False, "gap_centre_y": None, "gap_width": None},
    "A1": {"has_gap": True,  "gap_centre_y": 10.0, "gap_width": 5.0},   # 5 cm gap at y = 10 cm
    "A2": {"has_gap": True,  "gap_centre_y": 20.0, "gap_width": 5.0},   # 5 cm gap at y = 20 cm
}

config = gap_configs[args.config]
print(f"Configuration: {args.config}")
if config["has_gap"]:
    print(f"  Gap: {config['gap_width']} cm wide, centred at y = {config['gap_centre_y']} cm")
else:
    print("  Gap: None (solid slab)")


# =============================================================================
# Materials
# =============================================================================

# --- Natural iron ---
# The TUD slab was built from iron blocks. We use natural iron at the
# standard density of 7.874 g/cm3. The SINBAD data package includes a
# detailed chemical analysis showing essentially pure iron with negligible
# trace impurities (C, Mn, Si, etc.) that we omit here.
iron = openmc.Material(name="Natural iron")
iron.add_element("Fe", 1.0, "wo")   # 100 wt% iron (natural isotopic mix)
iron.set_density("g/cm3", 7.874)     # standard iron density

# Collect all materials and export to XML
materials = openmc.Materials([iron])
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# Coordinate system:
#   x-axis: beam axis (source propagation direction, perpendicular to slab)
#   y-axis: horizontal transverse (direction of gap offset)
#   z-axis: vertical
#   Origin: centre of the front face of the iron slab
#
# Key positions along x-axis:
#   Source:     x = -19 cm
#   Front face: x = 0 cm
#   Rear face:  x = 30 cm (slab is 30 cm thick)
#   Detector:   x = 330 cm (300 cm behind rear face)

# --- Slab dimensions ---
slab_half_width = 50.0    # cm (100 cm / 2) -- half-extent in y and z
slab_thickness = 30.0     # cm -- extent in x
source_distance = 19.0    # cm -- source-to-front-face distance

# --- Bounding surfaces for the iron slab ---
# Front face of slab (closest to source)
front_face = openmc.XPlane(x0=0.0, name="Slab front face")

# Rear face of slab (furthest from source, detector side)
rear_face = openmc.XPlane(x0=slab_thickness, name="Slab rear face")

# Lateral boundaries of the slab (y-direction: horizontal transverse)
slab_left = openmc.YPlane(y0=-slab_half_width, name="Slab left edge (y=-50)")
slab_right = openmc.YPlane(y0=slab_half_width, name="Slab right edge (y=+50)")

# Top and bottom boundaries of the slab (z-direction: vertical)
slab_bottom = openmc.ZPlane(z0=-slab_half_width, name="Slab bottom edge (z=-50)")
slab_top = openmc.ZPlane(z0=slab_half_width, name="Slab top edge (z=+50)")

# --- Vacuum boundary ---
# A large box or sphere to enclose the entire problem. We use a sphere
# centred roughly at the midpoint between source and detector.
# Source at x=-19, detector at x=330 => midpoint ~155 cm.
# Radius of 400 cm is sufficient to enclose everything.
boundary_sphere = openmc.Sphere(
    x0=slab_thickness / 2.0,   # centre at midpoint of slab for simplicity
    y0=0.0,
    z0=0.0,
    r=400.0,
    boundary_type="vacuum",
    name="Vacuum boundary sphere",
)

# --- Define the full slab region (before subtracting any gap) ---
# The slab occupies the rectangular box: 0 < x < 30, -50 < y < 50, -50 < z < 50
full_slab_region = (
    +front_face & -rear_face &
    +slab_left & -slab_right &
    +slab_bottom & -slab_top
)

# --- Build cells depending on configuration ---
cells = []

if config["has_gap"]:
    # For A1 and A2: create a vertical gap (void channel) through the slab.
    # The gap is a rectangular void extending the full height (z) and depth (x)
    # of the slab, but only 5 cm wide in y, offset from centre.
    #
    # Gap boundaries in y:
    #   left edge:  gap_centre_y - gap_width/2
    #   right edge: gap_centre_y + gap_width/2
    gap_y_min = config["gap_centre_y"] - config["gap_width"] / 2.0
    gap_y_max = config["gap_centre_y"] + config["gap_width"] / 2.0

    # Planes defining the gap in y-direction
    gap_left_plane = openmc.YPlane(
        y0=gap_y_min,
        name=f"Gap left edge (y={gap_y_min} cm)",
    )
    gap_right_plane = openmc.YPlane(
        y0=gap_y_max,
        name=f"Gap right edge (y={gap_y_max} cm)",
    )

    # Gap region: full slab depth and height, but narrow in y
    gap_region = (
        +front_face & -rear_face &
        +gap_left_plane & -gap_right_plane &
        +slab_bottom & -slab_top
    )

    # Gap cell: void (no material fill) -- represents the air channel
    gap_cell = openmc.Cell(
        name=f"Gap ({args.config}: {config['gap_width']} cm at y={config['gap_centre_y']} cm)",
        region=gap_region,
        # No fill => void
    )
    cells.append(gap_cell)

    # Iron slab cell: full slab region minus the gap
    iron_slab_region = full_slab_region & ~gap_region
    iron_cell = openmc.Cell(
        name=f"Iron slab ({args.config}, with gap subtracted)",
        fill=iron,
        region=iron_slab_region,
    )
    cells.append(iron_cell)

else:
    # A0: Solid slab, no gap -- simple rectangular iron block
    iron_cell = openmc.Cell(
        name="Iron slab (A0, solid -- no gap)",
        fill=iron,
        region=full_slab_region,
    )
    cells.append(iron_cell)

# --- Void region (everything outside the slab but inside boundary) ---
# This includes the space between source and slab, and behind the slab
# out to the detector position and beyond.
void_region = -boundary_sphere & ~full_slab_region
void_cell = openmc.Cell(
    name="Void (air/vacuum surrounding the slab assembly)",
    region=void_region,
)
cells.append(void_cell)

# --- Build the root universe and geometry ---
root_universe = openmc.Universe(name="Root universe", cells=cells)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The TUD D-T source is approximated as an isotropic point source of 14.1 MeV
# neutrons located 19 cm upstream of the slab front face, on the beam axis.
#
# In our coordinate system, the source is at x = -19 cm, y = 0, z = 0.
#
# The actual D-T generator had a pulsed time structure with a Gaussian profile
# (FWHM ~ 1.4 ns), but since we are computing time-integrated spectra, the
# pulse shape is not modelled here. The beam angle of 74 degrees causes a
# slight angular and energy variation of the source neutrons due to D-T
# kinematics, but we use the standard isotropic 14.1 MeV approximation.

source = openmc.IndependentSource()
source.space = openmc.stats.Point(xyz=(-source_distance, 0.0, 0.0))  # 19 cm from front face
source.angle = openmc.stats.Isotropic()                                # isotropic emission
source.energy = openmc.stats.Discrete([14.1e6], [1.0])                # 14.1 MeV mono-energetic
source.particle = "neutron"                                            # neutron source


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"           # fixed-source mode (not eigenvalue)
settings.source = source                      # D-T point source defined above
settings.particles = args.particles           # particles per batch (from CLI)
settings.batches = args.batches               # number of batches (from CLI)
settings.photon_transport = True              # enable photon transport (experiment measured photons too)
settings.output = {"tallies": True}           # write tally results to statepoint

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# We define tallies to extract the transmitted neutron and photon spectra
# on the rear face of the iron slab. The experiment measured:
#   - Neutron spectra above ~1 MeV (NE213)
#   - Photon spectra above ~0.2 MeV (NE213)
#   - Neutron spectra from 30 keV to 2.3 MeV (proportional counters)
#   - Time-of-arrival spectra (not modelled here -- requires time binning)

# --- Energy bin structure ---
# Logarithmic bins from 0.1 MeV to 15 MeV, covering the main range of
# interest for the transmitted spectrum. 100 bins gives ~7 bins per decade.
energy_bins = np.logspace(
    np.log10(0.1e6),    # 0.1 MeV in eV (lower bound)
    np.log10(15.0e6),   # 15 MeV in eV (upper bound, above D-T peak)
    101,                 # 101 edges => 100 bins
)

# Energy filter shared by neutron and photon tallies
energy_filter = openmc.EnergyFilter(energy_bins)

# Particle filters to separate neutron and photon contributions
neutron_filter = openmc.ParticleFilter(["neutron"])
photon_filter = openmc.ParticleFilter(["photon"])

# Surface filter on the rear face of the slab
rear_surface_filter = openmc.SurfaceFilter([rear_face])

tallies = openmc.Tallies()

# --- Tally 1: Transmitted neutron spectrum (surface current on rear face) ---
# This tally captures all neutrons crossing the rear surface of the slab.
# The current score counts particles per unit source particle.
neutron_tally = openmc.Tally(name="transmitted_neutron_spectrum")
neutron_tally.filters = [rear_surface_filter, energy_filter, neutron_filter]
neutron_tally.scores = ["current"]   # particles crossing the surface
tallies.append(neutron_tally)

# --- Tally 2: Transmitted photon spectrum (surface current on rear face) ---
# Photons are generated by (n,gamma) reactions in iron and by inelastic
# scattering. The photon spectrum behind the slab is also of interest.
photon_tally = openmc.Tally(name="transmitted_photon_spectrum")
photon_tally.filters = [rear_surface_filter, energy_filter, photon_filter]
photon_tally.scores = ["current"]   # particles crossing the surface
tallies.append(photon_tally)

# --- Tally 3: Total transmitted neutron current (energy-integrated) ---
# Useful for quick comparison of total attenuation across configurations.
total_neutron_tally = openmc.Tally(name="total_transmitted_neutrons")
total_neutron_tally.filters = [rear_surface_filter, neutron_filter]
total_neutron_tally.scores = ["current"]
tallies.append(total_neutron_tally)

# --- Tally 4: Total transmitted photon current (energy-integrated) ---
total_photon_tally = openmc.Tally(name="total_transmitted_photons")
total_photon_tally.filters = [rear_surface_filter, photon_filter]
total_photon_tally.scores = ["current"]
tallies.append(total_photon_tally)

tallies.export_to_xml()


# =============================================================================
# Summary
# =============================================================================
print("=" * 70)
print("TUD Iron Slab Benchmark -- OpenMC Model")
print("=" * 70)
print(f"  Configuration: {args.config}")
if config["has_gap"]:
    print(f"  Gap:           {config['gap_width']} cm wide at y = {config['gap_centre_y']} cm")
else:
    print(f"  Gap:           None (solid slab)")
print(f"  Slab:          {2*slab_half_width:.0f} x {2*slab_half_width:.0f} x {slab_thickness:.0f} cm")
print(f"  Material:      Natural iron, rho = 7.874 g/cm3")
print(f"  Source:        14.1 MeV isotropic point at x = -{source_distance} cm")
print(f"  Particles:     {args.particles:,} per batch x {args.batches} batches")
print(f"  Energy bins:   {len(energy_bins)-1} logarithmic bins, "
      f"{energy_bins[0]/1e6:.1f} - {energy_bins[-1]/1e6:.1f} MeV")
print(f"  Tallies:       Neutron spectrum (1), photon spectrum (1), "
      f"totals (2)")
print(f"  Photon transport: Enabled")
print("=" * 70)
print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
print("Run with: openmc")
