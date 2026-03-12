#!/usr/bin/env python3
"""
KANT Beryllium Spherical Shell Benchmark - OpenMC Model
========================================================

Source: SINBAD database, KANT experiment at Forschungszentrum Karlsruhe (FZK)
URL: https://www.oecd-nea.org/science/wprs/shielding/sinbad/kant/fzk-be_a.htm

The Karlsruhe Neutron Transmission (KANT) experiment measured neutron leakage
spectra through concentric spherical beryllium shells with a central T(d,n)
neutron source. The experiment exploits the Be-9(n,2n) reaction, which has a
threshold around 1.85 MeV. Because 14.1 MeV neutrons are well above this
threshold, the total neutron leakage through beryllium shells EXCEEDS the
source strength -- i.e., more neutrons come out than go in. This "neutron
multiplication" is a key physics feature relevant to fusion blanket design,
where beryllium is used as a neutron multiplier to breed tritium.

Three shell configurations were tested:
  Shell 1:  5 cm thick  (inner radius 5 cm, outer radius 10 cm)
  Shell 2: 10 cm thick  (inner radius 5 cm, outer radius 15 cm)
  Shell 3: 17 cm thick  (inner radius 5 cm, outer radius 22 cm)

The T(d,n)He-4 source produces 14.1 MeV neutrons. Measurements used NE-213
scintillators, proton recoil proportional counters, and time-of-flight
techniques covering 50 keV to 15 MeV in 178 energy groups. Combined
uncertainty is 8% (1-sigma).
"""

import argparse
import numpy as np
import openmc

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="KANT Beryllium Spherical Shell benchmark (SINBAD)"
)
parser.add_argument(
    "--shell", type=int, choices=[1, 2, 3], default=1,
    help="Shell configuration: 1 (5 cm), 2 (10 cm), or 3 (17 cm) thickness"
)
parser.add_argument(
    "--particles", type=int, default=1_000_000,
    help="Number of source particles (default: 1,000,000)"
)
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Shell geometry parameters (all dimensions in cm)
# ---------------------------------------------------------------------------
INNER_RADIUS = 5.0  # Same for all three configurations

# Map shell number to thickness and outer radius
SHELL_CONFIG = {
    1: {"thickness": 5.0,  "outer_radius": 10.0},
    2: {"thickness": 10.0, "outer_radius": 15.0},
    3: {"thickness": 17.0, "outer_radius": 22.0},
}

config = SHELL_CONFIG[args.shell]
thickness = config["thickness"]
outer_radius = config["outer_radius"]

print(f"=== KANT Be Shell Configuration {args.shell} ===")
print(f"    Inner radius:  {INNER_RADIUS} cm")
print(f"    Outer radius:  {outer_radius} cm")
print(f"    Thickness:     {thickness} cm")
print(f"    Particles:     {args.particles}")

# ---------------------------------------------------------------------------
# Materials
# ---------------------------------------------------------------------------
# Beryllium metal -- the neutron multiplier material
# Density ~1.85 g/cm3, essentially pure Be-9
beryllium = openmc.Material(name="Beryllium")
beryllium.add_nuclide("Be9", 1.0)  # 100% Be-9
beryllium.set_density("g/cm3", 1.85)

materials = openmc.Materials([beryllium])
materials.export_to_xml()

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
# Inner void sphere (contains the point source)
inner_sphere = openmc.Sphere(r=INNER_RADIUS, name="Inner sphere (void boundary)")

# Outer sphere of the beryllium shell
outer_sphere = openmc.Sphere(r=outer_radius, name="Outer sphere (Be shell boundary)")

# Vacuum boundary -- far enough to capture all leakage
boundary_sphere = openmc.Sphere(
    r=outer_radius + 10.0,
    boundary_type="vacuum",
    name="Vacuum boundary"
)

# Region 1: Central void (source location)
inner_void_cell = openmc.Cell(name="Central void (source)")
inner_void_cell.region = -inner_sphere  # Inside inner sphere

# Region 2: Beryllium shell
be_shell_cell = openmc.Cell(name=f"Be shell ({thickness} cm thick)")
be_shell_cell.fill = beryllium
be_shell_cell.region = +inner_sphere & -outer_sphere  # Between inner and outer

# Region 3: Outer void (between shell and vacuum boundary)
outer_void_cell = openmc.Cell(name="Outer void")
outer_void_cell.region = +outer_sphere & -boundary_sphere

# Build the universe
root_universe = openmc.Universe(
    cells=[inner_void_cell, be_shell_cell, outer_void_cell]
)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# ---------------------------------------------------------------------------
# Source: 14.1 MeV D-T neutrons, isotropic point source at origin
# ---------------------------------------------------------------------------
# The T(d,n)He-4 reaction produces monoenergetic 14.1 MeV neutrons.
# We model this as a point source at the center of the sphere.
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))  # Center of sphere
source.angle = openmc.stats.Isotropic()               # Isotropic emission
source.energy = openmc.stats.Discrete([14.1e6], [1.0]) # 14.1 MeV mono-energetic
source.particle = "neutron"

# ---------------------------------------------------------------------------
# Settings: fixed source mode
# ---------------------------------------------------------------------------
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = [source]
settings.particles = args.particles
settings.batches = 10  # 10 batches for statistics
settings.output = {"tallies": True}
settings.export_to_xml()

# ---------------------------------------------------------------------------
# Tallies: neutron current on outer shell surface
# ---------------------------------------------------------------------------
# We tally the surface current on the outer sphere of the Be shell.
# This gives us the neutron leakage spectrum.

# Energy bins: ~100 logarithmic bins from 50 keV to 15 MeV
# (matching the experimental range of 50 keV to 15 MeV)
energy_bins = np.logspace(np.log10(0.05e6), np.log10(15.0e6), 101)  # 100 bins

# Energy filter for spectral resolution
energy_filter = openmc.EnergyFilter(energy_bins)

# Surface filter on the outer sphere of the beryllium shell
surface_filter = openmc.SurfaceFilter([outer_sphere])

# Tally 1: Leakage spectrum (surface current on outer shell surface)
leakage_tally = openmc.Tally(name="leakage_spectrum")
leakage_tally.filters = [surface_filter, energy_filter]
leakage_tally.scores = ["current"]  # Surface current = net leakage

# Tally 2: Total leakage (no energy binning, for quick multiplication check)
total_leakage_tally = openmc.Tally(name="total_leakage")
total_leakage_tally.filters = [surface_filter]
total_leakage_tally.scores = ["current"]

tallies = openmc.Tallies([leakage_tally, total_leakage_tally])
tallies.export_to_xml()

print("\n=== Model exported to XML ===")
print(f"Running OpenMC with {args.particles} particles in {settings.batches} batches...")

# ---------------------------------------------------------------------------
# Run the simulation
# ---------------------------------------------------------------------------
openmc.run()
