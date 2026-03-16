#!/usr/bin/env python3
"""
ARC Reactor Fusion Neutronics Benchmark -- Full Toroidal Sector Model
======================================================================

Facility / Design
-----------------
The ARC (Affordable, Robust, Compact) reactor is a compact, high-field
tokamak concept developed at MIT. It exploits high-temperature
superconducting (HTS) magnets (REBCO) to achieve very high on-axis
toroidal fields (~9.2 T) in a physically small device (major radius
R = 3.3 m, minor radius a = 1.13 m). The blanket uses molten FLiBe
(Li2BeF4), which simultaneously serves as:

  - Tritium breeder (via Li-6(n,t)He-4 and Li-7(n,n't)He-4)
  - Neutron moderator and multiplier (via Be-9)
  - Primary coolant

The ARC design targets ~525 MW of D-T fusion power.

Geometry
--------
This model builds a 20-degree toroidal sector (360/18 = 20 degrees,
corresponding to one of 18 TF coil periods) with reflective boundary
conditions on the sector-cutting planes. The cross-section is
constructed from concentric torus surfaces (ZTorus in OpenMC):

  1. Plasma region (void):    minor radius 113 cm
  2. First wall (Inconel 718): 113 - 116 cm (3 cm thick)
  3. FLiBe blanket:            116 - 196 cm (80 cm thick)
  4. Vacuum vessel (Inconel):  196 - 206 cm (10 cm thick)
  5. Neutron shield (borated steel): 206 - 216 cm (10 cm thick)
  6. TF coil (simplified composite): 216 - 246 cm (30 cm thick)

All torus surfaces share a major radius of 330 cm and are centred on
the z-axis.

Source
------
D-T fusion neutrons at 14.1 MeV, distributed uniformly in a cylindrical
annulus around the magnetic axis: R in [300, 360] cm, Z in [-50, 50] cm,
phi in [0, 20 degrees]. This approximates the peaked fusion reaction
rate near the plasma core.

Tallies
-------
  - Tritium Breeding Ratio (TBR): (n,Xt) reactions in FLiBe blanket
  - Nuclear heating in FLiBe blanket
  - Neutron flux at the TF coil region
  - Energy-dependent neutron spectrum leaking past the shield

Weight Windows
--------------
Deep penetration through 80 cm of FLiBe + 10 cm steel requires variance
reduction for good statistics at the TF coil. The --generate-ww flag
activates FW-CADIS weight window generation via the random ray solver.

References
----------
  B. N. Sorbom et al., "ARC: A compact, high-field, fusion nuclear
  science facility and target technology for demonstrating net
  fusion gain," Fusion Eng. Des., vol. 100, pp. 378-405, 2015.

  A. Q. Kuang et al., "Conceptual design study for heat exhaust
  management in the ARC fusion pilot plant," Fusion Eng. Des.,
  vol. 137, pp. 221-242, 2018.

  J. H. Bae, G. P. Peterson, T. Shimwell, "Benchmarking of tritium
  breeding ratio calculations for the ARC reactor," Fusion Sci.
  Technol., 2022.
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
    description="ARC Reactor Fusion Neutronics Benchmark -- OpenMC model"
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
    help="Generate weight windows via FW-CADIS random ray before the CE run",
)
args = parser.parse_args()


# =============================================================================
# ARC device parameters (all lengths in cm)
# =============================================================================
# Major radius of the plasma magnetic axis
R_MAJOR = 330.0  # cm (3.3 m)

# Minor radius of the last closed flux surface (LCFS)
A_PLASMA = 113.0  # cm (1.13 m)

# Plasma elongation and triangularity (used for reference; this model
# uses circular cross-section tori for simplicity)
KAPPA = 1.84   # elongation
DELTA = 0.52   # triangularity

# Radial build: each layer defined by its outer minor radius
R_FW_OUTER = 116.0    # first wall outer surface (3 cm Inconel 718)
R_BLANKET_OUTER = 196.0  # FLiBe blanket outer surface (80 cm blanket)
R_VV_OUTER = 206.0    # vacuum vessel outer surface (10 cm Inconel 718)
R_SHIELD_OUTER = 216.0  # neutron shield outer surface (10 cm borated steel)
R_COIL_OUTER = 246.0  # TF coil outer surface (30 cm composite)

# Number of TF coils and sector angle
N_TF_COILS = 18
SECTOR_ANGLE_DEG = 360.0 / N_TF_COILS  # 20 degrees

# Operating temperature of FLiBe blanket (for cross section Doppler broadening)
TEMPERATURE_K = 600.0  # K


# =============================================================================
# Materials
# =============================================================================

# --- FLiBe (Li2BeF4) molten salt blanket ---
# Li2BeF4: 2 Li atoms, 1 Be atom, 4 F atoms per formula unit
# Natural lithium: 7.5 at% Li-6, 92.5 at% Li-7
# Density at ~600 C: approximately 1.94 g/cm^3
flibe = openmc.Material(name="FLiBe (Li2BeF4)")
flibe.set_density("g/cm3", 1.94)
flibe.add_nuclide("Li6", 2 * 0.075, "ao")   # 2 Li atoms * 7.5% Li-6
flibe.add_nuclide("Li7", 2 * 0.925, "ao")   # 2 Li atoms * 92.5% Li-7
flibe.add_nuclide("Be9", 1.0, "ao")          # 1 Be atom
flibe.add_nuclide("F19", 4.0, "ao")          # 4 F atoms
flibe.temperature = TEMPERATURE_K

# --- Inconel 718 (first wall and vacuum vessel) ---
# A nickel-chromium superalloy with high strength at elevated temperatures.
# Widely used in aerospace and nuclear applications.
# Density: 8.19 g/cm^3
# Composition (weight percent):
#   Ni: 52.5%, Cr: 19.0%, Fe: 18.5%, Nb: 5.1%, Mo: 3.0%,
#   Ti: 0.9%, Al: 0.5%, C: 0.04%, others: ~0.46%
inconel718 = openmc.Material(name="Inconel 718")
inconel718.set_density("g/cm3", 8.19)
inconel718.add_element("Ni", 0.525, "wo")   # nickel (balance element)
inconel718.add_element("Cr", 0.190, "wo")   # chromium
inconel718.add_element("Fe", 0.185, "wo")   # iron
inconel718.add_element("Nb", 0.051, "wo")   # niobium (columbium)
inconel718.add_element("Mo", 0.030, "wo")   # molybdenum
inconel718.add_element("Ti", 0.009, "wo")   # titanium
inconel718.add_element("Al", 0.005, "wo")   # aluminium
inconel718.add_nuclide("C12", 0.0004 * 0.9893, "wo")  # carbon-12
inconel718.add_nuclide("C13", 0.0004 * 0.0107, "wo")  # carbon-13
# Remaining 0.46% approximated as Mn + Si
inconel718.add_element("Mn", 0.003, "wo")   # manganese
inconel718.add_element("Si", 0.0016, "wo")  # silicon

# --- Borated steel neutron shield ---
# SS316-type stainless steel with 2 wt% natural boron added for thermal
# neutron absorption. Boron-10 has a very large (n,alpha) cross section.
# Density: approximately 7.65 g/cm^3 (slightly less than pure SS316 due
# to lower-density boron inclusions).
borated_steel = openmc.Material(name="Borated Steel Shield")
borated_steel.set_density("g/cm3", 7.65)
borated_steel.add_element("Fe", 0.6370, "wo")  # iron (balance)
borated_steel.add_element("Cr", 0.1666, "wo")  # chromium 17% * 0.98
borated_steel.add_element("Ni", 0.1176, "wo")  # nickel 12% * 0.98
borated_steel.add_element("Mo", 0.0245, "wo")  # molybdenum 2.5% * 0.98
borated_steel.add_element("Mn", 0.0196, "wo")  # manganese 2% * 0.98
borated_steel.add_element("Si", 0.0098, "wo")  # silicon 1% * 0.98
borated_steel.add_element("N", 0.0049, "wo")   # nitrogen 0.5% * 0.98
borated_steel.add_element("B", 0.0200, "wo")   # 2 wt% natural boron

# --- TF coil composite ---
# Simplified homogenised mixture representing the HTS winding pack:
# copper stabiliser, REBCO tape, and stainless steel structure.
# Approximate composition: 40% Cu, 50% SS316, 10% Hastelloy (approx as Ni)
# Density: ~8.0 g/cm^3 (volume-weighted average)
coil_composite = openmc.Material(name="TF Coil Composite")
coil_composite.set_density("g/cm3", 8.0)
coil_composite.add_element("Cu", 0.40, "wo")  # copper stabiliser
coil_composite.add_element("Fe", 0.325, "wo") # SS316 iron fraction
coil_composite.add_element("Cr", 0.085, "wo") # SS316 chromium fraction
coil_composite.add_element("Ni", 0.16, "wo")  # SS316 nickel + Hastelloy
coil_composite.add_element("Mo", 0.015, "wo") # SS316 molybdenum fraction
coil_composite.add_element("Mn", 0.015, "wo") # SS316 manganese fraction

# Collect all materials
materials = openmc.Materials([flibe, inconel718, borated_steel, coil_composite])
materials.cross_sections = "/data/endfb-viii.0-hdf5/cross_sections.xml"


# =============================================================================
# Geometry
# =============================================================================
# The tokamak is built from concentric ZTorus surfaces, all sharing the same
# major radius (R_MAJOR) and centred at the origin. The axis of revolution
# is the z-axis. Each torus surface has circular cross-section (b = c = minor
# radius). While the real ARC plasma is elongated (kappa ~ 1.84), this
# simplified CSG model uses circular cross-sections for all layers. A more
# detailed model would use elliptical cross-sections or CAD geometry.

# --- Torus surfaces (concentric shells) ---
# ZTorus parameters: a = major radius, b = minor radius (z-direction),
#                    c = minor radius (radial direction)
# For circular cross-section: b = c = minor_radius

plasma_boundary = openmc.ZTorus(
    a=R_MAJOR, b=A_PLASMA, c=A_PLASMA,
    name="Plasma boundary (LCFS)"
)
fw_outer = openmc.ZTorus(
    a=R_MAJOR, b=R_FW_OUTER, c=R_FW_OUTER,
    name="First wall outer surface"
)
blanket_outer = openmc.ZTorus(
    a=R_MAJOR, b=R_BLANKET_OUTER, c=R_BLANKET_OUTER,
    name="FLiBe blanket outer surface"
)
vv_outer = openmc.ZTorus(
    a=R_MAJOR, b=R_VV_OUTER, c=R_VV_OUTER,
    name="Vacuum vessel outer surface"
)
shield_outer = openmc.ZTorus(
    a=R_MAJOR, b=R_SHIELD_OUTER, c=R_SHIELD_OUTER,
    name="Neutron shield outer surface"
)
coil_outer = openmc.ZTorus(
    a=R_MAJOR, b=R_COIL_OUTER, c=R_COIL_OUTER,
    name="TF coil outer surface"
)

# --- Sector-cutting planes ---
# We model a 20-degree toroidal sector using two planes that pass through
# the z-axis. Both planes have reflective boundary conditions, simulating
# a full torus via symmetry.
#
# Plane 1: the xz-plane (y = 0). Particles reflecting off this plane
#          stay within the sector.
# Plane 2: a plane rotated 20 degrees about the z-axis from the xz-plane.
#          Its normal vector is (sin(phi), -cos(phi), 0) where phi = 20 deg.
phi_rad = math.radians(SECTOR_ANGLE_DEG)

sector_plane_1 = openmc.YPlane(
    y0=0.0,
    boundary_type="reflective",
    name="Sector plane 1 (xz-plane)"
)
sector_plane_2 = openmc.Plane(
    a=math.sin(phi_rad),
    b=-math.cos(phi_rad),
    c=0.0,
    d=0.0,
    boundary_type="reflective",
    name="Sector plane 2 (rotated 20 deg)"
)

# --- Outer bounding sphere ---
# Enclose the entire geometry. The outermost torus extends to
# R_MAJOR + R_COIL_OUTER = 330 + 246 = 576 cm from the z-axis.
# Use a sphere of radius 600 cm as a vacuum boundary.
bounding_sphere = openmc.Sphere(
    r=600.0,
    boundary_type="vacuum",
    name="Outer vacuum boundary"
)

# --- Sector region ---
# The sector is defined as the region between the two cutting planes
# where y >= 0 (plane 1) and on the correct side of plane 2.
# For a sector starting at phi=0 and ending at phi=20 deg,
# the region is: +sector_plane_1 AND -sector_plane_2
# (positive half-space of y=0 plane AND negative half-space of the
#  rotated plane)
sector_region = +sector_plane_1 & +sector_plane_2

# --- Cells ---
# Each cell is the region between consecutive torus surfaces, intersected
# with the sector region and bounded by the outer sphere.

# 1. Plasma (void) -- inside the plasma boundary torus
plasma_cell = openmc.Cell(name="Plasma (void)")
plasma_cell.region = -plasma_boundary & sector_region & -bounding_sphere

# 2. First wall (Inconel 718) -- between plasma boundary and FW outer
fw_cell = openmc.Cell(name="First Wall (Inconel 718)")
fw_cell.fill = inconel718
fw_cell.region = +plasma_boundary & -fw_outer & sector_region & -bounding_sphere

# 3. FLiBe blanket -- between FW outer and blanket outer
blanket_cell = openmc.Cell(name="FLiBe Blanket")
blanket_cell.fill = flibe
blanket_cell.region = +fw_outer & -blanket_outer & sector_region & -bounding_sphere

# 4. Vacuum vessel (Inconel 718) -- between blanket outer and VV outer
vv_cell = openmc.Cell(name="Vacuum Vessel (Inconel 718)")
vv_cell.fill = inconel718
vv_cell.region = +blanket_outer & -vv_outer & sector_region & -bounding_sphere

# 5. Neutron shield (borated steel) -- between VV outer and shield outer
shield_cell = openmc.Cell(name="Neutron Shield (Borated Steel)")
shield_cell.fill = borated_steel
shield_cell.region = +vv_outer & -shield_outer & sector_region & -bounding_sphere

# 6. TF coil (composite) -- between shield outer and coil outer
coil_cell = openmc.Cell(name="TF Coil (Composite)")
coil_cell.fill = coil_composite
coil_cell.region = +shield_outer & -coil_outer & sector_region & -bounding_sphere

# 7. External void -- outside the coil but inside the bounding sphere,
#    within the sector. This is the region between the outermost torus
#    and the bounding sphere (and also the region inside the torus hole).
external_void = openmc.Cell(name="External void")
external_void.region = +coil_outer & sector_region & -bounding_sphere

# Build the root universe from all cells
root_universe = openmc.Universe(
    cells=[
        plasma_cell, fw_cell, blanket_cell, vv_cell,
        shield_cell, coil_cell, external_void
    ]
)
geometry = openmc.Geometry(root_universe)


# =============================================================================
# Source definition: D-T fusion neutrons (14.1 MeV)
# =============================================================================
# The D-T fusion reaction produces a 14.1 MeV neutron and a 3.5 MeV alpha
# particle. In a tokamak, the fusion rate peaks near the magnetic axis and
# decreases toward the plasma edge. We approximate this with a uniform
# cylindrical annulus source centred on the magnetic axis.
#
# Spatial distribution:
#   R (cylindrical radius): uniform in [300, 360] cm (centred on R_MAJOR=330)
#   phi (toroidal angle):   uniform in [0, 20 deg] (within sector)
#   Z (vertical):           uniform in [-50, 50] cm (near midplane)
#
# Energy: monoenergetic at 14.1 MeV
# Direction: isotropic (fusion neutrons are born isotropically in the lab frame
#            to a very good approximation)

source = openmc.IndependentSource()
source.space = openmc.stats.CylindricalIndependent(
    r=openmc.stats.Uniform(a=300.0, b=360.0),
    phi=openmc.stats.Uniform(a=0.0, b=phi_rad),
    z=openmc.stats.Uniform(a=-50.0, b=50.0),
)
source.energy = openmc.stats.Discrete([14.1e6], [1.0])
source.angle = openmc.stats.Isotropic()


# =============================================================================
# Settings
# =============================================================================
settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.source = [source]
settings.particles = args.particles
settings.batches = args.batches
settings.temperature = {"default": TEMPERATURE_K}

# Photon transport is disabled by default. Enable if nuclear heating
# from photons is needed (increases runtime significantly).
settings.photon_transport = False


# =============================================================================
# Tallies
# =============================================================================
tallies = openmc.Tallies()

# --- Tally 1: Tritium Breeding Ratio (TBR) ---
# The TBR is defined as the number of tritium atoms produced per source
# neutron (i.e., per fusion reaction). For a self-sustaining reactor,
# TBR > 1 is required. The (n,Xt) score counts all reactions that produce
# at least one triton, weighted by the number of tritons produced.
tbr_tally = openmc.Tally(name="TBR")
tbr_tally.filters = [openmc.MaterialFilter([flibe])]
tbr_tally.scores = ["(n,Xt)"]
tallies.append(tbr_tally)

# --- Tally 2: Nuclear heating in FLiBe blanket ---
# The 'heating' score gives the total nuclear heating (kinetic energy
# deposited by all reaction products) in eV per source particle.
# This can be converted to watts using the known fusion power.
heating_tally = openmc.Tally(name="Blanket Heating")
heating_tally.filters = [openmc.MaterialFilter([flibe])]
heating_tally.scores = ["heating"]
tallies.append(heating_tally)

# --- Tally 3: Neutron flux at TF coil ---
# Radiation damage and heating in the HTS magnets is a critical design
# constraint. This tally gives the neutron flux (track-length estimator)
# in the coil region per source neutron.
coil_flux_tally = openmc.Tally(name="Coil Neutron Flux")
coil_flux_tally.filters = [openmc.CellFilter([coil_cell])]
coil_flux_tally.scores = ["flux"]
tallies.append(coil_flux_tally)

# --- Tally 4: Energy-dependent neutron spectrum past the shield ---
# 175-group VITAMIN-J-like energy structure (simplified to 20 equal
# lethargy bins from 1e-5 eV to 15 MeV) to capture the spectrum shape.
energy_bins = np.logspace(np.log10(1.0e-5), np.log10(15.0e6), 21)
spectrum_tally = openmc.Tally(name="Shield Leakage Spectrum")
spectrum_tally.filters = [
    openmc.CellFilter([coil_cell]),
    openmc.EnergyFilter(energy_bins),
]
spectrum_tally.scores = ["flux"]
tallies.append(spectrum_tally)

# --- Tally 5: Neutron flux in shield region ---
# Useful for evaluating the shielding performance.
shield_flux_tally = openmc.Tally(name="Shield Neutron Flux")
shield_flux_tally.filters = [openmc.CellFilter([shield_cell])]
shield_flux_tally.scores = ["flux"]
tallies.append(shield_flux_tally)


# =============================================================================
# Build the model
# =============================================================================
def build_model():
    """Construct and return the OpenMC Model object.

    Returns
    -------
    openmc.Model
        The complete ARC reactor model ready for simulation.
    """
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )
    return model


def generate_weight_windows(model):
    """Generate weight windows using FW-CADIS via the random ray solver.

    This function creates a multigroup random ray version of the model,
    runs it to generate weight windows, and then loads those weight
    windows into the original continuous-energy model.

    The FW-CADIS (Forward-Weighted Consistent Adjoint Driven Importance
    Sampling) method computes importance maps that are used to set weight
    window bounds. This is essential for deep-penetration problems like
    neutron transport through the thick FLiBe blanket and steel shield
    to reach the TF coil region.

    Parameters
    ----------
    model : openmc.Model
        The continuous-energy model to generate weight windows for.

    Returns
    -------
    openmc.Model
        The original model with weight windows loaded from the
        weight_windows.h5 file.
    """
    print("=" * 70)
    print("GENERATING WEIGHT WINDOWS VIA FW-CADIS RANDOM RAY")
    print("=" * 70)

    # Deep copy the model to avoid modifying the original
    model_rr = copy.deepcopy(model)

    # Step 1: Convert to multigroup using CASMO-4 energy groups
    # This generates a multigroup cross section library (mgxs.h5) by
    # running material-wise infinite medium calculations.
    print("\n[Step 1] Converting to multigroup (CASMO-4 groups)...")
    model_rr.convert_to_multigroup(
        method="material_wise",
        groups="CASMO-4",
        nparticles=3000,
        mgxs_path="mgxs.h5",
    )

    # Step 2: Convert to random ray solver
    # This sets appropriate ray tracing parameters based on the geometry.
    print("[Step 2] Converting to random ray solver...")
    model_rr.convert_to_random_ray()

    # Step 3: Set up source region mesh for random ray
    # A 10 cm mesh covers the bounding box of the geometry.
    print("[Step 3] Setting up source region mesh (10 cm resolution)...")
    bbox = model_rr.geometry.bounding_box
    ll = bbox.lower_left
    ur = bbox.upper_right

    n_x = max(1, int((ur[0] - ll[0]) / 10.0))
    n_y = max(1, int((ur[1] - ll[1]) / 10.0))
    n_z = max(1, int((ur[2] - ll[2]) / 10.0))

    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (n_x, n_y, n_z)
    ww_mesh.lower_left = list(ll)
    ww_mesh.upper_right = list(ur)

    root = model_rr.geometry.root_universe
    model_rr.settings.random_ray["source_region_meshes"] = [
        (ww_mesh, [root])
    ]

    # Step 4: Configure weight window generator
    # FW-CADIS will compute adjoint fluxes and forward fluxes to derive
    # optimal weight window bounds.
    print("[Step 4] Configuring FW-CADIS weight window generator...")
    model_rr.settings.batches = 200
    model_rr.settings.inactive = 100
    model_rr.settings.particles = 500

    wwg = openmc.WeightWindowGenerator(
        mesh=ww_mesh,
        method="fw_cadis",
        max_realizations=model_rr.settings.batches,
    )
    model_rr.settings.weight_window_generators = wwg

    # Step 5: Run the random ray model
    print("[Step 5] Running random ray model for WW generation...")
    model_rr.run()

    # Step 6: Load weight windows into the CE model
    print("\n[Step 6] Loading weight windows into CE model...")
    ww = openmc.hdf5_to_wws("weight_windows.h5")
    model.settings.weight_windows = ww

    print("Weight window generation complete.")
    print("=" * 70)

    return model


# =============================================================================
# Main entry point
# =============================================================================
if __name__ == "__main__":
    model = build_model()

    if args.generate_ww:
        # Generate weight windows first, then run with them
        model = generate_weight_windows(model)

    # Export the model to XML and run
    model.export_to_model_xml()
    print(f"\nModel exported. Running with {args.particles:,} particles/batch "
          f"x {args.batches} batches...")
    model.run()
    print("Simulation complete.")
