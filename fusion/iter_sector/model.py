#!/usr/bin/env python3
"""
Simplified ITER Tokamak 40-Degree Sector Benchmark
=====================================================

Device Overview
---------------
ITER (International Thermonuclear Experimental Reactor) is the world's largest
tokamak, currently under construction at the CEA Cadarache site near
Saint-Paul-les-Durance in southern France. It is designed to demonstrate the
scientific and technical feasibility of fusion energy by achieving Q >= 10
(500 MW fusion power from 50 MW of auxiliary heating).

The tokamak confines a deuterium-tritium plasma in a toroidal magnetic field
produced by 18 superconducting toroidal field (TF) coils. The D-T fusion
reaction produces 14.1 MeV neutrons and 3.5 MeV alpha particles. Of the
500 MW total fusion power, approximately 400 MW is carried by neutrons
(80%) and 100 MW by alpha particles which heat the plasma.

This model creates a simplified 40-degree sector (1/9 of the full torus,
corresponding to 2 TF coil periods) with reflective boundary conditions on
the sector faces. The geometry uses concentric toroidal shells (ZTorus
surfaces) to represent the major radial build components.

Geometry (Simplified CSG)
-------------------------
The model builds concentric toroidal shells around the plasma, cut by two
sector planes at phi=0 and phi=40 degrees. Each shell represents a major
ITER component:

    Component               Minor radius range (cm)    Thickness (cm)
    ---------               ----------------------     --------------
    Plasma (void)           0 -- 200                   200
    Be armor (first wall)   200 -- 201                 1
    Cu heat sink            201 -- 202                 1
    SS316LN structure       202 -- 204                 2
    Blanket (SS+H2O)        204 -- 249                 45
    Gap (void)              249 -- 253                 4
    VV inner shell          253 -- 259                 6
    VV fill (SS+H2O)        259 -- 289                 30
    VV outer shell          289 -- 295                 6
    Gap (void)              295 -- 300                 5
    Thermal shield (SS304)  300 -- 302                 2
    Gap (void)              302 -- 310                 8
    TF coil case (SS316LN)  310 -- 315                 5
    TF winding pack (hom.)  315 -- 340                 25
    TF coil case (SS316LN)  340 -- 345                 5

All toroidal shells share the major radius R0 = 620 cm and use circular
cross-sections (b = c in ZTorus) for simplicity. The actual ITER plasma is
elongated (kappa = 1.85) and triangular (delta = 0.48), but circular
cross-sections are adequate for this shielding benchmark.

Materials
---------
Six materials are defined:

  1. Beryllium armor (1.85 g/cc): Plasma-facing armor tiles on the first
     wall. Be is chosen for its low atomic number, which minimizes plasma
     contamination from sputtering.

  2. Copper heat sink (8.96 g/cc): CuCrZr alloy (approximated as pure Cu)
     bonded to the Be armor, providing thermal conductivity to remove the
     surface heat flux (~0.5 MW/m2 average, up to 5 MW/m2 at divertor).

  3. SS316LN (7.93 g/cc): Low-nitrogen variant of austenitic stainless
     steel, the primary structural material for the blanket, vacuum vessel,
     and TF coil case. Composition: Fe 62.5%, Cr 17.5%, Ni 12.5%, Mo 2.5%,
     Mn 2%, Si 1%, N 0.15%, C 0.03% by weight.

  4. Blanket/VV fill -- homogenized SS316LN + H2O (60/40 vol%): The ITER
     baseline shielding blanket is non-breeding and consists of SS316LN
     blocks with water cooling channels. The 60/40 steel/water mixture is a
     standard approximation in ITER neutronics analyses.

  5. SS304 thermal shield (7.93 g/cc): Simplified as SS316LN for this model.
     The thermal shield sits between the vacuum vessel and TF coils,
     intercepting thermal radiation to protect the cryogenic coils.

  6. TF coil winding pack (homogenized): The winding pack contains Nb3Sn
     superconductor in a cable-in-conduit conductor (CICC) design with Cu
     stabilizer, SS316LN conduit, and epoxy insulation. Simplified here as
     a 50/30/15/5 vol% mixture of SS316LN/Cu/Nb/Sn.

D-T Source
----------
The neutron source is distributed uniformly within the plasma volume using
cylindrical coordinates. The R distribution spans 460--780 cm (R0 +/- 160 cm,
using 80% of the plasma minor radius to avoid sampling near the edge where
the density drops). The Z distribution spans -160 to +160 cm. The phi
distribution covers the 40-degree sector. All neutrons are emitted at
14.1 MeV (monoenergetic D-T fusion neutrons) with isotropic angular
distribution.

Weight Windows
--------------
The ~140 cm of shielding between the plasma and the TF coils attenuates the
neutron flux by many orders of magnitude. Weight windows are essential for
obtaining meaningful statistics at the TF coil location. This model supports
automatic weight window generation via the FW-CADIS random ray method.

Reference
---------
ITER Organization, "ITER Research Plan within the Staged Approach,"
ITR-18-003, 2018.

M.J. Loughlin et al., "ITER Nuclear Analysis Strategy and Requirements,"
Fusion Sci. Technol., vol. 56, no. 2, pp. 566-572, 2009.

R. Juarez et al., "Shutdown dose rate calculations at ITER equatorial
port using MCNP," Fusion Eng. Des., vol. 100, pp. 501-506, 2015.

Usage
-----
    python model.py [--particles N] [--batches N] [--generate-ww]
                    [--use-ww FILE] [--small-tallies]
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
    description="Simplified ITER 40-degree sector -- OpenMC fixed-source model"
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
    help="Generate weight windows via FW-CADIS random ray before running",
)
parser.add_argument(
    "--use-ww",
    type=str,
    default=None,
    help="Path to weight_windows.h5 file to load",
)
parser.add_argument(
    "--small-tallies",
    action="store_true",
    help="Use reduced tally set (total heating and coil flux only)",
)
args = parser.parse_args()


# =============================================================================
# Device parameters
# =============================================================================
R0 = 620.0         # Major radius [cm]
a_plasma = 200.0   # Plasma minor radius [cm]
sector_angle = 40.0  # Sector angle [degrees] (= 360/9 for 18 TF coils)


# =============================================================================
# Materials
# =============================================================================

# --- 1. Beryllium armor ---
# Beryllium tiles cover the ITER first wall, providing a low-Z plasma-facing
# surface that minimizes radiation losses from impurity line emission. The
# tiles are ~10 mm thick and bonded to the Cu heat sink via hot isostatic
# pressing (HIP). Be is also an excellent neutron multiplier via the
# Be-9(n,2n) reaction, which is important for tritium breeding in future
# reactor blankets (but not exploited in ITER's non-breeding baseline).
beryllium = openmc.Material(name="Beryllium armor")
beryllium.add_element("Be", 1.0, "ao")
beryllium.set_density("g/cm3", 1.85)

# --- 2. Copper heat sink ---
# The CuCrZr alloy heat sink (approximated as pure Cu) provides the thermal
# path between the Be armor and the actively-cooled steel structure. CuCrZr
# retains adequate strength and thermal conductivity up to ~350 C under
# irradiation. Natural Cu composition: 69.17% Cu-63, 30.83% Cu-65.
copper = openmc.Material(name="Copper heat sink")
copper.add_nuclide("Cu63", 0.6917, "ao")
copper.add_nuclide("Cu65", 0.3083, "ao")
copper.set_density("g/cm3", 8.96)

# --- 3. SS316LN structural steel ---
# SS316L(N) is the primary structural alloy for ITER's in-vessel components
# and vacuum vessel. The "LN" designation indicates low carbon (< 0.03%)
# and controlled nitrogen (0.06--0.08%) for improved mechanical properties
# after welding and irradiation. The composition below is simplified.
#
# Density: 7.93 g/cm^3 at room temperature. The actual operating temperature
# varies: ~200 C for the blanket structure, ~100 C for the VV (water-cooled),
# and ~4 K for the TF coil case (cryogenic). Density variation with
# temperature is neglected in this simplified model.
ss316ln = openmc.Material(name="SS316LN")
ss316ln.add_element("Fe", 0.6252, "wo")   # iron, balance
ss316ln.add_element("Cr", 0.1750, "wo")   # chromium
ss316ln.add_element("Ni", 0.1250, "wo")   # nickel
ss316ln.add_element("Mo", 0.0250, "wo")   # molybdenum
ss316ln.add_element("Mn", 0.0200, "wo")   # manganese
ss316ln.add_element("Si", 0.0100, "wo")   # silicon
ss316ln.add_element("N",  0.0015, "wo")   # nitrogen
ss316ln.add_element("C",  0.0003, "wo")   # carbon
ss316ln.set_density("g/cm3", 7.93)

# --- 4. Blanket / VV fill: homogenized SS316LN + H2O (60/40 vol%) ---
# The ITER shielding blanket modules consist of SS316LN blocks with water
# cooling channels. A homogenized 60% steel / 40% water mixture by volume
# is the standard representation in ITER neutronics analyses. The same
# composition is used for the vacuum vessel internal fill, which also
# contains SS plates and water for shielding and cooling.
#
# Volume fractions:
#   SS316LN: 60% at 7.93 g/cc -> partial density 4.758 g/cc
#   H2O:     40% at 1.00 g/cc -> partial density 0.400 g/cc
#   Effective density: 5.158 g/cc
#
# We define this as a single material using weight fractions derived from
# the volume fractions and individual densities.
blanket_fill = openmc.Material(name="Blanket fill (SS316LN+H2O)")
eff_density = 0.60 * 7.93 + 0.40 * 1.00  # 5.158 g/cc
ss_wt_frac = (0.60 * 7.93) / eff_density  # weight fraction of SS316LN
h2o_wt_frac = (0.40 * 1.00) / eff_density  # weight fraction of water

# Add SS316LN components scaled by their weight fraction in the mixture
blanket_fill.add_element("Fe", 0.6252 * ss_wt_frac, "wo")
blanket_fill.add_element("Cr", 0.1750 * ss_wt_frac, "wo")
blanket_fill.add_element("Ni", 0.1250 * ss_wt_frac, "wo")
blanket_fill.add_element("Mo", 0.0250 * ss_wt_frac, "wo")
blanket_fill.add_element("Mn", 0.0200 * ss_wt_frac, "wo")
blanket_fill.add_element("Si", 0.0100 * ss_wt_frac, "wo")
blanket_fill.add_element("N",  0.0015 * ss_wt_frac, "wo")
blanket_fill.add_element("C",  0.0003 * ss_wt_frac, "wo")
# Add water components
blanket_fill.add_element("H",  (2 * 1.008 / 18.015) * h2o_wt_frac, "wo")
blanket_fill.add_element("O",  (15.999 / 18.015) * h2o_wt_frac, "wo")
blanket_fill.set_density("g/cm3", eff_density)

# --- 5. Thermal shield (SS304, simplified as SS316LN) ---
# The thermal shield is a stainless steel panel cooled by helium gas,
# positioned between the vacuum vessel and the TF coils. Its function is to
# intercept thermal radiation from the ~100 C vacuum vessel to protect the
# ~4 K superconducting coils. For neutronics purposes, SS304 and SS316LN are
# nearly identical, so we reuse SS316LN.
thermal_shield_mat = ss316ln

# --- 6. TF coil winding pack (homogenized) ---
# The ITER TF coils use Nb3Sn superconductor in a cable-in-conduit conductor
# (CICC) geometry: ~1000 Nb3Sn strands twisted into cables, enclosed in a
# circular SS316LN conduit cooled by supercritical helium at 4.5 K. The
# conduit is embedded in copper stabilizer and glass-epoxy insulation.
#
# A simplified homogenization by volume:
#   SS316LN: 50% at 7.93 g/cc -> 3.965 g/cc
#   Copper:  30% at 8.96 g/cc -> 2.688 g/cc
#   Niobium: 15% at 8.57 g/cc -> 1.286 g/cc
#   Tin:      5% at 7.31 g/cc -> 0.366 g/cc
#   Effective density: 8.305 g/cc
#
# The Nb and Sn represent the Nb3Sn intermetallic compound. In reality,
# the stoichiometry is Nb3Sn (3:1 atomic ratio), but we treat them as
# separate elemental contributions for simplicity.
tf_winding = openmc.Material(name="TF winding pack (homogenized)")
tf_eff_density = 0.50 * 7.93 + 0.30 * 8.96 + 0.15 * 8.57 + 0.05 * 7.31
tf_ss_wf = (0.50 * 7.93) / tf_eff_density
tf_cu_wf = (0.30 * 8.96) / tf_eff_density
tf_nb_wf = (0.15 * 8.57) / tf_eff_density
tf_sn_wf = (0.05 * 7.31) / tf_eff_density

# SS316LN components in the winding pack
tf_winding.add_element("Fe", 0.6252 * tf_ss_wf, "wo")
tf_winding.add_element("Cr", 0.1750 * tf_ss_wf, "wo")
tf_winding.add_element("Ni", 0.1250 * tf_ss_wf, "wo")
tf_winding.add_element("Mo", 0.0250 * tf_ss_wf, "wo")
tf_winding.add_element("Mn", 0.0200 * tf_ss_wf, "wo")
tf_winding.add_element("Si", 0.0100 * tf_ss_wf, "wo")
tf_winding.add_element("N",  0.0015 * tf_ss_wf, "wo")
tf_winding.add_element("C",  0.0003 * tf_ss_wf, "wo")
# Copper in the winding pack
tf_winding.add_element("Cu", tf_cu_wf, "wo")
# Niobium and tin (Nb3Sn superconductor)
tf_winding.add_element("Nb", tf_nb_wf, "wo")
tf_winding.add_element("Sn", tf_sn_wf, "wo")
tf_winding.set_density("g/cm3", tf_eff_density)

# --- Cross section library ---
materials = openmc.Materials([
    beryllium, copper, ss316ln, blanket_fill, tf_winding,
])
materials.cross_sections = "/data/endfb-viii.0-hdf5/cross_sections.xml"
materials.export_to_xml()


# =============================================================================
# Geometry
# =============================================================================
# The geometry is built from concentric ZTorus surfaces centred at the origin
# (x0=0, y0=0, z0=0) with major radius a=R0=620 cm. Each torus has circular
# cross-section (b = c = minor_radius). The space between successive tori
# defines a toroidal shell for each component.
#
# The 40-degree sector is cut by two planes passing through the z-axis:
#   - Plane at phi=0:   the y=0 plane (normal in the y-direction), with the
#                        half-space x > 0 being inside the sector
#   - Plane at phi=40:  a plane making angle 40 degrees with the x-axis
#
# Both sector planes have reflective boundary conditions to simulate the
# full 360-degree torus via symmetry.

# --- Toroidal shell boundaries ---
# Each entry: (name, minor_radius_cm, description)
# The surfaces are listed from innermost (plasma) to outermost (TF coil).
shell_boundaries = [
    ("plasma",           200.0,  "Plasma boundary"),
    ("be_armor",         201.0,  "Be armor outer surface"),
    ("cu_heatsink",      202.0,  "Cu heat sink outer surface"),
    ("fw_structure",     204.0,  "First wall SS316LN outer surface"),
    ("blanket_outer",    249.0,  "Blanket outer surface"),
    ("vv_gap_inner",     253.0,  "VV inner gap / blanket-VV gap"),
    ("vv_inner_shell",   259.0,  "VV inner shell outer surface"),
    ("vv_fill",          289.0,  "VV fill outer surface"),
    ("vv_outer_shell",   295.0,  "VV outer shell outer surface"),
    ("ts_gap",           300.0,  "Thermal shield inner gap"),
    ("thermal_shield",   302.0,  "Thermal shield outer surface"),
    ("tf_gap",           310.0,  "TF coil inner gap"),
    ("tf_case_inner",    315.0,  "TF inner case outer surface"),
    ("tf_winding",       340.0,  "TF winding pack outer surface"),
    ("tf_case_outer",    345.0,  "TF outer case outer surface"),
]

# Create ZTorus surfaces for each boundary
torus_surfaces = {}
for name, minor_r, desc in shell_boundaries:
    torus_surfaces[name] = openmc.ZTorus(
        x0=0.0, y0=0.0, z0=0.0,
        a=R0, b=minor_r, c=minor_r,
        name=desc,
    )

# --- Sector boundary planes ---
# Plane at phi = 0 degrees: the y = 0 plane
# We want the sector to span from phi=0 to phi=40 degrees, where phi is
# measured from the +x axis in the x-y plane. The y=0 plane divides space
# into y>0 (phi in [0,180]) and y<0 (phi in [180,360]). For our sector,
# we need y >= 0 on this boundary.
phi_rad = math.radians(sector_angle)

plane_phi0 = openmc.YPlane(
    y0=0.0,
    boundary_type="reflective",
    name="Sector plane at phi=0 (y=0)",
)

# Plane at phi = 40 degrees:
# A plane through the z-axis at angle phi from the x-axis has the equation:
#   x * sin(phi) - y * cos(phi) = 0
# Points inside the sector (phi < 40 deg) satisfy x*sin(40)-y*cos(40) > 0,
# i.e. they are on the positive side of this plane.
plane_phi40 = openmc.Plane(
    a=math.sin(phi_rad),
    b=-math.cos(phi_rad),
    c=0.0,
    d=0.0,
    boundary_type="reflective",
    name=f"Sector plane at phi={sector_angle} deg",
)

# --- Outer bounding sphere ---
# The outermost TF coil surface has minor radius 345 cm on a major radius
# of 620 cm. The farthest point from the origin is at R0 + 345 = 965 cm
# in the midplane, and z = +/- 345 cm above/below. A bounding sphere of
# radius 1000 cm comfortably encloses the entire geometry.
boundary_sphere = openmc.Sphere(
    r=1000.0,
    boundary_type="vacuum",
    name="Outer vacuum boundary",
)

# --- Sector region (common to all cells) ---
# The sector region selects the angular wedge between phi=0 and phi=40.
# Points in the sector have y >= 0 (above the phi=0 plane) AND satisfy
# x*sin(40) - y*cos(40) >= 0 (below the phi=40 plane).
sector_region = +plane_phi0 & +plane_phi40

# --- Build cells for each toroidal shell ---
cells = []

# Helper to get torus surfaces by index
torus_list = list(torus_surfaces.values())
shell_names = [name for name, _, _ in shell_boundaries]

# 1. Plasma (void): inside the innermost torus
plasma_cell = openmc.Cell(
    name="Plasma (void)",
    region=-torus_list[0] & sector_region & -boundary_sphere,
)
cells.append(plasma_cell)

# 2-14. Toroidal shells between successive torus surfaces
# Each shell is the region between two consecutive tori, within the sector.
shell_materials = {
    # (inner_surface_name, outer_surface_name): (material, cell_name)
    ("plasma",         "be_armor"):        (beryllium,    "Be armor (first wall)"),
    ("be_armor",       "cu_heatsink"):     (copper,       "Cu heat sink"),
    ("cu_heatsink",    "fw_structure"):    (ss316ln,      "First wall SS316LN structure"),
    ("fw_structure",   "blanket_outer"):   (blanket_fill, "Shielding blanket (SS316LN+H2O)"),
    ("blanket_outer",  "vv_gap_inner"):    (None,         "Blanket-VV gap (void)"),
    ("vv_gap_inner",   "vv_inner_shell"):  (ss316ln,      "VV inner shell (SS316LN)"),
    ("vv_inner_shell", "vv_fill"):         (blanket_fill, "VV fill (SS316LN+H2O)"),
    ("vv_fill",        "vv_outer_shell"):  (ss316ln,      "VV outer shell (SS316LN)"),
    ("vv_outer_shell", "ts_gap"):          (None,         "VV-thermal shield gap (void)"),
    ("ts_gap",         "thermal_shield"):  (ss316ln,      "Thermal shield (SS304/SS316LN)"),
    ("thermal_shield", "tf_gap"):          (None,         "Thermal shield-TF coil gap (void)"),
    ("tf_gap",         "tf_case_inner"):   (ss316ln,      "TF coil inner case (SS316LN)"),
    ("tf_case_inner",  "tf_winding"):      (tf_winding,   "TF winding pack (homogenized)"),
    ("tf_winding",     "tf_case_outer"):   (ss316ln,      "TF coil outer case (SS316LN)"),
}

# Store cell references for tallies
cell_by_name = {}

for (inner_name, outer_name), (mat, cell_name) in shell_materials.items():
    inner_surf = torus_surfaces[inner_name]
    outer_surf = torus_surfaces[outer_name]

    shell_cell = openmc.Cell(
        name=cell_name,
        fill=mat,
        region=+inner_surf & -outer_surf & sector_region & -boundary_sphere,
    )
    cells.append(shell_cell)
    cell_by_name[cell_name] = shell_cell

# 15. Void region: outside the outermost TF coil but inside the bounding sphere
void_cell = openmc.Cell(
    name="External void",
    region=+torus_list[-1] & sector_region & -boundary_sphere,
)
cells.append(void_cell)

# --- Build universe and geometry ---
root_universe = openmc.Universe(name="Root universe", cells=cells)
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()


# =============================================================================
# Source definition
# =============================================================================
# The D-T fusion neutron source is distributed within the plasma volume.
# In a real tokamak, the fusion reaction rate peaks at the magnetic axis
# and decreases roughly as (1 - (r/a)^2)^2 where r is the distance from
# the magnetic axis and a is the plasma minor radius. For this simplified
# model, we use a uniform distribution in R and Z within reduced plasma
# bounds (80% of the minor radius) to avoid sampling near the separatrix
# where the reaction rate is negligible.
#
# Cylindrical coordinates (R, phi, Z):
#   R: 460 to 780 cm (R0 +/- 160 cm = R0 +/- 0.8*a)
#   phi: 0 to 40 degrees
#   Z: -160 to +160 cm (+/- 0.8*a)

source = openmc.IndependentSource()

# Radial distribution: uniform between R0-160 and R0+160 cm
r_dist = openmc.stats.Uniform(a=R0 - 0.8 * a_plasma, b=R0 + 0.8 * a_plasma)

# Azimuthal distribution: uniform across the 40-degree sector
phi_dist = openmc.stats.Uniform(a=0.0, b=phi_rad)

# Vertical distribution: uniform between -160 and +160 cm
z_dist = openmc.stats.Uniform(a=-0.8 * a_plasma, b=0.8 * a_plasma)

source.space = openmc.stats.CylindricalIndependent(
    r=r_dist, phi=phi_dist, z=z_dist,
)

# D-T fusion neutrons: monoenergetic at 14.1 MeV
source.energy = openmc.stats.Discrete([14.1e6], [1.0])

# Isotropic emission (fusion neutrons are emitted isotropically in the
# centre-of-mass frame; the small correction from beam-target kinematics
# is negligible for this benchmark)
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
settings.photon_transport = False
settings.temperature = {"default": 400.0}
settings.output = {"tallies": True}

# --- Weight windows ---
if args.use_ww:
    wws = openmc.hdf5_to_wws(args.use_ww)
    settings.weight_windows = wws

settings.export_to_xml()


# =============================================================================
# Tallies
# =============================================================================
# Four categories of tallies for ITER-relevant quantities:
#
#   1. Nuclear heating in blanket, VV, and TF coil regions
#   2. Fast neutron flux (E > 0.1 MeV) at the TF winding pack
#   3. Neutron flux spectra at multiple radial depths
#   4. Neutron current (wall loading) on the first wall surface

tallies = openmc.Tallies()

# --- Energy bin structures ---
# Coarse energy groups for spectral tallies: 200 logarithmic bins
energy_bins = np.logspace(
    np.log10(1.0e-3),   # 1 meV (thermal)
    np.log10(15.0e6),    # 15 MeV (above D-T peak)
    201,
)
energy_filter = openmc.EnergyFilter(energy_bins)

# Fast neutron threshold filter (E > 0.1 MeV = 1e5 eV)
fast_energy_filter = openmc.EnergyFilter([1.0e5, 15.0e6])

if args.small_tallies:
    # ----- Reduced tally set -----

    # Total nuclear heating in all material cells
    heating_cells = [
        cell_by_name[n] for n in cell_by_name
        if "void" not in n.lower() and "gap" not in n.lower()
    ]
    heating_filter = openmc.CellFilter(heating_cells)
    tally_heat = openmc.Tally(name="total_heating")
    tally_heat.filters = [heating_filter]
    tally_heat.scores = ["heating"]
    tallies.append(tally_heat)

    # Fast flux at TF winding pack
    tf_wp_cell = cell_by_name["TF winding pack (homogenized)"]
    tf_filter = openmc.CellFilter([tf_wp_cell])
    tally_tf_fast = openmc.Tally(name="tf_coil_fast_flux")
    tally_tf_fast.filters = [tf_filter, fast_energy_filter]
    tally_tf_fast.scores = ["flux"]
    tallies.append(tally_tf_fast)

else:
    # ----- Full tally set -----

    # 1. Nuclear heating by component
    # Heating is tallied in each major component to assess power deposition.
    # The dominant energy deposition is in the blanket (~85% of neutron
    # energy), with residual heating in the VV and TF coils.
    heating_components = [
        ("Be armor (first wall)",              "heating_be_armor"),
        ("Cu heat sink",                       "heating_cu_heatsink"),
        ("First wall SS316LN structure",       "heating_fw_structure"),
        ("Shielding blanket (SS316LN+H2O)",    "heating_blanket"),
        ("VV inner shell (SS316LN)",           "heating_vv_inner"),
        ("VV fill (SS316LN+H2O)",              "heating_vv_fill"),
        ("VV outer shell (SS316LN)",           "heating_vv_outer"),
        ("Thermal shield (SS304/SS316LN)",     "heating_thermal_shield"),
        ("TF coil inner case (SS316LN)",       "heating_tf_case_inner"),
        ("TF winding pack (homogenized)",      "heating_tf_winding"),
        ("TF coil outer case (SS316LN)",       "heating_tf_case_outer"),
    ]

    for cell_name, tally_name in heating_components:
        cell = cell_by_name[cell_name]
        cf = openmc.CellFilter([cell])
        tally = openmc.Tally(name=tally_name)
        tally.filters = [cf]
        tally.scores = ["heating"]
        tallies.append(tally)

    # 2. Fast neutron flux (E > 0.1 MeV) at TF winding pack
    # The fast neutron fluence at the TF coils is a critical design parameter.
    # ITER requires the cumulative fast fluence to remain below ~10^22 n/m^2
    # over the device lifetime to limit radiation damage to the Nb3Sn
    # superconductor (which degrades critical current density Jc under
    # irradiation) and to the epoxy insulation.
    tf_wp_cell = cell_by_name["TF winding pack (homogenized)"]
    tf_filter = openmc.CellFilter([tf_wp_cell])

    tally_tf_fast = openmc.Tally(name="tf_coil_fast_flux")
    tally_tf_fast.filters = [tf_filter, fast_energy_filter]
    tally_tf_fast.scores = ["flux"]
    tallies.append(tally_tf_fast)

    # 3. Neutron flux spectra at key radial locations
    # Spectra at four depths through the shielding provide detailed
    # validation data and insight into the neutron moderation process.
    spectrum_locations = [
        ("Be armor (first wall)",              "spectrum_first_wall"),
        ("Shielding blanket (SS316LN+H2O)",    "spectrum_blanket"),
        ("VV fill (SS316LN+H2O)",              "spectrum_vv"),
        ("TF winding pack (homogenized)",      "spectrum_tf_coil"),
    ]

    for cell_name, tally_name in spectrum_locations:
        cell = cell_by_name[cell_name]
        cf = openmc.CellFilter([cell])
        tally = openmc.Tally(name=tally_name)
        tally.filters = [cf, energy_filter]
        tally.scores = ["flux"]
        tallies.append(tally)

    # 4. Total neutron flux at first wall (proxy for wall loading)
    # The neutron wall loading (NWL) is the neutron power per unit area
    # crossing the first wall surface. In ITER, the average NWL is about
    # 0.56 MW/m^2, peaking at ~1.0 MW/m^2 at the outboard midplane.
    # We tally the total flux in the Be armor cell as a proxy.
    fw_cell = cell_by_name["Be armor (first wall)"]
    fw_filter = openmc.CellFilter([fw_cell])
    tally_nwl = openmc.Tally(name="first_wall_flux")
    tally_nwl.filters = [fw_filter]
    tally_nwl.scores = ["flux", "current"]
    tallies.append(tally_nwl)

tallies.export_to_xml()


# =============================================================================
# Weight window generation (optional)
# =============================================================================
def generate_weight_windows(model):
    """Generate weight windows using FW-CADIS via the random ray solver.

    This function creates a deep copy of the model, converts it to a
    multigroup random ray calculation, runs it to generate adjoint-informed
    weight windows, and returns the path to the weight window file.

    The FW-CADIS (Forward-Weighted Consistent Adjoint Driven Importance
    Sampling) method uses an adjoint flux solution to define importance
    maps that guide particle transport toward the tally regions of interest.
    This is essential for deep-penetration shielding problems like ITER,
    where the flux attenuates by 8+ orders of magnitude between the plasma
    and the TF coils.

    Parameters
    ----------
    model : openmc.Model
        The base model to generate weight windows for.

    Returns
    -------
    str
        Path to the generated weight_windows.h5 file.
    """
    model_rr = copy.deepcopy(model)

    # Convert continuous-energy model to multigroup using CASMO-4 group
    # structure (70 energy groups). The material_wise method generates
    # multigroup cross sections for each material independently.
    model_rr.convert_to_multigroup(
        method="material_wise",
        groups="CASMO-4",
        nparticles=3000,
    )

    # Convert to random ray solver, which determines appropriate settings
    # for the deterministic flux solution
    model_rr.convert_to_random_ray()

    # Create a regular mesh covering the bounding box for the weight windows
    bbox = model_rr.geometry.bounding_box
    lower = bbox.lower_left
    upper = bbox.upper_right

    # Mesh spacing of ~10 cm per element
    n_x = max(1, int((upper[0] - lower[0]) / 10.0))
    n_y = max(1, int((upper[1] - lower[1]) / 10.0))
    n_z = max(1, int((upper[2] - lower[2]) / 10.0))

    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (n_x, n_y, n_z)
    ww_mesh.lower_left = list(lower)
    ww_mesh.upper_right = list(upper)

    # Configure the random ray source regions on the mesh
    root = model_rr.geometry.root_universe
    model_rr.settings.random_ray['source_region_meshes'] = [(ww_mesh, [root])]

    # Set up the weight window generator
    wwg = openmc.WeightWindowGenerator(
        method="fw_cadis",
        mesh=ww_mesh,
        max_realizations=model_rr.settings.batches,
    )
    model_rr.settings.weight_window_generators = wwg

    # Export and run
    model_rr.export_to_xml()
    openmc.run()

    return "weight_windows.h5"


# =============================================================================
# Build Model object and optionally generate weight windows
# =============================================================================
def build_model():
    """Build and return the ITER sector OpenMC model.

    Returns
    -------
    openmc.Model
        The complete model with geometry, materials, settings, and tallies.
    """
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies,
    )
    return model


if args.generate_ww:
    print("Generating weight windows via FW-CADIS random ray...")
    model = build_model()
    ww_file = generate_weight_windows(model)
    print(f"Weight windows written to: {ww_file}")
    print("Re-run with: python model.py --use-ww weight_windows.h5")
else:
    model = build_model()
    model.export_to_xml()

    print("=" * 72)
    print("Simplified ITER 40-Degree Sector -- OpenMC Model")
    print("=" * 72)
    print(f"  Major radius R0:   {R0} cm")
    print(f"  Plasma minor rad:  {a_plasma} cm")
    print(f"  Sector angle:      {sector_angle} deg (1/{int(360/sector_angle)} torus)")
    print(f"  Radial build:      {len(shell_boundaries)} surfaces, "
          f"outermost at {shell_boundaries[-1][1]} cm minor radius")
    print(f"  Materials:         {len(materials)} defined")
    print(f"  Source:            14.1 MeV D-T, uniform in plasma volume")
    print(f"  Particles:         {args.particles:,}/batch x {args.batches} batches"
          f" = {args.particles * args.batches:,} total")
    if args.use_ww:
        print(f"  Weight windows:    {args.use_ww}")
    else:
        print("  Weight windows:    NONE (consider --generate-ww for deep penetration)")
    print(f"  Small tallies:     {args.small_tallies}")
    print("=" * 72)
    print("XML files written: materials.xml, geometry.xml, settings.xml, tallies.xml")
    print("Run with: openmc")
