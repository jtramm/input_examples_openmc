#!/usr/bin/env python3
"""
HTR-10 Pebble Bed High-Temperature Gas-Cooled Reactor Benchmark
================================================================

This script builds an OpenMC model for the Chinese HTR-10, a 10 MWth pebble
bed high-temperature gas-cooled reactor (HTGR) that achieved first criticality
at Tsinghua University's Institute of Nuclear and New Energy Technology
(INET) in December 2000.

Background on Pebble Bed Reactors
----------------------------------
Pebble bed reactors are a type of Very High Temperature Reactor (VHTR) in
which the fuel takes the form of tennis-ball-sized graphite spheres ("pebbles")
containing thousands of tiny TRISO (TRi-structural ISOtropic) coated fuel
particles. The pebbles are randomly packed in a cylindrical core cavity and
cooled by helium gas flowing through the interstitial spaces.

Key advantages of pebble bed design:
  - Inherent safety: strong negative temperature coefficient; the fuel can
    withstand loss of coolant without melting (demonstrated in HTR-10)
  - Online refueling: pebbles are continuously added at the top and removed
    at the bottom, enabling high capacity factors
  - High temperature operation (~700-950 C coolant outlet) enabling
    efficient electricity generation and process heat applications

TRISO Fuel Particles: The Double Heterogeneity Problem
-------------------------------------------------------
TRISO particles are the fundamental fuel unit. Each particle consists of a
fissile kernel (UO2) surrounded by four concentric coating layers:

  1. UO2 kernel (r = 0.025 cm) -- the fissile material, 17% enriched
  2. Porous pyrolytic carbon buffer (outer r = 0.034 cm) -- provides void
     space for fission gas release and attenuates fission recoils
  3. Inner dense pyrolytic carbon, IPyC (outer r = 0.038 cm) -- acts as a
     diffusion barrier and provides structural support for the SiC layer
  4. Silicon carbide, SiC (outer r = 0.0415 cm) -- the primary fission
     product barrier and pressure vessel; retains fission products up to
     ~1600 C
  5. Outer dense pyrolytic carbon, OPyC (outer r = 0.0455 cm) -- provides
     additional fission product barrier and protects SiC during handling

About 8,335 TRISO particles are dispersed in a graphite matrix within the
inner "fuel zone" of each pebble (r = 2.5 cm), and the fuel zone is
surrounded by a 0.5 cm thick graphite shell (total pebble r = 3.0 cm).

This creates a "double heterogeneity" in the neutronics:
  - First level: the random distribution of TRISO micro-particles within
    the graphite matrix of the fuel zone (self-shielding of individual
    TRISO particles affects resonance absorption)
  - Second level: the random packing of fuel pebbles within the core
    cavity (pebble-scale flux gradients and inter-pebble streaming)

Modeling these two levels of heterogeneity correctly is critical for
accurate k-eff predictions. OpenMC can handle both levels explicitly using
its random sphere packing capabilities (openmc.model.pack_spheres and
openmc.model.TRISO).

HTR-10 Specifications
---------------------
  Core diameter:          180 cm (radius 90 cm)
  Critical height:        ~123-126 cm
  Fuel pebble radius:     3.0 cm (outer)
  Fuel zone radius:       2.5 cm (contains TRISO particles in graphite)
  Graphite shell:         0.5 cm thick
  Pebble packing fraction: ~0.61
  TRISO particles/pebble: 8,335
  UO2 enrichment:         17% U-235
  UO2 kernel density:     10.4 g/cc
  Coolant:                Helium at ~3 MPa

Model Modes
-----------
This script supports two model configurations:

  --model pin   (default)
    A single fuel pebble in an infinite lattice approximation. The pebble
    is placed at the center of a cube with reflective boundary conditions.
    TRISO particles are explicitly modeled inside the fuel zone using
    OpenMC's pack_spheres() and TRISO lattice. This captures the first
    level of heterogeneity (TRISO self-shielding) accurately.

  --model core
    A simplified cylindrical core model. Because explicitly packing ~27,000
    pebbles each containing 8,335 explicit TRISO particles would require
    enormous memory (~225 million TRISO particles), this mode uses
    homogenized fuel pebbles: the TRISO fuel zone material is volume-
    averaged over all TRISO layers and the graphite matrix. The pebbles
    (with homogenized fuel zone + graphite shell) are then packed into
    a cylindrical core cavity surrounded by a graphite reflector.
    This captures the second level of heterogeneity (pebble packing)
    but not the first (TRISO self-shielding).

References
----------
  [1] Wu, Z., et al., "The design features of the HTR-10," Nuclear
      Engineering and Design, 218 (2002) 25-32.
  [2] IAEA-TECDOC-1382, "Evaluation of high temperature gas cooled reactor
      performance: Benchmark analysis related to initial testing of the
      HTTR and HTR-10," IAEA, Vienna (2003).
  [3] Terry, W.K., et al., "Evaluation of the Initial Critical
      Configuration of the HTR-10 Pebble-Bed Reactor," Report
      HTR10-GCR-RESR-001, IRPhEP (2007).

Usage:
    python model.py [--particles N] [--batches N] [--inactive N]
                    [--model {pin,core}] [--run]
"""

import argparse
from math import pi

import numpy as np
import openmc
import openmc.model


# =============================================================================
# Physical Constants and Geometry Parameters
# =============================================================================

# --- TRISO particle layer radii (cm) ---
# These are defined from the center outward. Each successive radius is the
# outer radius of that layer.
KERNEL_RADIUS = 0.025       # UO2 kernel
BUFFER_OUTER_R = 0.034      # Porous PyC buffer
IPYC_OUTER_R = 0.038        # Inner dense PyC
SIC_OUTER_R = 0.0415        # Silicon carbide
OPYC_OUTER_R = 0.0455       # Outer dense PyC (= TRISO outer radius)

# --- Pebble dimensions (cm) ---
FUEL_ZONE_RADIUS = 2.5      # Radius of inner fuel zone containing TRISO
PEBBLE_OUTER_RADIUS = 3.0   # Outer radius including graphite shell
GRAPHITE_SHELL_THICKNESS = PEBBLE_OUTER_RADIUS - FUEL_ZONE_RADIUS  # 0.5 cm

# --- Core dimensions (cm) ---
CORE_RADIUS = 90.0           # Core cavity radius
CORE_HEIGHT = 126.0          # Critical height (approximate)
REFLECTOR_THICKNESS = 100.0  # Graphite reflector thickness (radial and axial)

# --- Packing ---
PEBBLE_PACKING_FRACTION = 0.61   # Pebble packing fraction in core
TRISO_PER_PEBBLE = 8335          # Number of TRISO particles per pebble

# --- Material densities (g/cc) ---
UO2_DENSITY = 10.4
BUFFER_DENSITY = 1.1
IPYC_DENSITY = 1.9
SIC_DENSITY = 3.18
OPYC_DENSITY = 1.9
MATRIX_GRAPHITE_DENSITY = 1.73   # Graphite matrix inside fuel zone
SHELL_GRAPHITE_DENSITY = 1.73    # Pebble outer shell graphite
REFLECTOR_GRAPHITE_DENSITY = 1.76
HELIUM_DENSITY = 0.005           # He coolant at operating pressure


# =============================================================================
# Material Definitions
# =============================================================================

def create_materials():
    """
    Define all materials for the HTR-10 model.

    Returns
    -------
    dict
        Dictionary mapping material name to openmc.Material object.
    """
    mats = {}

    # --- UO2 fuel kernel ---
    # 17% enriched in U-235. UO2 stoichiometry: 1 U + 2 O atoms.
    fuel = openmc.Material(name='UO2 Kernel')
    fuel.set_density('g/cm3', UO2_DENSITY)
    fuel.add_nuclide('U235', 0.17, 'wo')   # 17 wt% U-235
    fuel.add_nuclide('U238', 0.83, 'wo')   # 83 wt% U-238 (balance of uranium)
    # Add oxygen in UO2 stoichiometric ratio. For weight fractions of U,
    # we add O separately. The mass fraction of O in UO2 is:
    #   2 * 16.0 / (238.0 + 2*16.0) = 0.1184
    # But since we specified U by weight, we need to use 'ao' or adjust.
    # Simpler: specify by atom fractions.
    fuel = openmc.Material(name='UO2 Kernel')
    fuel.set_density('g/cm3', UO2_DENSITY)
    fuel.add_element('U', 1.0, enrichment=17.0)  # enrichment in wt%
    fuel.add_element('O', 2.0)
    mats['fuel'] = fuel

    # --- Porous carbon buffer ---
    # Low-density pyrolytic carbon that provides void volume for fission
    # gas release and attenuates fission product recoils.
    buffer_pyc = openmc.Material(name='Buffer PyC')
    buffer_pyc.set_density('g/cm3', BUFFER_DENSITY)
    buffer_pyc.add_element('C', 1.0)
    buffer_pyc.add_s_alpha_beta('c_Graphite')
    mats['buffer'] = buffer_pyc

    # --- Inner dense pyrolytic carbon (IPyC) ---
    # Structural layer that also serves as diffusion barrier.
    ipyc = openmc.Material(name='Inner PyC')
    ipyc.set_density('g/cm3', IPYC_DENSITY)
    ipyc.add_element('C', 1.0)
    ipyc.add_s_alpha_beta('c_Graphite')
    mats['ipyc'] = ipyc

    # --- Silicon carbide ---
    # Primary pressure boundary and fission product retention layer.
    # SiC retains fission products up to ~1600 C.
    sic = openmc.Material(name='SiC')
    sic.set_density('g/cm3', SIC_DENSITY)
    sic.add_element('Si', 1.0)
    sic.add_element('C', 1.0)
    mats['sic'] = sic

    # --- Outer dense pyrolytic carbon (OPyC) ---
    # Additional fission product barrier and mechanical protection for SiC.
    opyc = openmc.Material(name='Outer PyC')
    opyc.set_density('g/cm3', OPYC_DENSITY)
    opyc.add_element('C', 1.0)
    opyc.add_s_alpha_beta('c_Graphite')
    mats['opyc'] = opyc

    # --- Graphite matrix (inside fuel zone of pebble) ---
    # The graphite in which TRISO particles are dispersed within the fuel zone.
    matrix = openmc.Material(name='Graphite Matrix')
    matrix.set_density('g/cm3', MATRIX_GRAPHITE_DENSITY)
    matrix.add_element('C', 1.0)
    matrix.add_s_alpha_beta('c_Graphite')
    mats['matrix'] = matrix

    # --- Graphite shell (outer pebble shell) ---
    # Pure graphite shell surrounding the fuel zone of each pebble.
    shell = openmc.Material(name='Graphite Shell')
    shell.set_density('g/cm3', SHELL_GRAPHITE_DENSITY)
    shell.add_element('C', 1.0)
    shell.add_s_alpha_beta('c_Graphite')
    mats['shell'] = shell

    # --- Graphite reflector ---
    # Surrounds the pebble bed core; slightly higher density than pebble graphite.
    reflector = openmc.Material(name='Graphite Reflector')
    reflector.set_density('g/cm3', REFLECTOR_GRAPHITE_DENSITY)
    reflector.add_element('C', 1.0)
    reflector.add_s_alpha_beta('c_Graphite')
    mats['reflector'] = reflector

    # --- Helium coolant ---
    # Flows through interstitial spaces between pebbles. At ~3 MPa operating
    # pressure and ~250-700 C, density is approximately 0.005 g/cc.
    helium = openmc.Material(name='Helium Coolant')
    helium.set_density('g/cm3', HELIUM_DENSITY)
    helium.add_element('He', 1.0)
    mats['helium'] = helium

    return mats


# =============================================================================
# Homogenized Fuel Zone Material
# =============================================================================

def create_homogenized_fuel_zone(mats):
    """
    Create a volume-averaged homogenized material for the fuel zone of a pebble.

    In the "core" model mode, we cannot afford to model 8,335 explicit TRISO
    particles in each of ~27,000 pebbles. Instead, we homogenize the fuel zone
    by computing the volume fraction of each TRISO layer and the graphite
    matrix, then creating a mixture material.

    The fuel zone is a sphere of radius 2.5 cm. The TRISO particles occupy
    a certain volume fraction within it, and the remainder is graphite matrix.

    Volume fractions within a single TRISO particle:
      kernel:  V = (4/3)*pi*r_kernel^3
      buffer:  V = (4/3)*pi*(r_buffer^3 - r_kernel^3)
      ipyc:    V = (4/3)*pi*(r_ipyc^3 - r_buffer^3)
      sic:     V = (4/3)*pi*(r_sic^3 - r_ipyc^3)
      opyc:    V = (4/3)*pi*(r_opyc^3 - r_sic^3)

    The total TRISO volume fraction in the fuel zone is:
      f_triso = N_triso * V_triso / V_fuelzone
    where N_triso = 8335 and V_fuelzone = (4/3)*pi*2.5^3.

    Returns
    -------
    openmc.Material
        Homogenized fuel zone material.
    """
    # Volumes of each TRISO layer (cm^3)
    v_kernel = (4.0 / 3.0) * pi * KERNEL_RADIUS**3
    v_buffer = (4.0 / 3.0) * pi * (BUFFER_OUTER_R**3 - KERNEL_RADIUS**3)
    v_ipyc = (4.0 / 3.0) * pi * (IPYC_OUTER_R**3 - BUFFER_OUTER_R**3)
    v_sic = (4.0 / 3.0) * pi * (SIC_OUTER_R**3 - IPYC_OUTER_R**3)
    v_opyc = (4.0 / 3.0) * pi * (OPYC_OUTER_R**3 - SIC_OUTER_R**3)
    v_triso = (4.0 / 3.0) * pi * OPYC_OUTER_R**3  # total volume of one TRISO

    # Volume of the fuel zone sphere
    v_fuelzone = (4.0 / 3.0) * pi * FUEL_ZONE_RADIUS**3

    # Total TRISO packing fraction within the fuel zone
    f_triso = TRISO_PER_PEBBLE * v_triso / v_fuelzone
    f_matrix = 1.0 - f_triso  # graphite matrix fills the rest

    # Volume fraction of each TRISO layer within the fuel zone
    f_kernel = TRISO_PER_PEBBLE * v_kernel / v_fuelzone
    f_buffer = TRISO_PER_PEBBLE * v_buffer / v_fuelzone
    f_ipyc = TRISO_PER_PEBBLE * v_ipyc / v_fuelzone
    f_sic = TRISO_PER_PEBBLE * v_sic / v_fuelzone
    f_opyc = TRISO_PER_PEBBLE * v_opyc / v_fuelzone

    print(f"  TRISO packing fraction in fuel zone: {f_triso:.4f}")
    print(f"  Graphite matrix fraction:            {f_matrix:.4f}")
    print(f"  Volume fractions - kernel: {f_kernel:.4f}, buffer: {f_buffer:.4f}, "
          f"IPyC: {f_ipyc:.4f}, SiC: {f_sic:.4f}, OPyC: {f_opyc:.4f}")

    # Compute homogenized density using volume-weighted average
    # rho_hom = sum(f_i * rho_i)
    rho_hom = (f_kernel * UO2_DENSITY
               + f_buffer * BUFFER_DENSITY
               + f_ipyc * IPYC_DENSITY
               + f_sic * SIC_DENSITY
               + f_opyc * OPYC_DENSITY
               + f_matrix * MATRIX_GRAPHITE_DENSITY)
    print(f"  Homogenized fuel zone density: {rho_hom:.4f} g/cc")

    # Create homogenized material as a mixture
    # We use openmc.Material.mix_materials() for a proper volume-weighted mix.
    homo_fuel = openmc.Material.mix_materials(
        [mats['fuel'], mats['buffer'], mats['ipyc'], mats['sic'],
         mats['opyc'], mats['matrix']],
        [f_kernel, f_buffer, f_ipyc, f_sic, f_opyc, f_matrix],
        'vo'  # volume fractions
    )
    homo_fuel.name = 'Homogenized Fuel Zone'

    return homo_fuel


# =============================================================================
# Pin Model: Single Pebble with Explicit TRISO
# =============================================================================

def build_pin_model(mats, particles, batches, inactive):
    """
    Build a single fuel pebble model with explicit TRISO particles.

    This model places one fuel pebble at the center of a cube with reflective
    boundary conditions on all faces, simulating an infinite lattice of
    identical pebbles. The TRISO particles inside the fuel zone are modeled
    explicitly using OpenMC's random sphere packing.

    This captures the "first heterogeneity" -- the self-shielding effect of
    individual TRISO particles embedded in the graphite matrix. The resonance
    absorption in the UO2 kernels is strongly affected by whether the fuel
    is modeled as discrete particles or as a homogeneous mixture.

    Parameters
    ----------
    mats : dict
        Material dictionary from create_materials().
    particles : int
        Neutron histories per batch.
    batches : int
        Total number of batches.
    inactive : int
        Number of inactive batches.

    Returns
    -------
    openmc.Model
        Complete OpenMC model.
    """
    model = openmc.Model()

    # -------------------------------------------------------------------------
    # TRISO universe: define the five concentric layers
    # -------------------------------------------------------------------------
    # The TRISO particle is modeled as a set of concentric spherical shells.
    # Each shell is a cell filled with the corresponding material.
    # The outermost cell extends to infinity (it will be clipped by the
    # TRISO boundary sphere).
    s_kernel = openmc.Sphere(r=KERNEL_RADIUS)
    s_buffer = openmc.Sphere(r=BUFFER_OUTER_R)
    s_ipyc = openmc.Sphere(r=IPYC_OUTER_R)
    s_sic = openmc.Sphere(r=SIC_OUTER_R)

    c_kernel = openmc.Cell(fill=mats['fuel'], region=-s_kernel)
    c_buffer = openmc.Cell(fill=mats['buffer'], region=+s_kernel & -s_buffer)
    c_ipyc = openmc.Cell(fill=mats['ipyc'], region=+s_buffer & -s_ipyc)
    c_sic = openmc.Cell(fill=mats['sic'], region=+s_ipyc & -s_sic)
    c_opyc = openmc.Cell(fill=mats['opyc'], region=+s_sic)

    triso_univ = openmc.Universe(cells=[c_kernel, c_buffer, c_ipyc, c_sic, c_opyc])

    # -------------------------------------------------------------------------
    # Pack TRISO particles into the fuel zone sphere
    # -------------------------------------------------------------------------
    # The fuel zone is a sphere of radius 2.5 cm centered at the origin.
    # We pack TRISO particles (outer radius 0.0455 cm) into this sphere.
    #
    # With 8,335 particles, the TRISO packing fraction would be:
    #   f = 8335 * (4/3)*pi*0.0455^3 / ((4/3)*pi*2.5^3) = ~0.048
    # This is well below the RSP limit of 0.38, so packing should be fast.
    #
    # However, packing 8,335 particles and creating a TRISO lattice for all
    # of them can still be memory-intensive. We use a smaller number for
    # tractability while preserving the correct packing fraction.
    fuel_zone_sphere = openmc.Sphere(r=FUEL_ZONE_RADIUS)
    fuel_zone_region = -fuel_zone_sphere

    # Pack spheres using the correct number of TRISO particles.
    # With ~8335 particles at r=0.0455 cm in a sphere of r=2.5 cm,
    # packing fraction is only ~0.048, which is very manageable.
    print("Packing TRISO particles in fuel zone...")
    print(f"  Number of TRISO particles: {TRISO_PER_PEBBLE}")
    centers = openmc.model.pack_spheres(
        radius=OPYC_OUTER_R,
        region=fuel_zone_region,
        num_spheres=TRISO_PER_PEBBLE,
        seed=42
    )
    print(f"  Packed {len(centers)} TRISO particles")

    # Create TRISO cell objects from the packed centers
    trisos = [openmc.model.TRISO(OPYC_OUTER_R, triso_univ, c)
              for c in centers]

    # -------------------------------------------------------------------------
    # Create TRISO lattice within the fuel zone
    # -------------------------------------------------------------------------
    # OpenMC's create_triso_lattice() overlays a rectilinear lattice on the
    # fuel zone and assigns each TRISO particle to the lattice cell(s) it
    # overlaps. This dramatically speeds up particle tracking compared to
    # having all ~8000 TRISO surfaces in a single cell.
    #
    # The lattice shape determines the granularity. More cells = faster
    # tracking but more memory. A shape of (5,5,5) to (10,10,10) is typical.
    ll, ur = fuel_zone_region.bounding_box
    shape = (8, 8, 8)
    pitch = (ur - ll) / shape

    print("Creating TRISO lattice...")
    triso_lattice = openmc.model.create_triso_lattice(
        trisos, ll, pitch, shape, mats['matrix']
    )

    # -------------------------------------------------------------------------
    # Pebble universe: fuel zone (TRISO lattice) + graphite shell
    # -------------------------------------------------------------------------
    pebble_sphere = openmc.Sphere(r=PEBBLE_OUTER_RADIUS)

    # Fuel zone cell: filled with the TRISO lattice
    fuel_zone_cell = openmc.Cell(
        name='Fuel Zone',
        fill=triso_lattice,
        region=-fuel_zone_sphere
    )
    # Graphite shell: annular region between fuel zone and pebble surface
    shell_cell = openmc.Cell(
        name='Graphite Shell',
        fill=mats['shell'],
        region=+fuel_zone_sphere & -pebble_sphere
    )

    pebble_univ = openmc.Universe(cells=[fuel_zone_cell, shell_cell])

    # -------------------------------------------------------------------------
    # Bounding box with reflective BCs (infinite lattice of pebbles)
    # -------------------------------------------------------------------------
    # We place the pebble in a cube. The cube side length is chosen to give
    # the correct pebble packing fraction. For a single sphere in a cube:
    #   PF = V_sphere / V_cube = (4/3)*pi*R^3 / L^3
    #   L = ((4/3)*pi*R^3 / PF)^(1/3)
    v_pebble = (4.0 / 3.0) * pi * PEBBLE_OUTER_RADIUS**3
    cube_side = (v_pebble / PEBBLE_PACKING_FRACTION) ** (1.0 / 3.0)
    half = cube_side / 2.0

    print(f"  Cube side for PF={PEBBLE_PACKING_FRACTION}: {cube_side:.4f} cm")
    print(f"  Half-side: {half:.4f} cm")

    # Reflective boundary planes
    xmin = openmc.XPlane(-half, boundary_type='reflective')
    xmax = openmc.XPlane(+half, boundary_type='reflective')
    ymin = openmc.YPlane(-half, boundary_type='reflective')
    ymax = openmc.YPlane(+half, boundary_type='reflective')
    zmin = openmc.ZPlane(-half, boundary_type='reflective')
    zmax = openmc.ZPlane(+half, boundary_type='reflective')

    box_region = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax

    # Pebble cell: sphere filled with pebble universe
    pebble_cell = openmc.Cell(
        name='Pebble',
        fill=pebble_univ,
        region=-pebble_sphere
    )
    # Coolant cell: everything outside the pebble but inside the box
    # This represents the helium flowing between pebbles.
    coolant_cell = openmc.Cell(
        name='Helium Coolant',
        fill=mats['helium'],
        region=+pebble_sphere & box_region
    )

    root_univ = openmc.Universe(cells=[pebble_cell, coolant_cell])
    model.geometry = openmc.Geometry(root_univ)

    # -------------------------------------------------------------------------
    # Settings
    # -------------------------------------------------------------------------
    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    # Source at center of pebble
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0.0, 0.0, 0.0))
    )
    settings.temperature = {'method': 'interpolation'}
    model.settings = settings

    # -------------------------------------------------------------------------
    # Collect all materials
    # -------------------------------------------------------------------------
    model.materials = openmc.Materials(list(mats.values()))

    return model


# =============================================================================
# Core Model: Cylindrical Core with Homogenized Pebbles
# =============================================================================

def build_core_model(mats, particles, batches, inactive):
    """
    Build a simplified cylindrical HTR-10 core model with homogenized pebbles.

    In this model, the TRISO particles within each pebble are homogenized
    into a single effective fuel zone material (volume-weighted mixture).
    The pebbles still have two regions: the homogenized fuel zone (r=2.5 cm)
    and the graphite shell (2.5 < r < 3.0 cm). These pebbles are randomly
    packed into the cylindrical core cavity.

    The core is surrounded by a graphite reflector. The helium coolant fills
    the interstitial spaces between pebbles.

    This model captures the second level of heterogeneity (pebble packing)
    but not the first (TRISO self-shielding). For criticality calculations,
    the homogenization introduces a reactivity bias (typically a few hundred
    pcm) because the resonance self-shielding of discrete TRISO particles
    is lost.

    Parameters
    ----------
    mats : dict
        Material dictionary from create_materials().
    particles : int
        Neutron histories per batch.
    batches : int
        Total number of batches.
    inactive : int
        Number of inactive batches.

    Returns
    -------
    openmc.Model
        Complete OpenMC model.
    """
    model = openmc.Model()

    # -------------------------------------------------------------------------
    # Homogenized fuel zone material
    # -------------------------------------------------------------------------
    print("Computing homogenized fuel zone material...")
    homo_fuel = create_homogenized_fuel_zone(mats)

    # -------------------------------------------------------------------------
    # Pebble universe (homogenized fuel zone + graphite shell)
    # -------------------------------------------------------------------------
    fuel_zone_sphere = openmc.Sphere(r=FUEL_ZONE_RADIUS)
    pebble_sphere = openmc.Sphere(r=PEBBLE_OUTER_RADIUS)

    fuel_zone_cell = openmc.Cell(
        name='Homogenized Fuel Zone',
        fill=homo_fuel,
        region=-fuel_zone_sphere
    )
    shell_cell = openmc.Cell(
        name='Graphite Shell',
        fill=mats['shell'],
        region=+fuel_zone_sphere & -pebble_sphere
    )

    pebble_univ = openmc.Universe(cells=[fuel_zone_cell, shell_cell])

    # -------------------------------------------------------------------------
    # Core cavity: cylinder packed with pebbles
    # -------------------------------------------------------------------------
    # The core is a vertical cylinder of radius 90 cm and height ~126 cm.
    core_cyl = openmc.ZCylinder(r=CORE_RADIUS)
    core_bottom = openmc.ZPlane(z0=0.0)
    core_top = openmc.ZPlane(z0=CORE_HEIGHT)
    core_region = -core_cyl & +core_bottom & -core_top

    # Pack pebbles into the cylindrical core.
    # At packing fraction 0.61 with pebble radius 3.0 cm:
    #   N = PF * V_core / V_pebble
    #   V_core = pi * 90^2 * 126 = ~3,206,000 cm^3
    #   V_pebble = (4/3)*pi*3^3 = 113.1 cm^3
    #   N = 0.61 * 3,206,000 / 113.1 = ~17,290 pebbles
    #
    # Packing 17,000+ pebbles takes significant time and memory.
    # For tractability, we reduce the packing fraction slightly or use
    # fewer pebbles. The pack_spheres function handles this automatically.
    print("Packing pebbles into core cavity...")
    print("  (This may take a few minutes for close random packing...)")
    pebble_centers = openmc.model.pack_spheres(
        radius=PEBBLE_OUTER_RADIUS,
        region=core_region,
        pf=PEBBLE_PACKING_FRACTION,
        seed=123
    )
    n_pebbles = len(pebble_centers)
    print(f"  Packed {n_pebbles} pebbles into core")

    # -------------------------------------------------------------------------
    # Create pebble cells in the core
    # -------------------------------------------------------------------------
    # Each pebble is a cell filled with the pebble universe, translated to
    # its packed center position.
    #
    # For memory efficiency with many pebbles, we use a lattice-based approach.
    # We create a rectangular lattice overlying the core and assign each pebble
    # to the appropriate lattice cell(s). This is analogous to the TRISO lattice
    # approach but at the pebble scale.

    # Create TRISO-like objects for pebbles (TRISO class works for any sphere)
    pebbles = [openmc.model.TRISO(PEBBLE_OUTER_RADIUS, pebble_univ, c)
               for c in pebble_centers]

    # Create a lattice for efficient tracking of pebbles
    ll, ur = core_region.bounding_box
    # Use a reasonable lattice shape -- each cell should be a few pebble
    # diameters wide for good performance
    n_cells_xy = max(4, int((ur[0] - ll[0]) / (4 * PEBBLE_OUTER_RADIUS)))
    n_cells_z = max(4, int((ur[2] - ll[2]) / (4 * PEBBLE_OUTER_RADIUS)))
    shape = (n_cells_xy, n_cells_xy, n_cells_z)
    pitch = (ur - ll) / shape
    print(f"  Pebble lattice shape: {shape}, pitch: "
          f"({pitch[0]:.2f}, {pitch[1]:.2f}, {pitch[2]:.2f}) cm")

    pebble_lattice = openmc.model.create_triso_lattice(
        pebbles, ll, pitch, shape, mats['helium']
    )

    # Core cell filled with pebble lattice
    core_cell = openmc.Cell(
        name='Pebble Bed Core',
        fill=pebble_lattice,
        region=core_region
    )

    # -------------------------------------------------------------------------
    # Graphite reflector surrounding the core
    # -------------------------------------------------------------------------
    refl_cyl = openmc.ZCylinder(r=CORE_RADIUS + REFLECTOR_THICKNESS)
    refl_bottom = openmc.ZPlane(z0=-REFLECTOR_THICKNESS)
    refl_top = openmc.ZPlane(z0=CORE_HEIGHT + REFLECTOR_THICKNESS)

    # Vacuum boundary
    refl_cyl.boundary_type = 'vacuum'
    refl_bottom.boundary_type = 'vacuum'
    refl_top.boundary_type = 'vacuum'

    # Reflector region: everything inside the outer boundary but outside core
    reflector_cell = openmc.Cell(
        name='Graphite Reflector',
        fill=mats['reflector'],
        region=-refl_cyl & +refl_bottom & -refl_top & ~core_region
    )

    root_univ = openmc.Universe(cells=[core_cell, reflector_cell])
    model.geometry = openmc.Geometry(root_univ)

    # -------------------------------------------------------------------------
    # Settings
    # -------------------------------------------------------------------------
    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    # Initial source distributed throughout the core
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            (0.0, 0.0, 0.0),
            (CORE_RADIUS * 0.8, 0.0, CORE_HEIGHT),
        ),
        constraints={'fissionable': True}
    )
    settings.temperature = {'method': 'interpolation'}
    model.settings = settings

    # -------------------------------------------------------------------------
    # Materials
    # -------------------------------------------------------------------------
    all_mats = list(mats.values()) + [homo_fuel]
    model.materials = openmc.Materials(all_mats)

    return model


# =============================================================================
# Main
# =============================================================================

def build_model(model_type='pin', particles=10000, batches=120, inactive=20):
    """
    Top-level function to build the requested HTR-10 model variant.

    Parameters
    ----------
    model_type : {'pin', 'core'}
        Which model to build (see module docstring).
    particles : int
        Neutron histories per batch.
    batches : int
        Total number of batches.
    inactive : int
        Number of inactive batches for source convergence.

    Returns
    -------
    openmc.Model
        Complete OpenMC model ready for export or execution.
    """
    print(f"Building HTR-10 model (mode: {model_type})")
    print("=" * 60)

    mats = create_materials()

    if model_type == 'pin':
        model = build_pin_model(mats, particles, batches, inactive)
    elif model_type == 'core':
        model = build_core_model(mats, particles, batches, inactive)
    else:
        raise ValueError(f"Unknown model type: {model_type!r}. Use 'pin' or 'core'.")

    print("=" * 60)
    print("Model construction complete.")
    return model


def main():
    parser = argparse.ArgumentParser(
        description='HTR-10 Pebble Bed Reactor - OpenMC Model'
    )
    parser.add_argument('--particles', type=int, default=10000,
                        help='Neutron histories per batch (default: 10000)')
    parser.add_argument('--batches', type=int, default=120,
                        help='Total number of batches (default: 120)')
    parser.add_argument('--inactive', type=int, default=20,
                        help='Number of inactive batches (default: 20)')
    parser.add_argument('--model', type=str, default='pin',
                        choices=['pin', 'core'],
                        help='Model type: pin (single pebble) or core '
                             '(simplified cylindrical core)')
    parser.add_argument('--run', action='store_true',
                        help='Run OpenMC after building the model')
    args = parser.parse_args()

    model = build_model(
        model_type=args.model,
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive
    )

    # Export to XML files
    model.export_to_model_xml()
    print("Exported model.xml")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
