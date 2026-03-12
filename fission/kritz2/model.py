#!/usr/bin/env python3
"""
KRITZ-2 Critical Experiments -- OpenMC Pin Cell Model
=====================================================

This script builds an infinite-lattice (pin cell) model for the KRITZ-2
critical experiments performed at Studsvik, Sweden in the early 1970s.

Background
----------
The KRITZ reactor was a zero-power critical facility designed to measure
lattice physics parameters for light-water reactor fuel.  Three KRITZ-2
configurations used regular square lattices of identical fuel rods in
light water, with criticality achieved by adjusting both the boron
concentration and the water level:

    KRITZ-2:1   UO2 fuel, 1.86 wt% U-235, 14.85 mm pitch, 44x44 rods
    KRITZ-2:13  UO2 fuel, 1.86 wt% U-235, 16.35 mm pitch, 40x40 rods
    KRITZ-2:19  MOX fuel (weapons-grade Pu), 18.00 mm pitch, 25x24 rods

Each configuration was measured at room temperature (~20 C) and at an
elevated temperature (~245 C).  The temperature reactivity effect is the
principal quantity of interest -- it tests the accuracy of Doppler
broadening, thermal scattering S(alpha,beta) data, and water density
changes.

Physics notes
-------------
* At room temperature the moderator is liquid water near 1.0 g/cc.  At
  the hot condition (~245 C, ~3.5 MPa) the water density drops to about
  0.8 g/cc, hardening the neutron spectrum and reducing moderation.
  Criticality at the hot condition therefore requires a much taller water
  column and lower boron concentration.

* The pin cell calculation with reflective boundary conditions represents
  an infinite lattice.  This eliminates leakage effects and isolates the
  lattice-level neutron balance, making it the standard test for cross-
  section validation.

* Boron in the moderator acts as a thermal neutron absorber (B-10 has a
  very large 1/v absorption cross section, ~3840 barns at 0.025 eV).
  The boron concentration varies between configurations and temperature
  conditions to tune criticality.

* Zircaloy-2 cladding is nearly transparent to thermal neutrons
  (Zr has a very small absorption cross section), which is why it is
  the standard cladding material in BWRs and CANDU reactors.

Reference
---------
OECD/NEA, "Benchmark on the KRITZ-2 LEU and MOX Critical Experiments,"
NEA/NSC/DOC(2005)24, 2006.

Reference k-eff values (average of 13 participant solutions):
    KRITZ-2:1  cold  0.99463 +/- 0.00330
    KRITZ-2:1  hot   0.99290 +/- 0.00346
    KRITZ-2:13 cold  0.99736 +/- 0.00204
    KRITZ-2:13 hot   0.99605 +/- 0.00280
    KRITZ-2:19 cold  0.99908
    KRITZ-2:19 hot   0.99822
"""

import argparse
import numpy as np
import openmc


# ---------------------------------------------------------------------------
# Configuration data from NEA/NSC/DOC(2005)24, Appendix A, Tables 11-12
# ---------------------------------------------------------------------------
# All dimensions in cm (converted from mm in the specification).
# Atom densities in units of 1e24 atoms/cm^3 (i.e. atoms/barn-cm).
# Weight fractions are dimensionless.

CONFIGS = {
    # -----------------------------------------------------------------------
    # KRITZ-2:1  --  UO2, tight pitch
    # -----------------------------------------------------------------------
    '2:1-cold': {
        'description': 'KRITZ-2:1 at 19.7 C (room temperature)',
        'temperature_C': 19.7,
        'temperature_K': 293.0,      # ~19.7 C, use 293 K for cross sections
        'fuel_type': 'UO2',
        # Geometry (Table 11 -- at 19.7 C)
        'fuel_radius': 0.529,        # DF/2 = 10.58/2 mm = 0.529 cm
        'clad_inner_radius': 0.529,  # No gap specified; clad sits on fuel
        'clad_outer_radius': 0.6125, # DR/2 = 12.25/2 mm = 0.6125 cm
        'pitch': 1.485,              # 14.85 mm = 1.485 cm
        'n_rods': '44 x 44',
        # UO2 fuel (Table 12, Appendix A.1)
        # Density: 10.145 g/cc (homogenised pellet + dish volume)
        # Atom densities (1e24 atoms/cm^3):
        #   U-235: 4.22966e-4   U-234: 3.40280e-6   U-236: 0.0
        #   U-238: 2.20319e-2   O:     4.49165e-2
        'fuel_density': 10.145,
        'fuel_atom_densities': {
            'U235': 4.22966e-4,
            'U234': 3.40280e-6,
            'U238': 2.20319e-2,
            'O16':  4.49165e-2,
        },
        # Zircaloy-2 cladding (Table 12, generic composition)
        # Density: 6.506 g/cc
        # Weight fractions: Zr 98.225%, Sn 1.5%, Fe 0.125%, Cr 0.1%, N 0.05%
        'clad_density': 6.506,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator: light water + boron at 19.7 C
        # Boron: 217.9 ppm;  Water density: 0.99822 g/cc
        # Atom densities (1e24 atoms/cm^3):
        #   H:    6.67223e-2   O:    3.33612e-2
        #   B-10: 2.41113e-6   B-11: 9.70509e-6
        'mod_density': 0.99822,
        'mod_boron_ppm': 217.9,
        'mod_atom_densities': {
            'H1':  6.67223e-2,
            'O16': 3.33612e-2,
            'B10': 2.41113e-6,
            'B11': 9.70509e-6,
        },
        # Reference k-eff (average of 13 solutions)
        'ref_keff': 0.99463,
        'ref_keff_unc': 0.00330,
    },

    '2:1-hot': {
        'description': 'KRITZ-2:1 at 248.5 C (elevated temperature)',
        'temperature_C': 248.5,
        'temperature_K': 521.65,     # 248.5 C
        'fuel_type': 'UO2',
        # Geometry at 248.5 C (thermally expanded, Table 11)
        'fuel_radius': 0.53033,      # 10.6066/2 mm
        'clad_inner_radius': 0.53033,
        'clad_outer_radius': 0.61348,# 12.2696/2 mm
        'pitch': 1.49112,            # 14.9112 mm
        'n_rods': '44 x 44',
        # Fuel atom densities at 248.5 C (Table 12, Appendix A.1)
        'fuel_density': 10.06879,
        'fuel_atom_densities': {
            'U235': 4.26167e-4,
            'U234': 3.42855e-6,
            'U238': 2.21987e-2,
            'O16':  4.52565e-2,
        },
        # Zircaloy-2 at elevated temperature
        'clad_density': 6.47484,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator at 248.5 C, 26.2 ppm boron
        # Density: 0.80111 g/cc (hot pressurised water -- much lower than cold)
        'mod_density': 0.80111,
        'mod_boron_ppm': 26.2,
        'mod_atom_densities': {
            'H1':  5.35575e-2,
            'O16': 2.67788e-2,
            'B10': 2.32665e-7,
            'B11': 9.36504e-7,
        },
        'ref_keff': 0.99290,
        'ref_keff_unc': 0.00346,
    },

    # -----------------------------------------------------------------------
    # KRITZ-2:13  --  UO2, wider pitch (same fuel rods as 2:1)
    # -----------------------------------------------------------------------
    '2:13-cold': {
        'description': 'KRITZ-2:13 at 22.1 C (room temperature)',
        'temperature_C': 22.1,
        'temperature_K': 295.25,
        'fuel_type': 'UO2',
        # Geometry at 22.1 C (Table 11, Appendix A.1 continued)
        'fuel_radius': 0.529,        # Same fuel rod as 2:1
        'clad_inner_radius': 0.529,
        'clad_outer_radius': 0.6125,
        'pitch': 1.635,              # 16.35 mm -- wider lattice
        'n_rods': '40 x 40',
        # Fuel -- same UO2 as KRITZ-2:1
        'fuel_density': 10.145,
        'fuel_atom_densities': {
            'U235': 4.23076e-4,
            'U234': 3.40368e-6,
            'U238': 2.20376e-2,
            'O16':  4.49282e-2,
        },
        'clad_density': 6.506,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator at 22.1 C, 451.9 ppm boron (higher boron because
        # the wider pitch gives more moderation and higher reactivity)
        'mod_density': 0.99771,
        'mod_boron_ppm': 451.9,
        'mod_atom_densities': {
            'H1':  6.66726e-2,
            'O16': 3.33363e-2,
            'B10': 4.99785e-6,
            'B11': 2.01170e-5,
        },
        'ref_keff': 0.99736,
        'ref_keff_unc': 0.00204,
    },

    '2:13-hot': {
        'description': 'KRITZ-2:13 at 243.0 C (elevated temperature)',
        'temperature_C': 243.0,
        'temperature_K': 516.15,
        'fuel_type': 'UO2',
        # Geometry at 243.0 C (Table 11 continued)
        'fuel_radius': 0.530285,     # 10.6057/2 mm
        'clad_inner_radius': 0.530285,
        'clad_outer_radius': 0.61345,# 12.2689/2 mm
        'pitch': 1.6415,             # 16.415 mm (thermally expanded)
        'n_rods': '40 x 40',
        'fuel_density': 10.06879,
        'fuel_atom_densities': {
            'U235': 4.26167e-4,
            'U234': 3.42855e-6,
            'U238': 2.21987e-2,
            'O16':  4.52565e-2,
        },
        'clad_density': 6.47484,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator at 243.0 C, 280.1 ppm boron
        'mod_density': 0.80910,
        'mod_boron_ppm': 280.1,
        'mod_atom_densities': {
            'H1':  5.39000e-2,
            'O16': 2.69500e-2,
            'B10': 2.55000e-6,
            'B11': 1.02600e-5,
        },
        'ref_keff': 0.99605,
        'ref_keff_unc': 0.00280,
    },

    # -----------------------------------------------------------------------
    # KRITZ-2:19  --  MOX fuel, weapons-grade Pu
    # -----------------------------------------------------------------------
    '2:19-cold': {
        'description': 'KRITZ-2:19 at 21.1 C (room temperature, MOX)',
        'temperature_C': 21.1,
        'temperature_K': 294.25,
        'fuel_type': 'MOX',
        # Geometry at 21.1 C (Table 11, Appendix A.2)
        # MOX rods are smaller than the UO2 rods
        'fuel_radius': 0.4725,       # 9.45/2 mm = 0.4725 cm
        'clad_inner_radius': 0.4725, # No gap; fuel fills cladding
        'clad_outer_radius': 0.5395, # 10.79/2 mm = 0.5395 cm
        'pitch': 1.800,              # 18.00 mm = 1.800 cm
        'n_rods': '25 x 24',
        # MOX fuel: vibrocompacted PuO2-UO2
        # Density: 9.58 g/cc
        # Pu composition (decayed to measurement date):
        #   Pu-239: 91.41 at%, Pu-240: 7.83 at%, Pu-241: 0.73 at%, Pu-242: 0.03 at%
        # (Pu-241 decayed by factor 0.6070 from 1962 to 1972; Am-241 produced)
        # Atom densities from Table 11, Appendix A.2:
        'fuel_density': 9.58,
        'fuel_atom_densities': {
            'O16':   4.27252e-2,
            'U235':  3.40995e-5,
            'U238':  2.10093e-2,
            'Pu239': 2.91742e-4,
            'Pu240': 2.49901e-5,
            'Pu241': 1.41422e-6,
            'Pu242': 9.57474e-8,
            'Am241': 9.15633e-7,
        },
        'clad_density': 6.506,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator at 21.1 C, 4.8 ppm boron (very low -- MOX has lower
        # reactivity than enriched UO2 per unit fissile content because
        # of the Pu-240 resonance absorption and the harder fission spectrum)
        'mod_density': 0.99793,
        'mod_boron_ppm': 4.8,
        'mod_atom_densities': {
            'H1':  6.67172e-2,
            'O16': 3.33586e-2,
            'B10': 5.30980e-8,
            'B11': 2.13726e-7,
        },
        'ref_keff': 0.99908,
        'ref_keff_unc': None,  # Not reported in the benchmark document
    },

    '2:19-hot': {
        'description': 'KRITZ-2:19 at 235.9 C (elevated temperature, MOX)',
        'temperature_C': 235.9,
        'temperature_K': 509.05,
        'fuel_type': 'MOX',
        # Geometry at 235.9 C (Table 11, Appendix A.2)
        'fuel_radius': 0.47385,      # 9.477/2 mm (thermally expanded)
        'clad_inner_radius': 0.47385,
        'clad_outer_radius': 0.54075,# 10.815/2 mm
        'pitch': 1.80696,            # 18.0696 mm
        'n_rods': '25 x 24',
        'fuel_density': 9.55126,
        'fuel_atom_densities': {
            'O16':   4.25970e-2,
            'U235':  3.39972e-5,
            'U238':  2.09463e-2,
            'Pu239': 2.90867e-4,
            'Pu240': 2.49151e-5,
            'Pu241': 1.40998e-6,
            'Pu242': 9.54602e-8,
            'Am241': 9.12885e-7,
        },
        'clad_density': 6.47484,
        'clad_weight_fractions': {
            'Zr': 0.982250,
            'Sn': 0.015000,
            'Fe': 0.001250,
            'Cr': 0.001000,
            'N14': 0.000500,
        },
        # Moderator at 235.9 C, 5.2 ppm boron
        'mod_density': 0.82000,
        'mod_boron_ppm': 5.2,
        'mod_atom_densities': {
            'H1':  5.46000e-2,
            'O16': 2.73000e-2,
            'B10': 4.60000e-8,
            'B11': 1.85000e-7,
        },
        'ref_keff': 0.99822,
        'ref_keff_unc': None,
    },
}


def build_model(config_name, particles, batches, inactive):
    """
    Build an OpenMC pin cell model for the specified KRITZ-2 configuration.

    The model consists of three concentric regions:
      1. Fuel pellet (UO2 or MOX)
      2. Zircaloy-2 cladding
      3. Moderator (borated light water)

    Reflective boundary conditions on all six faces make this an infinite
    lattice calculation, which is the standard approach for benchmarking
    lattice-level nuclear data and methods.

    Parameters
    ----------
    config_name : str
        Configuration key, e.g. '2:1-cold'
    particles : int
        Number of neutron histories per batch
    batches : int
        Total number of batches
    inactive : int
        Number of inactive (source convergence) batches

    Returns
    -------
    openmc.Model
        Fully configured OpenMC model
    """
    cfg = CONFIGS[config_name]
    print(f"Building model: {cfg['description']}")
    print(f"  Fuel type:     {cfg['fuel_type']}")
    print(f"  Temperature:   {cfg['temperature_C']} C ({cfg['temperature_K']} K)")
    print(f"  Lattice pitch: {cfg['pitch']} cm")
    print(f"  Boron:         {cfg['mod_boron_ppm']} ppm")
    print(f"  Reference keff: {cfg['ref_keff']}")

    # ===================================================================
    # Materials
    # ===================================================================

    # --- Fuel ---
    # For UO2: enriched uranium dioxide at ~10.15 g/cc
    # For MOX: vibrocompacted PuO2-UO2 at ~9.58 g/cc
    # We use explicit atom densities from the benchmark specification
    # to ensure exact reproduction of the reference problem.
    fuel = openmc.Material(name=f'Fuel ({cfg["fuel_type"]})')
    fuel.temperature = cfg['temperature_K']
    # Set atom densities directly from the benchmark tables.
    # These were computed by the benchmark organisers accounting for
    # the exact isotopic composition and pellet/dish homogenisation.
    for nuclide, density in cfg['fuel_atom_densities'].items():
        fuel.add_nuclide(nuclide, density)
    fuel.set_density('atom/b-cm', sum(cfg['fuel_atom_densities'].values()))

    # --- Cladding (Zircaloy-2) ---
    # Zircaloy-2 is the standard BWR cladding material:
    #   ~98.2% Zr, ~1.5% Sn, ~0.125% Fe, ~0.1% Cr, ~0.05% N
    # Zirconium has a very low thermal neutron absorption cross section
    # (~0.18 barns), making it nearly transparent to thermal neutrons.
    # The small amounts of Sn, Fe, Cr provide corrosion resistance.
    clad = openmc.Material(name='Zircaloy-2 Cladding')
    clad.temperature = cfg['temperature_K']
    for element, wt_frac in cfg['clad_weight_fractions'].items():
        if element == 'N14':
            clad.add_nuclide('N14', wt_frac, 'wo')
        else:
            clad.add_element(element, wt_frac, 'wo')
    clad.set_density('g/cm3', cfg['clad_density'])

    # --- Moderator (borated light water) ---
    # Light water is the primary moderator and coolant.  Hydrogen's large
    # scattering cross section (~20 barns) and low mass number make it
    # the most effective moderator per atom for thermalising neutrons.
    # Boron (specifically B-10 with its 3840 barn thermal absorption
    # cross section) is dissolved in the water as boric acid to control
    # reactivity.  The B-10/B-11 ratio follows natural abundance (19.9/80.1).
    moderator = openmc.Material(name=f'Borated Water ({cfg["mod_boron_ppm"]} ppm B)')
    moderator.temperature = cfg['temperature_K']
    for nuclide, density in cfg['mod_atom_densities'].items():
        moderator.add_nuclide(nuclide, density)
    moderator.set_density('atom/b-cm', sum(cfg['mod_atom_densities'].values()))
    # Add S(alpha,beta) thermal scattering data for hydrogen bound in
    # water.  This is critical for accurate thermalisation -- free-gas
    # treatment of hydrogen would be grossly incorrect below ~4 eV
    # because it ignores the molecular binding and translational modes
    # of the water molecule.
    moderator.add_s_alpha_beta('c_H_in_H2O')

    materials = openmc.Materials([fuel, clad, moderator])

    # ===================================================================
    # Geometry
    # ===================================================================
    # The pin cell is the fundamental repeating unit of the square lattice.
    # It consists of concentric cylinders (fuel, clad) inside a square
    # region filled with moderator.
    #
    #    +---------------------------+
    #    |          moderator        |
    #    |     +--clad outer--+      |
    #    |     | +--fuel----+ |      |
    #    |     | |          | |      |
    #    |     | +----------+ |      |
    #    |     +--------------+      |
    #    +---------------------------+
    #           <--- pitch --->

    # Cylindrical surface for the fuel pellet outer radius.
    # In this benchmark there is no explicit gap between fuel and
    # cladding -- the specification states fuel diameter equals
    # cladding inner diameter (no gap gas is modelled).
    fuel_or = openmc.ZCylinder(
        r=cfg['fuel_radius'],
        name='Fuel pellet outer radius'
    )

    # Cylindrical surface for the cladding outer radius.
    # Cladding wall thickness = clad_outer_radius - fuel_radius
    # For KRITZ-2:1: 0.6125 - 0.529 = 0.0835 cm = 0.835 mm (= 0.74 mm
    # from the spec -- the slight difference is because fuel_radius is
    # given as pellet radius, not clad inner radius; they are equal here).
    clad_or = openmc.ZCylinder(
        r=cfg['clad_outer_radius'],
        name='Cladding outer radius'
    )

    # Square boundary for the pin cell -- reflective BCs to simulate
    # an infinite lattice.  The half-width is pitch/2.
    half_pitch = cfg['pitch'] / 2.0
    # The pin cell boundary is a rectangular prism (box).
    # Reflective boundary conditions on all surfaces mean that any
    # neutron reaching the boundary is specularly reflected back,
    # which is equivalent to an infinite repeating lattice.
    left = openmc.XPlane(-half_pitch, boundary_type='reflective',
                         name='Pin cell left boundary')
    right = openmc.XPlane(half_pitch, boundary_type='reflective',
                          name='Pin cell right boundary')
    bottom = openmc.YPlane(-half_pitch, boundary_type='reflective',
                           name='Pin cell bottom boundary')
    top = openmc.YPlane(half_pitch, boundary_type='reflective',
                        name='Pin cell top boundary')

    # --- Cells ---

    # Cell 1: Fuel region (inside the fuel pellet)
    # This is where fission occurs.  For UO2 fuel, U-235 undergoes
    # thermal fission (sigma_f ~ 585 barns at 0.025 eV) and U-238
    # undergoes fast fission (threshold ~ 1 MeV) and resonance capture.
    # For MOX fuel, Pu-239 is the primary fissile isotope (sigma_f ~
    # 748 barns at 0.025 eV), with contributions from Pu-241.
    fuel_cell = openmc.Cell(name='Fuel pellet')
    fuel_cell.fill = fuel
    fuel_cell.region = -fuel_or  # Inside the fuel cylinder

    # Cell 2: Cladding region (annular ring between fuel and clad OD)
    # The Zircaloy-2 is thin (~0.74 mm) and nearly transparent to
    # neutrons, but it does contribute a small amount of parasitic
    # absorption, particularly from the Hf impurity (not modelled
    # in this generic Zircaloy-2 composition).
    clad_cell = openmc.Cell(name='Cladding')
    clad_cell.fill = clad
    clad_cell.region = +fuel_or & -clad_or  # Between fuel OR and clad OR

    # Cell 3: Moderator region (everything outside the clad, inside the box)
    # The water fills the space between fuel rods.  The moderator-to-fuel
    # ratio (proportional to pitch^2 - pi*R_clad^2) determines the
    # neutron spectrum -- a higher ratio means more thermalisation and
    # typically higher k-inf (up to a point, then over-moderation occurs).
    mod_cell = openmc.Cell(name='Moderator')
    mod_cell.fill = moderator
    mod_cell.region = +clad_or & +left & -right & +bottom & -top

    # Create the root universe containing all three cells
    root_universe = openmc.Universe(cells=[fuel_cell, clad_cell, mod_cell])
    geometry = openmc.Geometry(root_universe)

    # ===================================================================
    # Settings
    # ===================================================================
    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles

    # Eigenvalue (criticality) mode -- this is the default for OpenMC
    # when no fixed source is specified.  The Monte Carlo process:
    # 1. Start with an initial source distribution
    # 2. Track neutrons through the geometry, scoring fission sites
    # 3. Use fission sites as the source for the next generation
    # 4. After inactive batches (source convergence), begin tallying k-eff
    settings.run_mode = 'eigenvalue'

    # Initial source: uniform distribution within the fuel pellet.
    # For a simple pin cell, a point source at the center also works,
    # but a spatial distribution converges the fission source faster.
    spatial_dist = openmc.stats.Box(
        lower_left=(-cfg['fuel_radius'], -cfg['fuel_radius'], -1.0),
        upper_right=(cfg['fuel_radius'], cfg['fuel_radius'], 1.0),
    )
    settings.source = openmc.IndependentSource(
        space=spatial_dist,
        constraints={'fissionable': True}
    )

    # Temperature treatment: when the temperature is specified on
    # materials, OpenMC will look for cross sections at that temperature
    # or interpolate if multipole data / temperature interpolation is
    # enabled.  For room temperature cases (~293 K), most libraries
    # have data at 293.6 K which is an excellent match.
    settings.temperature = {'method': 'interpolation'}

    # ===================================================================
    # Assemble the model
    # ===================================================================
    model = openmc.Model(geometry, materials, settings)
    return model


def main():
    parser = argparse.ArgumentParser(
        description='KRITZ-2 Critical Experiment Pin Cell Model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Configurations:
  2:1-cold   UO2, 14.85 mm pitch, 19.7 C,  217.9 ppm B  (default)
  2:1-hot    UO2, 14.85 mm pitch, 248.5 C,  26.2 ppm B
  2:13-cold  UO2, 16.35 mm pitch, 22.1 C,  451.9 ppm B
  2:13-hot   UO2, 16.35 mm pitch, 243.0 C,  280.1 ppm B
  2:19-cold  MOX, 18.00 mm pitch, 21.1 C,    4.8 ppm B
  2:19-hot   MOX, 18.00 mm pitch, 235.9 C,   5.2 ppm B

Reference: NEA/NSC/DOC(2005)24
        """
    )
    parser.add_argument(
        '--config', type=str, default='2:1-cold',
        choices=['2:1-cold', '2:1-hot', '2:13-cold', '2:13-hot',
                 '2:19-cold', '2:19-hot'],
        help='KRITZ-2 configuration to model (default: 2:1-cold)'
    )
    parser.add_argument(
        '--particles', type=int, default=10000,
        help='Number of neutron histories per batch (default: 10000)'
    )
    parser.add_argument(
        '--batches', type=int, default=150,
        help='Total number of batches (default: 150)'
    )
    parser.add_argument(
        '--inactive', type=int, default=50,
        help='Number of inactive batches for source convergence (default: 50)'
    )
    args = parser.parse_args()

    model = build_model(args.config, args.particles, args.batches,
                        args.inactive)

    # Export the model to XML files (geometry.xml, materials.xml, settings.xml)
    model.export_to_model_xml()
    print(f"\nModel exported to model.xml")
    print(f"Run with: openmc")


if __name__ == '__main__':
    main()
