#!/usr/bin/env python3
"""
===============================================================================
MSRE - Molten Salt Reactor Experiment
OpenMC Benchmark Model
===============================================================================

Historical Background
---------------------
The Molten Salt Reactor Experiment (MSRE) was designed, built, and operated at
Oak Ridge National Laboratory (ORNL) between 1964 and 1969. It was the world's
first molten-salt-fueled reactor to achieve sustained criticality and power
operation. The MSRE demonstrated the technical feasibility of the molten salt
reactor concept, a design in which the nuclear fuel is dissolved directly in a
circulating fluoride salt rather than being fabricated into solid fuel elements.

The MSRE was a thermal-spectrum reactor with a design power of 7.4 MWth. The
fuel salt -- a mixture of lithium fluoride, beryllium fluoride, zirconium
fluoride, and uranium fluoride (LiF-BeF2-ZrF4-UF4) -- flowed upward through
channels formed between vertical graphite moderator stringers. The graphite
served as the neutron moderator, thermalizing fast fission neutrons so they
could sustain the chain reaction in the dissolved uranium. After passing through
the core, the hot fuel salt exited the reactor vessel, flowed through a
primary heat exchanger where it transferred heat to a secondary coolant salt
(LiF-BeF2), and was then pumped back into the bottom of the reactor vessel.

Molten Salt Reactor Physics
---------------------------
Molten salt reactors (MSRs) present unique neutronics challenges compared to
solid-fueled reactors:

  1. LIQUID FUEL: The fissile material (here UF4 dissolved in the carrier salt)
     is a flowing liquid. This means delayed neutron precursors are carried out
     of the core by the circulating salt, reducing the effective delayed neutron
     fraction and making the reactor respond somewhat faster to reactivity
     changes. In the MSRE, roughly half the delayed neutron precursors decayed
     outside the core.

  2. LITHIUM-7 ENRICHMENT: Natural lithium contains ~7.5% Li-6, which has a
     very large thermal neutron absorption cross section (940 barns). Using
     natural lithium in the fuel salt would result in an unacceptably large
     parasitic neutron loss. The MSRE therefore used lithium enriched to
     99.99% Li-7. Even small amounts of Li-6 contamination significantly
     affect the neutron economy.

  3. FLUORINE MODERATION: The fluorine in the salt contributes some neutron
     moderation (fluorine has mass 19, intermediate between hydrogen and
     carbon). The F-19(n,alpha) reaction is also a neutron loss mechanism at
     higher energies.

  4. GRAPHITE MODERATOR: The CGB-grade graphite stringers provided the primary
     neutron moderation. The graphite was unclad and in direct contact with
     the fuel salt. Over time, fission products (especially xenon and krypton)
     could diffuse into the graphite pores, affecting reactivity.

  5. NEGATIVE TEMPERATURE COEFFICIENTS: The fuel salt has a negative temperature
     coefficient of reactivity due to thermal expansion -- as the salt heats
     up it becomes less dense, reducing the fuel density in the core and
     decreasing reactivity. This provides inherent passive safety.

Benchmark Reference
-------------------
This model is based on the IRPhEP benchmark MSRE-MSR-EXP-001 and the
comparative study published in Frontiers in Nuclear Engineering (2024):
  https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2024.1385478/full

Two model variants are provided:
  - "homogeneous": A fully homogenized cylindrical core where the graphite
    moderator and fuel salt are mixed by volume fraction into a single
    effective material. This is the simplest representation and is useful
    for scoping calculations.
  - "channel": A single fuel channel unit cell with reflective boundary
    conditions. A square graphite block surrounds a central rectangular
    fuel salt channel. This captures the heterogeneous self-shielding
    effects that the homogeneous model misses.

Reference k-eff Values:
  Experimental:            0.99978 +/- 0.00420
  OpenMC CSG (detailed):   1.00878 +/- 0.00032
  Serpent (ENDF/B-VII.1):  1.02132 +/- 0.00003

===============================================================================
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# Command-line argument parsing
# =============================================================================

def parse_args():
    """Parse command-line arguments for the MSRE benchmark model.

    The user can select the model type (homogeneous or channel unit cell)
    and adjust Monte Carlo simulation parameters.
    """
    parser = argparse.ArgumentParser(
        description="MSRE Molten Salt Reactor Experiment - OpenMC Benchmark Model"
    )
    parser.add_argument(
        "--model",
        choices=["homogeneous", "channel"],
        default="homogeneous",
        help=(
            "Model type: 'homogeneous' creates a cylindrical core with "
            "homogenized graphite+salt material; 'channel' creates a single "
            "fuel channel unit cell with reflective BCs. Default: homogeneous"
        ),
    )
    parser.add_argument(
        "--particles",
        type=int,
        default=10000,
        help="Number of neutron histories per batch (default: 10000)",
    )
    parser.add_argument(
        "--batches",
        type=int,
        default=150,
        help="Total number of batches (default: 150)",
    )
    parser.add_argument(
        "--inactive",
        type=int,
        default=50,
        help="Number of inactive (discarded) batches for source convergence (default: 50)",
    )
    return parser.parse_args()


# =============================================================================
# Material definitions
# =============================================================================

def create_fuel_salt():
    """Create the MSRE fuel salt material.

    The MSRE fuel salt composition is LiF-BeF2-ZrF4-UF4 with mole fractions
    of 65.0 - 29.1 - 5.0 - 0.9 mol%, respectively.

    The lithium is enriched to 99.99 at% Li-7 to minimize parasitic neutron
    absorption by Li-6 (thermal absorption cross section ~940 barns). Even
    a small fraction of Li-6 would significantly degrade the neutron economy
    and make the reactor subcritical.

    The uranium is enriched to approximately 33 at% U-235 in the initial
    U-235 fuel loading. The remainder is U-238, which provides some fertile
    conversion (breeding U-235 replacement via Pu-239, though at a very low
    conversion ratio in this thermal spectrum).

    The fuel salt density at the mean operating temperature (~922 K) is
    approximately 2.32 g/cm3.

    Composition breakdown by mole fraction of compound:
      LiF:   0.650  ->  contributes Li and F
      BeF2:  0.291  ->  contributes Be and 2*F per molecule
      ZrF4:  0.050  ->  contributes Zr and 4*F per molecule
      UF4:   0.009  ->  contributes U  and 4*F per molecule

    Total fluorine moles per mole of salt mixture:
      F = 0.650*1 + 0.291*2 + 0.050*4 + 0.009*4
      F = 0.650 + 0.582 + 0.200 + 0.036 = 1.468

    Total atom fractions (unnormalized):
      Li:  0.650
      Be:  0.291
      Zr:  0.050
      U:   0.009
      F:   1.468
    Total = 2.468

    Normalized atom fractions:
      Li:  0.650 / 2.468 = 0.26337
      Be:  0.291 / 2.468 = 0.11791
      Zr:  0.050 / 2.468 = 0.02026
      U:   0.009 / 2.468 = 0.00365
      F:   1.468 / 2.468 = 0.59481
    """
    fuel_salt = openmc.Material(name="MSRE Fuel Salt (LiF-BeF2-ZrF4-UF4)")

    # --- Lithium: 99.99% Li-7, 0.01% Li-6 ---
    # Li-6 has a thermal neutron absorption cross section of ~940 barns.
    # The extreme enrichment to Li-7 is absolutely essential for achieving
    # criticality. With natural lithium (~7.5% Li-6), the parasitic
    # absorption would be far too large.
    li_atom_frac = 0.26337
    li7_fraction = 0.9999  # 99.99 at% Li-7
    li6_fraction = 0.0001  # 0.01 at% Li-6
    fuel_salt.add_nuclide("Li7", li_atom_frac * li7_fraction)
    fuel_salt.add_nuclide("Li6", li_atom_frac * li6_fraction)

    # --- Beryllium-9 ---
    # Beryllium serves as an additional light moderator. Be-9 also
    # participates in (n,2n) reactions at higher energies, which can
    # be a small neutron source. BeF2 is the main component that gives
    # the salt mixture its favorable physical properties (low melting
    # point, good heat transfer, chemical stability).
    fuel_salt.add_nuclide("Be9", 0.11791)

    # --- Zirconium ---
    # ZrF4 is added to the salt to stabilize the UF4 in solution and
    # to prevent uranium oxide precipitation. Natural zirconium isotopic
    # composition is used.
    # Zr isotopic abundances (natural):
    #   Zr-90: 51.45%, Zr-91: 11.22%, Zr-92: 17.15%,
    #   Zr-94: 17.38%, Zr-96: 2.80%
    zr_atom_frac = 0.02026
    fuel_salt.add_nuclide("Zr90", zr_atom_frac * 0.5145)
    fuel_salt.add_nuclide("Zr91", zr_atom_frac * 0.1122)
    fuel_salt.add_nuclide("Zr92", zr_atom_frac * 0.1715)
    fuel_salt.add_nuclide("Zr94", zr_atom_frac * 0.1738)
    fuel_salt.add_nuclide("Zr96", zr_atom_frac * 0.0280)

    # --- Uranium: ~33 at% U-235, ~67 at% U-238 ---
    # The initial MSRE fuel loading used uranium enriched to approximately
    # 33% U-235. This relatively high enrichment was needed because the
    # uranium concentration in the salt is low (only 0.9 mol% UF4) and
    # the core has significant neutron leakage and parasitic absorption.
    # Later in the MSRE program, U-233 was substituted as fuel to
    # demonstrate the thorium fuel cycle, but this model represents the
    # initial U-235 loading.
    u_atom_frac = 0.00365
    u235_enrichment = 0.33  # 33 at% U-235
    u238_fraction = 1.0 - u235_enrichment
    fuel_salt.add_nuclide("U235", u_atom_frac * u235_enrichment)
    fuel_salt.add_nuclide("U238", u_atom_frac * u238_fraction)

    # --- Fluorine-19 ---
    # Fluorine is the most abundant element in the salt by atom fraction.
    # F-19 is 100% of natural fluorine (monoisotopic). Fluorine has a
    # relatively low thermal absorption cross section (~0.01 barns) so
    # it does not significantly impact the neutron economy. However, the
    # F-19(n,alpha)N-16 reaction at higher energies produces N-16 which
    # is a strong gamma emitter (6-7 MeV), contributing to the radiation
    # environment around the primary loop.
    fuel_salt.add_nuclide("F19", 0.59481)

    # --- Density ---
    # The fuel salt density at the mean operating temperature of ~922 K
    # (649 deg C) is approximately 2.32 g/cm3. The density decreases
    # with temperature (thermal expansion), which is the physical basis
    # for the negative fuel temperature coefficient of reactivity.
    fuel_salt.set_density("g/cm3", 2.32)

    # Set the temperature to the mean operating temperature.
    # The MSRE operated with an inlet temperature of ~905 K (632 deg C)
    # and an outlet temperature of ~940 K (667 deg C). We use the
    # average of ~922 K.
    fuel_salt.temperature = 922.0

    return fuel_salt


def create_graphite():
    """Create the MSRE graphite moderator material.

    The MSRE used CGB-grade graphite manufactured by Carbon Products Division
    of Union Carbide. This was a nuclear-grade graphite with:
      - High density (~1.86 g/cm3)
      - Low impurity content (especially low boron)
      - Fine grain structure for good dimensional stability under irradiation
      - Low porosity to minimize fuel salt and fission product infiltration

    The graphite was unclad -- it was in direct contact with the flowing fuel
    salt. Over time, xenon and krypton fission gases could diffuse into the
    open porosity of the graphite, affecting the xenon reactivity worth.
    ORNL developed special graphite sealing techniques to mitigate this.

    For nuclear modeling purposes, graphite is treated as pure C-12 (natural
    carbon). The S(alpha,beta) thermal scattering law for graphite is
    essential for correctly modeling the thermalization of neutrons, since
    the carbon atoms are bound in a crystalline lattice and their thermal
    motion is governed by the phonon spectrum rather than the free-gas model.
    """
    graphite = openmc.Material(name="CGB Graphite Moderator")

    # Natural carbon -- effectively 98.9% C-12 and 1.1% C-13, but we
    # use the C0 (natural carbon) designation in the cross section library.
    graphite.add_element("C", 1.0)

    # CGB-grade graphite density
    graphite.set_density("g/cm3", 1.86)

    # Add thermal scattering law for bound carbon in graphite.
    # This is critical for correct neutron thermalization. Without the
    # S(alpha,beta) treatment, the thermal spectrum would be wrong and
    # k-eff would be significantly off.
    graphite.add_s_alpha_beta("c_Graphite")

    # Same operating temperature as the fuel salt
    graphite.temperature = 922.0

    return graphite


def create_hastelloy_n():
    """Create the Hastelloy-N (INOR-8) reactor vessel material.

    Hastelloy-N (also called INOR-8 or Alloy N) was specifically developed
    at ORNL for molten salt reactor service. It is a nickel-molybdenum alloy
    with excellent corrosion resistance to fluoride salts at high temperature.

    Nominal composition (weight percent):
      Ni: balance (~71%)
      Mo: 16%
      Cr: 7%
      Fe: 4%
      Mn: 0.5%
      Si: 0.5%
      Others (C, W, Ti, etc.): ~1%

    For simplicity we use the major constituents and normalize.
    """
    hastelloy = openmc.Material(name="Hastelloy-N (INOR-8) Vessel")

    # Weight fractions of the major constituents
    hastelloy.add_element("Ni", 0.71, "wo")   # Nickel - base metal
    hastelloy.add_element("Mo", 0.16, "wo")   # Molybdenum - corrosion resistance
    hastelloy.add_element("Cr", 0.07, "wo")   # Chromium - oxidation resistance
    hastelloy.add_element("Fe", 0.04, "wo")   # Iron
    hastelloy.add_element("Mn", 0.01, "wo")   # Manganese
    hastelloy.add_element("Si", 0.01, "wo")   # Silicon

    # Hastelloy-N density
    hastelloy.set_density("g/cm3", 8.86)

    hastelloy.temperature = 922.0

    return hastelloy


def create_homogenized_core(fuel_salt, graphite):
    """Create a homogenized mixture of fuel salt and graphite for the core.

    In the MSRE core, the fuel salt occupies approximately 22.5% of the
    core volume and the graphite moderator occupies approximately 77.5%.
    These volume fractions come from the geometry of the graphite stringer
    array: the stringers are closely packed with narrow fuel channels
    between them.

    The homogenized material is created by mixing the fuel salt and graphite
    materials according to these volume fractions. This approach:
      - Preserves the overall atom densities of each component
      - Correctly represents the total macroscopic cross sections
      - Loses spatial self-shielding effects (heterogeneity effects)

    For the MSRE, the heterogeneity effect is relatively small because
    the fuel channels are narrow (~1 cm) and the thermal diffusion length
    in graphite (~59 cm) is much larger than the channel spacing. This
    means the thermal neutron flux is fairly uniform across the channel,
    and homogenization is a reasonable approximation.
    """
    # Volume fractions from the IRPhEP benchmark geometry
    fuel_volume_fraction = 0.225    # ~22.5% fuel salt by volume
    graphite_volume_fraction = 0.775  # ~77.5% graphite by volume

    # We cannot use openmc.Material.mix_materials because the graphite
    # contains an S(alpha,beta) thermal scattering table, which is not
    # supported by mix_materials. Instead, we manually compute the
    # homogenized atom densities.
    #
    # The effective number density of each nuclide is:
    #   N_eff = V_frac * N_original
    # where N_original is the atom density from the original material.
    # Since we specify atom fractions and a macroscopic density, we
    # scale the density contributions by volume fraction.

    core_mix = openmc.Material(name="Homogenized Core (Fuel Salt + Graphite)")

    # Compute an effective density as a volume-weighted average.
    # This is exact for the total atom count per unit volume.
    fuel_density = 2.32   # g/cm3
    graphite_density = 1.86  # g/cm3
    effective_density = (
        fuel_volume_fraction * fuel_density
        + graphite_volume_fraction * graphite_density
    )
    core_mix.set_density("g/cm3", effective_density)

    # Compute mass fractions for each nuclide.
    # Mass of fuel salt component per cm3 of mixture:
    fuel_mass = fuel_volume_fraction * fuel_density      # g/cm3
    graphite_mass = graphite_volume_fraction * graphite_density  # g/cm3
    total_mass = fuel_mass + graphite_mass  # = effective_density

    # Fuel salt nuclide atom fractions (from create_fuel_salt)
    # We add each nuclide scaled by the fuel salt mass fraction.
    fuel_wt = fuel_mass / total_mass  # weight fraction of fuel salt in mix

    # Get atom fractions from the fuel salt material and add them
    # weighted by the fuel mass fraction. We use "ao" (atom fraction)
    # and will normalize at the end.
    #
    # Actually, the simplest correct approach: add nuclides by their
    # partial densities (mass per unit volume). We'll use weight fractions
    # directly.

    # Fuel salt nuclides -- copy from the fuel salt definition above.
    # These are atom fractions within the fuel salt.
    fuel_nuclides = {
        "Li7":  0.26337 * 0.9999,
        "Li6":  0.26337 * 0.0001,
        "Be9":  0.11791,
        "Zr90": 0.02026 * 0.5145,
        "Zr91": 0.02026 * 0.1122,
        "Zr92": 0.02026 * 0.1715,
        "Zr94": 0.02026 * 0.1738,
        "Zr96": 0.02026 * 0.0280,
        "U235": 0.00365 * 0.33,
        "U238": 0.00365 * 0.67,
        "F19":  0.59481,
    }

    # Graphite: pure carbon (natural)
    # Graphite has 1 atom type: C (natural carbon, ~98.9% C-12, 1.1% C-13)
    # Atom fraction within graphite: 1.0

    # To combine, we need to convert atom fractions to a common basis.
    # Number density N = rho * N_A / M_avg, where M_avg is the average
    # atomic mass. We compute the average atomic mass for each material
    # and use that to get the atom density ratio.

    # Average atomic mass of fuel salt (from atom fractions and atomic masses)
    atomic_masses = {
        "Li7": 7.016, "Li6": 6.015, "Be9": 9.012, "F19": 18.998,
        "Zr90": 89.905, "Zr91": 90.906, "Zr92": 91.905,
        "Zr94": 93.906, "Zr96": 95.908,
        "U235": 235.044, "U238": 238.051, "C": 12.011,
    }

    fuel_M_avg = sum(
        frac * atomic_masses[nuc] for nuc, frac in fuel_nuclides.items()
    )
    graphite_M_avg = 12.011  # pure carbon

    # Number densities (relative, in units of N_A)
    # N = rho / M_avg  (atoms per cm3 in units of N_A/cm3)
    fuel_N = fuel_density / fuel_M_avg      # proportional to atom density
    graphite_N = graphite_density / graphite_M_avg

    # Effective number densities in the mixture
    fuel_N_eff = fuel_volume_fraction * fuel_N
    graphite_N_eff = graphite_volume_fraction * graphite_N

    total_N_eff = fuel_N_eff + graphite_N_eff

    # Add fuel salt nuclides with effective atom fractions
    for nuc, frac in fuel_nuclides.items():
        core_mix.add_nuclide(nuc, frac * fuel_N_eff / total_N_eff)

    # Add graphite carbon with effective atom fraction
    core_mix.add_element("C", 1.0 * graphite_N_eff / total_N_eff)

    # Add S(alpha,beta) thermal scattering for the graphite component.
    # Even in the homogenized mixture, the carbon atoms are still bound
    # in the graphite crystal lattice, so the graphite scattering kernel
    # applies to the carbon fraction.
    core_mix.add_s_alpha_beta("c_Graphite")

    core_mix.temperature = 922.0

    return core_mix


# =============================================================================
# Geometry definitions
# =============================================================================

def build_homogeneous_model(fuel_salt, graphite, hastelloy):
    """Build the homogeneous (fully homogenized) cylindrical core model.

    This model represents the MSRE as three concentric cylindrical regions:
      1. Core region: homogenized fuel salt + graphite mixture
      2. Reactor vessel wall: Hastelloy-N
      3. External void (vacuum boundary)

    Key dimensions from the IRPhEP benchmark MSRE-MSR-EXP-001:
      - Core active region radius:   ~68.58 cm (27 inches)
        This is the radius of the graphite-moderated lattice region.
      - Reactor vessel inner radius: ~73.66 cm (29 inches)
        The annular gap between core and vessel is the downcomer,
        which in this simplified model we fill with fuel salt.
      - Reactor vessel wall thickness: ~2.54 cm (1 inch)
      - Reactor vessel outer radius:  ~76.20 cm (30 inches)
      - Core active height:           ~162.56 cm (64 inches)
        This is the height of the graphite-moderated region.
      - Vessel total height:          ~238.76 cm (94 inches)
        Includes upper and lower plenums.

    For this simplified model, we use:
      - Homogenized core filling the vessel interior
      - Hastelloy-N vessel wall
      - Vacuum boundary conditions
    """
    # --- Create the homogenized core material ---
    core_material = create_homogenized_core(fuel_salt, graphite)

    materials = openmc.Materials([fuel_salt, graphite, hastelloy, core_material])

    # -------------------------------------------------------------------------
    # Define bounding surfaces
    # -------------------------------------------------------------------------

    # Core region (active lattice of graphite stringers + fuel channels)
    core_radius = openmc.ZCylinder(r=68.58, name="Core outer radius")

    # Downcomer annulus filled with fuel salt
    vessel_inner_radius = openmc.ZCylinder(r=73.66, name="Vessel inner wall")

    # Reactor vessel outer surface
    vessel_outer_radius = openmc.ZCylinder(r=76.20, name="Vessel outer wall")

    # Axial boundaries
    # The core active region sits within the vessel, with plenums above and below.
    # Lower plenum bottom (approximate vessel bottom)
    vessel_bottom = openmc.ZPlane(z0=-119.38, name="Vessel bottom")
    # Core bottom (lower grid plate elevation)
    core_bottom = openmc.ZPlane(z0=-81.28, name="Core bottom")
    # Core top (upper grid plate elevation)
    core_top = openmc.ZPlane(z0=81.28, name="Core top")
    # Vessel top
    vessel_top = openmc.ZPlane(z0=119.38, name="Vessel top")

    # -------------------------------------------------------------------------
    # Define cells and regions
    # -------------------------------------------------------------------------

    # Homogenized core region (graphite + fuel salt mixture)
    # This is where the fission chain reaction is sustained.
    core_cell = openmc.Cell(name="Homogenized Core")
    core_cell.region = -core_radius & +core_bottom & -core_top
    core_cell.fill = core_material

    # Downcomer region (annular gap between core and vessel wall)
    # In the real MSRE, fuel salt flows downward through this annulus
    # before entering the bottom plenum and flowing up through the core.
    downcomer_cell = openmc.Cell(name="Downcomer (fuel salt)")
    downcomer_cell.region = +core_radius & -vessel_inner_radius & +core_bottom & -core_top
    downcomer_cell.fill = fuel_salt

    # Upper plenum (fuel salt above the core)
    # The hot fuel salt collects here before exiting through the outlet nozzle.
    upper_plenum = openmc.Cell(name="Upper Plenum (fuel salt)")
    upper_plenum.region = -vessel_inner_radius & +core_top & -vessel_top
    upper_plenum.fill = fuel_salt

    # Lower plenum (fuel salt below the core)
    # The returning cooled fuel salt enters here from the downcomer.
    lower_plenum = openmc.Cell(name="Lower Plenum (fuel salt)")
    lower_plenum.region = -vessel_inner_radius & -core_bottom & +vessel_bottom
    lower_plenum.fill = fuel_salt

    # Reactor vessel wall (Hastelloy-N)
    # The vessel contains the fuel salt and provides the primary pressure boundary.
    vessel_wall = openmc.Cell(name="Reactor Vessel Wall (Hastelloy-N)")
    vessel_wall.region = (
        +vessel_inner_radius & -vessel_outer_radius & +vessel_bottom & -vessel_top
    )
    vessel_wall.fill = hastelloy

    # Vessel top and bottom caps (Hastelloy-N)
    vessel_top_cap = openmc.Cell(name="Vessel Top Cap")
    vessel_top_cap.region = -vessel_outer_radius & +vessel_top & -openmc.ZPlane(z0=121.92)
    vessel_top_cap.fill = hastelloy

    vessel_bottom_cap = openmc.Cell(name="Vessel Bottom Cap")
    vessel_bottom_cap.region = (
        -vessel_outer_radius
        & -vessel_bottom
        & +openmc.ZPlane(z0=-121.92)
    )
    vessel_bottom_cap.fill = hastelloy

    # External void -- everything outside the vessel is vacuum.
    # We set vacuum boundary conditions on a bounding box.
    outer_bound = openmc.model.RectangularParallelepiped(
        -90, 90, -90, 90, -135, 135,
        boundary_type="vacuum",
    )
    # The void fills all space inside the bounding box that is NOT occupied
    # by any vessel component. This includes:
    #   - The annular region outside the vessel (radially)
    #   - The space above and below the vessel caps
    cap_top = openmc.ZPlane(z0=121.92)
    cap_bottom = openmc.ZPlane(z0=-121.92)
    void_cell = openmc.Cell(name="External Void")
    void_cell.region = -outer_bound & (
        +vessel_outer_radius |  # outside vessel radially (at vessel height)
        +cap_top |              # above the top cap
        -cap_bottom             # below the bottom cap
    )

    # Build the universe with all cells including the void
    root_universe = openmc.Universe(
        cells=[
            core_cell,
            downcomer_cell,
            upper_plenum,
            lower_plenum,
            vessel_wall,
            vessel_top_cap,
            vessel_bottom_cap,
            void_cell,
        ]
    )

    geometry = openmc.Geometry(root_universe)

    return materials, geometry


def build_channel_model(fuel_salt, graphite):
    """Build a single fuel channel unit cell with reflective boundary conditions.

    This model represents one repeating unit of the MSRE graphite stringer
    lattice. The MSRE core contained ~1140 graphite stringers arranged in a
    regular pattern, with fuel salt flowing in the channels between them.

    Each stringer had a roughly rectangular cross section. The fuel channels
    between adjacent stringers were approximately 1.016 cm (0.4 inches) wide.
    The stringer pitch (center-to-center distance) was approximately
    5.08 cm (2 inches).

    For the unit cell model we use:
      - A square graphite block representing one stringer cell
      - A central square fuel channel carved out of the graphite
      - Reflective boundary conditions on all six faces to simulate
        an infinite lattice (no leakage)

    This captures the heterogeneous self-shielding effect: the spatial
    separation between fuel and moderator means thermal neutrons are
    preferentially born in the graphite (where they are moderated) and
    then stream into the fuel channel to cause fission. This spatial
    variation in the thermal flux is called "spatial self-shielding"
    and generally increases k-inf compared to a homogenized model.

    Unit cell dimensions:
      - Stringer pitch:     5.08 cm x 5.08 cm square
      - Fuel channel width: ~1.016 cm (centered in the cell)
      - Channel height:     162.56 cm (same as core active height)

    With reflective BCs on all faces, this is an eigenvalue problem for
    k-infinity (no leakage). The result will be higher than k-eff for
    the full core because it ignores radial and axial neutron leakage.
    """
    materials = openmc.Materials([fuel_salt, graphite])

    # -------------------------------------------------------------------------
    # Unit cell dimensions
    # -------------------------------------------------------------------------

    # Half-pitch of the square graphite stringer cell
    half_pitch = 5.08 / 2.0  # 2.54 cm

    # Half-width of the square fuel channel
    # The fuel channel gap between stringers is approximately 1.016 cm.
    # In this idealized model, we represent it as a central square channel.
    half_channel = 1.016 / 2.0  # 0.508 cm

    # Active core height
    half_height = 162.56 / 2.0  # 81.28 cm

    # -------------------------------------------------------------------------
    # Define surfaces
    # -------------------------------------------------------------------------

    # Outer boundary of the unit cell (reflective BCs to simulate infinite lattice)
    cell_boundary = openmc.model.RectangularParallelepiped(
        -half_pitch, half_pitch,
        -half_pitch, half_pitch,
        -half_height, half_height,
        boundary_type="reflective",
    )

    # Fuel channel boundaries (a square channel in the center of the cell)
    channel_xmin = openmc.XPlane(x0=-half_channel)
    channel_xmax = openmc.XPlane(x0=+half_channel)
    channel_ymin = openmc.YPlane(y0=-half_channel)
    channel_ymax = openmc.YPlane(y0=+half_channel)

    # The fuel channel region is the intersection of the four half-spaces
    fuel_channel_region = +channel_xmin & -channel_xmax & +channel_ymin & -channel_ymax

    # -------------------------------------------------------------------------
    # Define cells
    # -------------------------------------------------------------------------

    # Fuel salt channel -- this is where the molten salt flows upward
    # through the core, carrying the dissolved uranium fuel.
    fuel_cell = openmc.Cell(name="Fuel Salt Channel")
    fuel_cell.region = fuel_channel_region & cell_boundary.region
    fuel_cell.fill = fuel_salt

    # Graphite moderator block surrounding the fuel channel.
    # The graphite thermalizes fast neutrons via elastic scattering
    # from carbon-12 nuclei. It takes on average ~115 collisions to
    # thermalize a 2 MeV fission neutron in graphite (compared to ~18
    # in water), but graphite's very low absorption cross section means
    # neutrons survive many more collisions.
    graphite_cell = openmc.Cell(name="Graphite Moderator")
    graphite_cell.region = ~fuel_channel_region & cell_boundary.region
    graphite_cell.fill = graphite

    geometry = openmc.Geometry([fuel_cell, graphite_cell])

    return materials, geometry


# =============================================================================
# Settings
# =============================================================================

def create_settings(args, geometry_type):
    """Create OpenMC simulation settings.

    For the eigenvalue (k-eff) calculation, we need to specify:
      - Number of particles per batch (neutron histories)
      - Number of total batches
      - Number of inactive batches (for fission source convergence)
      - Initial source distribution

    The inactive batches allow the fission source distribution to converge
    from the initial guess to the true fundamental-mode distribution before
    we begin accumulating statistics. For this reactor, 50 inactive batches
    is generally sufficient.
    """
    settings = openmc.Settings()

    # Monte Carlo parameters
    settings.particles = args.particles
    settings.batches = args.batches
    settings.inactive = args.inactive

    # Initial neutron source distribution.
    # We start neutrons uniformly distributed in a box covering the fuel
    # region. The source will converge during the inactive batches.
    if geometry_type == "homogeneous":
        # Uniform source throughout the cylindrical core
        # Use a cylindrical spatial distribution that covers the core
        lower_left = (-68.0, -68.0, -80.0)
        upper_right = (68.0, 68.0, 80.0)
    else:
        # For the channel model, source in the fuel channel
        lower_left = (-0.5, -0.5, -80.0)
        upper_right = (0.5, 0.5, 80.0)

    uniform_dist = openmc.stats.Box(lower_left, upper_right)
    settings.source = openmc.IndependentSource(space=uniform_dist)

    # Temperature treatment
    # Use nearest-temperature method to handle the operating temperature.
    # The cross section libraries may not have data at exactly 922 K,
    # so OpenMC will interpolate or use the nearest available temperature.
    settings.temperature = {"method": "nearest", "tolerance": 1000}

    return settings


# =============================================================================
# Main execution
# =============================================================================

def main():
    """Build and export the MSRE OpenMC model.

    This function:
      1. Parses command-line arguments
      2. Creates all materials (fuel salt, graphite, Hastelloy-N)
      3. Builds the selected geometry model
      4. Configures simulation settings
      5. Exports the model to XML files for OpenMC execution

    After running this script, execute the model with:
      openmc
    or
      openmc --threads <N>
    """
    args = parse_args()

    print("=" * 70)
    print("MSRE - Molten Salt Reactor Experiment")
    print("OpenMC Benchmark Model Generator")
    print("=" * 70)
    print()

    # --- Create materials ---
    # These material definitions are shared between both model types.
    print("Creating materials...")
    print("  - Fuel salt: LiF-BeF2-ZrF4-UF4 (65-29.1-5-0.9 mol%)")
    print("    Li-7 enrichment: 99.99 at%")
    print("    U-235 enrichment: 33 at%")
    print("    Density: 2.32 g/cm3 at 922 K")
    print("  - Graphite moderator: CGB grade, 1.86 g/cm3")

    fuel_salt = create_fuel_salt()
    graphite = create_graphite()
    hastelloy = create_hastelloy_n()

    # --- Build geometry based on selected model type ---
    print()
    print(f"Building geometry: {args.model} model")

    if args.model == "homogeneous":
        print("  Homogenized core: 22.5% fuel salt, 77.5% graphite by volume")
        print("  Core radius: 68.58 cm")
        print("  Vessel inner radius: 73.66 cm")
        print("  Core active height: 162.56 cm")
        materials, geometry = build_homogeneous_model(fuel_salt, graphite, hastelloy)
    else:
        print("  Single fuel channel unit cell with reflective BCs")
        print("  Stringer pitch: 5.08 cm x 5.08 cm")
        print("  Fuel channel width: 1.016 cm")
        print("  Channel height: 162.56 cm")
        materials, geometry = build_channel_model(fuel_salt, graphite)

    # --- Configure settings ---
    print()
    print(f"Simulation settings:")
    print(f"  Particles per batch: {args.particles}")
    print(f"  Total batches:       {args.batches}")
    print(f"  Inactive batches:    {args.inactive}")

    settings = create_settings(args, args.model)

    # --- Export the model ---
    # OpenMC reads XML input files: geometry.xml, materials.xml, settings.xml
    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
    model.export_to_model_xml()

    print()
    print("Model exported to: model.xml")
    print()
    print("To run the simulation:")
    print("  openmc")
    print("  openmc --threads <N>  (for parallel execution)")
    print()
    print("Reference k-eff values:")
    print("  Experimental:          0.99978 +/- 0.00420")
    print("  OpenMC CSG (detailed): 1.00878 +/- 0.00032")
    print("  Serpent:               1.02132 +/- 0.00003")
    print()

    if args.model == "homogeneous":
        print("NOTE: The homogeneous model will differ from the detailed CSG")
        print("result because it removes spatial self-shielding effects and")
        print("simplifies the geometry (no control rods, simplified vessel).")
    else:
        print("NOTE: The channel model computes k-infinity (infinite lattice),")
        print("which will be HIGHER than the full-core k-eff because there is")
        print("no neutron leakage with reflective boundary conditions.")

    print()


if __name__ == "__main__":
    main()
