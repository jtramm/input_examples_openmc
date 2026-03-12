#!/usr/bin/env python3
"""
IAEA 10 MW MTR (Material Testing Reactor) Research Reactor Benchmark
=====================================================================

This script builds an OpenMC model for the IAEA 10 MW MTR benchmark as specified
in IAEA-TECDOC-643, "Research Reactor Core Conversion Guidebook" (1992).

BACKGROUND
----------
The IAEA 10 MW MTR benchmark is one of the most widely used research reactor
benchmarks in the world. It was developed under the Reduced Enrichment for
Research and Test Reactors (RERTR) program to provide a common calculational
benchmark for comparing neutronic analysis methods across international
institutions. Dozens of organizations have analyzed this benchmark using
diffusion theory, transport theory, and Monte Carlo methods.

Unlike power reactor benchmarks that use cylindrical fuel pins, the IAEA MTR
uses PLATE-TYPE fuel elements. Each standard fuel element (SFE) contains 23
thin fuel plates arranged in a rectangular housing. The fuel meat is a
uranium-aluminum (UAl) alloy clad in aluminum, with light water flowing
between the plates as both moderator and coolant.

GEOMETRY
--------
This is a 2D (XY) model with reflective boundary conditions on the top and
bottom (Z-direction), which effectively models an infinite-height core. This
is consistent with the standard 2D benchmark problem and avoids axial leakage
complications while preserving the essential physics of plate-type fuel.

The core consists of a 5-column x 6-row rectangular arrangement of fuel elements:
  - 23 Standard Fuel Elements (SFE): each with 23 fuel plates
  - 5 Control Fuel Elements (CFE): 17 plates each with wider water channel
    (control absorber withdrawn for the fresh-core eigenvalue calculation)
  - 2 Water irradiation positions at diagonally opposite corners
  - Beryllium reflector surrounding the core
  - Graphite reflector outside the beryllium

REFERENCE VALUES
----------------
Published k-eff values for the fresh HEU 3D core (all rods out) range from
approximately 1.05 to 1.08 depending on the calculational method and nuclear
data library used.

IMPORTANT: This 2D model uses reflective Z-boundaries, which eliminates axial
leakage entirely. The 3D core has only 60 cm active height, so axial leakage
is very significant. Therefore, the 2D k-eff will be substantially higher
than the 3D reference values (typically k_2D ~ 1.4-1.5 vs k_3D ~ 1.06).
This is expected and physically correct behavior.

DIMENSIONS (from IAEA-TECDOC-643 and widely published companion papers)
------------------------------------------------------------------------
  Fuel meat thickness:     0.051 cm   (0.51 mm)
  Clad thickness:          0.038 cm   (0.38 mm) per side
  Total plate thickness:   0.127 cm   (= 0.051 + 2*0.038)
  Water channel gap:       0.223 cm   (2.23 mm between plates)
  Plate-to-plate pitch:    0.350 cm   (= 0.127 + 0.223)
  Active fuel meat width:  6.3 cm
  Fuel plate width:        6.65 cm    (meat + side plates)
  Element cross-section:   7.7 cm x 8.1 cm (outer aluminum frame)
  Active fuel height:      60.0 cm    (modeled as infinite in 2D)
  U-235 per SFE:           280 g

Core layout (columns A-E, rows 1-6):
  Row 1:  H2O  SFE  SFE  SFE  SFE
  Row 2:  SFE  CFE  SFE  CFE  SFE
  Row 3:  SFE  SFE  SFE  SFE  SFE
  Row 4:  SFE  SFE  CFE  SFE  SFE
  Row 5:  SFE  CFE  SFE  CFE  SFE
  Row 6:  SFE  SFE  SFE  SFE  H2O
  (H2O = irradiation/water position; 23 SFE + 5 CFE + 2 H2O = 30)
"""

import argparse
import numpy as np
import openmc
import openmc.model


# =============================================================================
# Command-line arguments
# =============================================================================

parser = argparse.ArgumentParser(
    description="IAEA 10 MW MTR Benchmark - OpenMC Model"
)
parser.add_argument(
    "--particles", type=int, default=20000,
    help="Number of particles per batch (default: 20000)"
)
parser.add_argument(
    "--batches", type=int, default=200,
    help="Total number of batches (default: 200)"
)
parser.add_argument(
    "--inactive", type=int, default=50,
    help="Number of inactive batches (default: 50)"
)
parser.add_argument(
    "--small-tallies", action="store_true",
    help="Use smaller tally bins for faster runtime"
)
args = parser.parse_args()


# =============================================================================
# Geometric parameters (all in cm)
# =============================================================================

# --- Fuel plate dimensions ---
FUEL_MEAT_THICKNESS = 0.051       # UAl alloy fuel meat thickness
CLAD_THICKNESS = 0.038            # Aluminum clad thickness (each side)
PLATE_THICKNESS = 0.127           # Total plate thickness (meat + 2*clad)
WATER_GAP = 0.223                 # Water channel between adjacent plates
PLATE_PITCH = 0.350               # Plate-to-plate pitch (plate + water gap)

# --- Fuel plate width ---
FUEL_MEAT_WIDTH = 6.3             # Width of the active fuel meat region
PLATE_WIDTH = 6.65                # Total fuel plate width including edges

# --- Fuel element outer dimensions ---
ELEMENT_WIDTH_X = 7.7             # Element width in X (across plates)
ELEMENT_WIDTH_Y = 8.1             # Element width in Y (along plates)

# --- Number of plates ---
NUM_PLATES_SFE = 23               # Plates in a Standard Fuel Element
NUM_PLATES_CFE = 17               # Plates in a Control Fuel Element (fewer plates)
# CFE has a wider central water channel where control absorber sits (withdrawn here)

# --- Active height (used for 3D; infinite in our 2D model) ---
ACTIVE_HEIGHT = 60.0              # Active fuel height in cm

# --- Reflector thicknesses ---
BE_REFLECTOR_THICKNESS = 7.7      # Beryllium reflector thickness (one element width)
GRAPHITE_REFLECTOR_THICKNESS = 10.0  # Outer graphite reflector thickness

# --- Half-dimensions for the model (Z is arbitrary for 2D) ---
HALF_HEIGHT = 30.0                # Half of active height for Z bounds


# =============================================================================
# Materials
# =============================================================================

# --- Fuel meat: Uranium-Aluminum alloy (UAl) ---
# HEU fuel: 93 wt% U-235 enriched uranium in aluminum matrix
# The fuel meat density is determined by the U-235 loading requirement.
# For 280 g U-235 per SFE with 23 plates:
#   U-235 per plate = 280/23 = 12.17 g
#   Volume per plate = 6.3 * 0.051 * 60 = 19.278 cm^3
#   U density in meat ~ 0.68 g/cm^3 (U-235) -> total U ~ 0.731 g/cm^3
# The UAl alloy meat has a composite density accounting for both U and Al.
# A commonly used fuel meat density is ~3.44 g/cm^3 for HEU UAl alloy.

fuel = openmc.Material(name="UAl fuel meat (HEU 93%)")
fuel.set_density("g/cm3", 3.44)
# Uranium fraction: U-235 loading determines the uranium weight fraction
# 280 g U-235 / 23 plates = 12.17 g U-235 per plate
# Plate meat volume = 6.3 * 0.051 * 60 = 19.278 cm^3
# Mass of meat per plate = 19.278 * 3.44 = 66.32 g
# U-235 mass fraction = 12.17 / 66.32 = 0.1835
# Total U mass fraction = 0.1835 / 0.93 = 0.1973 (since 93% enriched)
# Al mass fraction = 1 - 0.1973 = 0.8027
fuel.add_nuclide("U235", 0.1835, "wo")      # 93% of uranium is U-235
fuel.add_nuclide("U238", 0.0138, "wo")      # 7% of uranium is U-238
fuel.add_element("Al", 0.8027, "wo")         # Aluminum matrix balance
fuel.temperature = 293.6                      # Room temperature in K

# --- Cladding: Aluminum alloy (similar to Al-6061) ---
# For simplicity, modeled as pure aluminum. The small alloying elements
# (Mg, Si, etc.) have negligible neutronic effect.
clad = openmc.Material(name="Aluminum cladding")
clad.set_density("g/cm3", 2.70)              # Standard Al density
clad.add_element("Al", 1.0)                  # Pure aluminum approximation
clad.temperature = 293.6

# --- Moderator/Coolant: Light water ---
# Room temperature water at atmospheric pressure
water = openmc.Material(name="Light water (moderator/coolant)")
water.set_density("g/cm3", 0.997)            # Water density at ~20 C
water.add_nuclide("H1", 2.0, "ao")          # Two hydrogen atoms
water.add_nuclide("O16", 1.0, "ao")         # One oxygen atom
water.add_s_alpha_beta("c_H_in_H2O")        # Thermal scattering for bound H
water.temperature = 293.6

# --- Beryllium reflector ---
# Solid beryllium metal used as the inner reflector
beryllium = openmc.Material(name="Beryllium reflector")
beryllium.set_density("g/cm3", 1.85)         # Be metal density
beryllium.add_element("Be", 1.0)             # Pure beryllium
beryllium.add_s_alpha_beta("c_Be")            # Thermal scattering for Be metal
beryllium.temperature = 293.6

# --- Graphite reflector ---
# Reactor-grade graphite used as the outer reflector
graphite = openmc.Material(name="Graphite reflector")
graphite.set_density("g/cm3", 1.70)           # Reactor-grade graphite density
graphite.add_element("C", 1.0)                # Pure carbon
graphite.add_s_alpha_beta("c_Graphite")       # Thermal scattering for graphite
graphite.temperature = 293.6

# Collect all materials into a Materials object
materials = openmc.Materials([fuel, clad, water, beryllium, graphite])


# =============================================================================
# Helper function: build a single fuel plate (fuel meat + clad)
# =============================================================================

def make_fuel_plate(x_center, plate_id_offset=0):
    """
    Create regions for a single fuel plate centered at x_center.

    A fuel plate consists of three layers (in X):
      - Left clad  (aluminum)
      - Fuel meat   (UAl alloy)
      - Right clad  (aluminum)

    The plate extends across the full Y-width and Z-height of the element.

    Parameters
    ----------
    x_center : float
        X-coordinate of the plate center
    plate_id_offset : int
        Offset for surface IDs to ensure uniqueness

    Returns
    -------
    list of (region, material) tuples
    """
    # Half-thicknesses
    half_plate = PLATE_THICKNESS / 2.0       # Half of total plate
    half_meat = FUEL_MEAT_THICKNESS / 2.0    # Half of fuel meat

    # Surfaces bounding the fuel meat (inner layer)
    meat_left = openmc.XPlane(x0=x_center - half_meat)
    meat_right = openmc.XPlane(x0=x_center + half_meat)

    # Surfaces bounding the full plate (outer aluminum clad)
    plate_left = openmc.XPlane(x0=x_center - half_plate)
    plate_right = openmc.XPlane(x0=x_center + half_plate)

    # Build the three regions:
    # 1) Left clad: between plate_left and meat_left
    left_clad_region = +plate_left & -meat_left
    # 2) Fuel meat: between meat_left and meat_right
    meat_region = +meat_left & -meat_right
    # 3) Right clad: between meat_right and plate_right
    right_clad_region = +meat_right & -plate_right

    return [
        (left_clad_region, clad),
        (meat_region, fuel),
        (right_clad_region, clad),
    ]


# =============================================================================
# Build a Standard Fuel Element (SFE) universe
# =============================================================================

def build_sfe_universe():
    """
    Build a Standard Fuel Element (SFE) containing 23 fuel plates.

    The SFE has 23 equally-spaced fuel plates arranged in the X-direction,
    with water channels between them. The entire assembly is contained within
    an outer aluminum frame (side plates).

    The element is centered at the origin of its local universe.

    Returns
    -------
    openmc.Universe
        Universe containing the SFE geometry
    """
    cells = []

    # --- Calculate plate positions ---
    # 23 plates with pitch of 0.350 cm, centered at x=0
    # Total span of plates: 22 * 0.350 = 7.70 cm (between first and last centers)
    # Centering: first plate at x = -22/2 * 0.350 = -3.85 cm
    first_plate_x = -(NUM_PLATES_SFE - 1) / 2.0 * PLATE_PITCH

    # Collect all plate surfaces for building water regions
    all_plate_surfaces = []  # List of (left_surface, right_surface) per plate

    for i in range(NUM_PLATES_SFE):
        x_center = first_plate_x + i * PLATE_PITCH

        # Half-thicknesses
        half_plate = PLATE_THICKNESS / 2.0

        # Create plate surfaces
        plate_left = openmc.XPlane(x0=x_center - half_plate)
        plate_right = openmc.XPlane(x0=x_center + half_plate)
        all_plate_surfaces.append((plate_left, plate_right))

        # Fuel meat surfaces
        half_meat = FUEL_MEAT_THICKNESS / 2.0
        meat_left = openmc.XPlane(x0=x_center - half_meat)
        meat_right = openmc.XPlane(x0=x_center + half_meat)

        # Left clad cell
        left_clad_cell = openmc.Cell(
            name=f"SFE plate {i+1} left clad",
            fill=clad,
            region=+plate_left & -meat_left
        )
        cells.append(left_clad_cell)

        # Fuel meat cell
        meat_cell = openmc.Cell(
            name=f"SFE plate {i+1} fuel meat",
            fill=fuel,
            region=+meat_left & -meat_right
        )
        cells.append(meat_cell)

        # Right clad cell
        right_clad_cell = openmc.Cell(
            name=f"SFE plate {i+1} right clad",
            fill=clad,
            region=+meat_right & -plate_right
        )
        cells.append(right_clad_cell)

    # --- Water channels between plates ---
    for i in range(NUM_PLATES_SFE - 1):
        # Water between plate i (right surface) and plate i+1 (left surface)
        water_cell = openmc.Cell(
            name=f"SFE water channel {i+1}",
            fill=water,
            region=+all_plate_surfaces[i][1] & -all_plate_surfaces[i+1][0]
        )
        cells.append(water_cell)

    # --- Outer element boundary ---
    # The element has an outer aluminum frame. Water fills the space between
    # the outermost plates and the frame.
    half_elem_x = ELEMENT_WIDTH_X / 2.0
    half_elem_y = ELEMENT_WIDTH_Y / 2.0

    elem_left = openmc.XPlane(x0=-half_elem_x)
    elem_right = openmc.XPlane(x0=half_elem_x)
    elem_front = openmc.YPlane(y0=-half_elem_y)
    elem_back = openmc.YPlane(y0=half_elem_y)

    # Water region to the left of the first plate (between element wall and first plate)
    water_left = openmc.Cell(
        name="SFE water left margin",
        fill=water,
        region=+elem_left & -all_plate_surfaces[0][0] & +elem_front & -elem_back
    )
    cells.append(water_left)

    # Water region to the right of the last plate
    water_right = openmc.Cell(
        name="SFE water right margin",
        fill=water,
        region=+all_plate_surfaces[-1][1] & -elem_right & +elem_front & -elem_back
    )
    cells.append(water_right)

    # --- Side plate regions (Y boundaries) ---
    # The fuel plates have a limited width in Y. The side regions (in Y beyond
    # the plate width) are aluminum side plates + water.
    # For simplicity, we fill the full Y extent within the element boundary
    # and the plates span the full Y range. The side plates are thin Al strips
    # along the Y edges. We approximate by filling the element uniformly in Y
    # for the plate region, then adding Y-boundary water/Al.
    #
    # Since this is a simplified model and the plates extend nearly the full
    # Y width, we treat the inter-plate regions as extending the full Y width
    # of the element. This is a standard simplification in 2D MTR models.

    # Apply Y bounds to all inter-plate water cells (they need Y limits)
    for cell in cells:
        if "water channel" in cell.name:
            cell.region = cell.region & +elem_front & -elem_back
        elif "clad" in cell.name or "meat" in cell.name:
            cell.region = cell.region & +elem_front & -elem_back

    # Create the universe with the element boundary
    sfe_universe = openmc.Universe(name="Standard Fuel Element (SFE)")
    sfe_universe.add_cells(cells)

    return sfe_universe, (elem_left, elem_right, elem_front, elem_back)


# =============================================================================
# Build a Control Fuel Element (CFE) universe
# =============================================================================

def build_cfe_universe():
    """
    Build a Control Fuel Element (CFE) with 17 fuel plates.

    The CFE is similar to the SFE but has fewer fuel plates (17 instead of 23)
    to accommodate the control absorber mechanism. The central region where the
    absorber would sit is filled with water (representing the rod-withdrawn
    condition for the fresh core eigenvalue calculation).

    Returns
    -------
    openmc.Universe
        Universe containing the CFE geometry
    """
    cells = []

    # --- Calculate plate positions ---
    # 17 plates centered at x=0 with the same pitch as SFE
    first_plate_x = -(NUM_PLATES_CFE - 1) / 2.0 * PLATE_PITCH

    all_plate_surfaces = []

    for i in range(NUM_PLATES_CFE):
        x_center = first_plate_x + i * PLATE_PITCH

        half_plate = PLATE_THICKNESS / 2.0
        plate_left = openmc.XPlane(x0=x_center - half_plate)
        plate_right = openmc.XPlane(x0=x_center + half_plate)
        all_plate_surfaces.append((plate_left, plate_right))

        half_meat = FUEL_MEAT_THICKNESS / 2.0
        meat_left = openmc.XPlane(x0=x_center - half_meat)
        meat_right = openmc.XPlane(x0=x_center + half_meat)

        # Left clad
        left_clad_cell = openmc.Cell(
            name=f"CFE plate {i+1} left clad",
            fill=clad,
            region=+plate_left & -meat_left
        )
        cells.append(left_clad_cell)

        # Fuel meat
        meat_cell = openmc.Cell(
            name=f"CFE plate {i+1} fuel meat",
            fill=fuel,
            region=+meat_left & -meat_right
        )
        cells.append(meat_cell)

        # Right clad
        right_clad_cell = openmc.Cell(
            name=f"CFE plate {i+1} right clad",
            fill=clad,
            region=+meat_right & -plate_right
        )
        cells.append(right_clad_cell)

    # --- Water channels between plates ---
    for i in range(NUM_PLATES_CFE - 1):
        water_cell = openmc.Cell(
            name=f"CFE water channel {i+1}",
            fill=water,
            region=+all_plate_surfaces[i][1] & -all_plate_surfaces[i+1][0]
        )
        cells.append(water_cell)

    # --- Element boundary ---
    half_elem_x = ELEMENT_WIDTH_X / 2.0
    half_elem_y = ELEMENT_WIDTH_Y / 2.0

    elem_left = openmc.XPlane(x0=-half_elem_x)
    elem_right = openmc.XPlane(x0=half_elem_x)
    elem_front = openmc.YPlane(y0=-half_elem_y)
    elem_back = openmc.YPlane(y0=half_elem_y)

    # Water margins (wider than SFE since fewer plates)
    water_left = openmc.Cell(
        name="CFE water left margin",
        fill=water,
        region=+elem_left & -all_plate_surfaces[0][0] & +elem_front & -elem_back
    )
    cells.append(water_left)

    water_right = openmc.Cell(
        name="CFE water right margin",
        fill=water,
        region=+all_plate_surfaces[-1][1] & -elem_right & +elem_front & -elem_back
    )
    cells.append(water_right)

    # Apply Y bounds to all cells
    for cell in cells:
        if "water channel" in cell.name:
            cell.region = cell.region & +elem_front & -elem_back
        elif "clad" in cell.name or "meat" in cell.name:
            cell.region = cell.region & +elem_front & -elem_back

    cfe_universe = openmc.Universe(name="Control Fuel Element (CFE)")
    cfe_universe.add_cells(cells)

    return cfe_universe, (elem_left, elem_right, elem_front, elem_back)


# =============================================================================
# Build reflector element universes
# =============================================================================

def build_be_element_universe():
    """
    Build a beryllium reflector element (same size as a fuel element).

    In the MTR core layout, certain positions at the corners are occupied
    by beryllium reflector elements instead of fuel elements.

    Returns
    -------
    openmc.Universe
        Universe filled with beryllium
    """
    half_x = ELEMENT_WIDTH_X / 2.0
    half_y = ELEMENT_WIDTH_Y / 2.0

    left = openmc.XPlane(x0=-half_x)
    right = openmc.XPlane(x0=half_x)
    front = openmc.YPlane(y0=-half_y)
    back = openmc.YPlane(y0=half_y)

    # Single cell filled with beryllium
    be_cell = openmc.Cell(
        name="Beryllium reflector element",
        fill=beryllium,
        region=+left & -right & +front & -back
    )

    be_universe = openmc.Universe(name="Be reflector element")
    be_universe.add_cell(be_cell)

    return be_universe


def build_water_element_universe():
    """
    Build a water-filled element (irradiation position).

    In the IAEA MTR core layout, two corner positions are not occupied by
    fuel elements. These are irradiation positions filled with water.

    Returns
    -------
    openmc.Universe
        Universe filled with water
    """
    half_x = ELEMENT_WIDTH_X / 2.0
    half_y = ELEMENT_WIDTH_Y / 2.0

    left = openmc.XPlane(x0=-half_x)
    right = openmc.XPlane(x0=half_x)
    front = openmc.YPlane(y0=-half_y)
    back = openmc.YPlane(y0=half_y)

    # Single cell filled with water (irradiation position)
    water_cell = openmc.Cell(
        name="Water irradiation position",
        fill=water,
        region=+left & -right & +front & -back
    )

    water_universe = openmc.Universe(name="Water element (irradiation)")
    water_universe.add_cell(water_cell)

    return water_universe


# =============================================================================
# Build the universes
# =============================================================================

print("Building Standard Fuel Element (SFE) universe...")
sfe_univ, _ = build_sfe_universe()

print("Building Control Fuel Element (CFE) universe...")
cfe_univ, _ = build_cfe_universe()

print("Building Beryllium reflector element universe...")
be_univ = build_be_element_universe()

print("Building Water irradiation element universe...")
w_univ = build_water_element_universe()


# =============================================================================
# Core layout using a rectangular lattice
# =============================================================================
# The core is a 5x6 rectangular array of fuel elements.
# Using the layout from IAEA-TECDOC-643 (HEU core, all rods withdrawn):
#
# Row 6 (top):    H2O  SFE  SFE  SFE  SFE     (1 water irradiation position)
# Row 5:          SFE  CFE  SFE  CFE  SFE
# Row 4:          SFE  SFE  CFE  SFE  SFE
# Row 3:          SFE  SFE  SFE  SFE  SFE
# Row 2:          SFE  CFE  SFE  CFE  SFE
# Row 1 (bottom): SFE  SFE  SFE  SFE  H2O     (1 water irradiation position)
#
# Total: 23 SFE + 5 CFE + 2 H2O = 30 positions
# The two water positions are at diagonally opposite corners (irradiation slots).
#
# Note: OpenMC lattice indexing has the first row at the top (positive Y).

print("Building core lattice...")

# Create the rectangular lattice for the 5x6 core
core_lattice = openmc.RectLattice(name="MTR Core Lattice (5x6)")
core_lattice.lower_left = [
    -5 * ELEMENT_WIDTH_X / 2.0,   # Left edge of 5-column core
    -6 * ELEMENT_WIDTH_Y / 2.0,   # Bottom edge of 6-row core
]
core_lattice.pitch = [ELEMENT_WIDTH_X, ELEMENT_WIDTH_Y]  # Element pitch

# Shorthand aliases for readability
S = sfe_univ   # Standard Fuel Element
C = cfe_univ   # Control Fuel Element
W = w_univ     # Water irradiation position

# Define the lattice universes (top row first in OpenMC convention)
# Row ordering: from top (positive Y) to bottom (negative Y)
core_lattice.universes = [
    [W, S, S, S, S],   # Row 6 (top,    positive Y) - H2O at top-left corner
    [S, C, S, C, S],   # Row 5
    [S, S, C, S, S],   # Row 4
    [S, S, S, S, S],   # Row 3
    [S, C, S, C, S],   # Row 2
    [S, S, S, S, W],   # Row 1 (bottom, negative Y) - H2O at bottom-right corner
]

# Tally: count SFE and CFE positions
num_sfe = sum(list(row).count(S) for row in core_lattice.universes)
num_cfe = sum(list(row).count(C) for row in core_lattice.universes)
num_w = sum(list(row).count(W) for row in core_lattice.universes)
print(f"  Core layout: {num_sfe} SFE + {num_cfe} CFE + {num_w} H2O = {num_sfe+num_cfe+num_w} positions")


# =============================================================================
# Root geometry: core + reflectors
# =============================================================================

print("Building root geometry with reflectors...")

# --- Core boundary ---
# The core lattice spans 5*7.7 = 38.5 cm in X and 6*8.1 = 48.6 cm in Y
core_half_x = 5 * ELEMENT_WIDTH_X / 2.0   # 19.25 cm
core_half_y = 6 * ELEMENT_WIDTH_Y / 2.0   # 24.30 cm

core_left = openmc.XPlane(x0=-core_half_x, name="Core left boundary")
core_right = openmc.XPlane(x0=core_half_x, name="Core right boundary")
core_bottom = openmc.YPlane(y0=-core_half_y, name="Core bottom boundary")
core_top = openmc.YPlane(y0=core_half_y, name="Core top boundary")

# --- Beryllium reflector boundary (surrounds the core) ---
be_half_x = core_half_x + BE_REFLECTOR_THICKNESS
be_half_y = core_half_y + BE_REFLECTOR_THICKNESS

be_left = openmc.XPlane(x0=-be_half_x, name="Be reflector left")
be_right = openmc.XPlane(x0=be_half_x, name="Be reflector right")
be_bottom = openmc.YPlane(y0=-be_half_y, name="Be reflector bottom")
be_top = openmc.YPlane(y0=be_half_y, name="Be reflector top")

# --- Graphite reflector boundary (surrounds the beryllium) ---
gr_half_x = be_half_x + GRAPHITE_REFLECTOR_THICKNESS
gr_half_y = be_half_y + GRAPHITE_REFLECTOR_THICKNESS

gr_left = openmc.XPlane(x0=-gr_half_x, name="Graphite reflector left", boundary_type="vacuum")
gr_right = openmc.XPlane(x0=gr_half_x, name="Graphite reflector right", boundary_type="vacuum")
gr_bottom = openmc.YPlane(y0=-gr_half_y, name="Graphite reflector bottom", boundary_type="vacuum")
gr_top = openmc.YPlane(y0=gr_half_y, name="Graphite reflector top", boundary_type="vacuum")

# --- Z boundaries (reflective for 2D model) ---
z_bottom = openmc.ZPlane(z0=-HALF_HEIGHT, name="Z bottom (reflective)", boundary_type="reflective")
z_top = openmc.ZPlane(z0=HALF_HEIGHT, name="Z top (reflective)", boundary_type="reflective")

# --- Core cell: contains the lattice ---
core_region = +core_left & -core_right & +core_bottom & -core_top & +z_bottom & -z_top
core_cell = openmc.Cell(
    name="Core (5x6 fuel element lattice)",
    fill=core_lattice,
    region=core_region
)

# --- Beryllium reflector cell: annular region between core and Be boundary ---
be_region = (
    +be_left & -be_right & +be_bottom & -be_top &  # Inside Be boundary
    ~(+core_left & -core_right & +core_bottom & -core_top) &  # Outside core
    +z_bottom & -z_top
)
be_cell = openmc.Cell(
    name="Beryllium reflector (surrounding core)",
    fill=beryllium,
    region=be_region
)

# --- Graphite reflector cell: annular region between Be and graphite boundary ---
gr_region = (
    +gr_left & -gr_right & +gr_bottom & -gr_top &  # Inside graphite boundary
    ~(+be_left & -be_right & +be_bottom & -be_top) &  # Outside Be boundary
    +z_bottom & -z_top
)
gr_cell = openmc.Cell(
    name="Graphite reflector (outer)",
    fill=graphite,
    region=gr_region
)

# --- Root universe ---
root_universe = openmc.Universe(name="Root universe")
root_universe.add_cells([core_cell, be_cell, gr_cell])

# --- Geometry object ---
geometry = openmc.Geometry(root_universe)
geometry.root_universe = root_universe


# =============================================================================
# Settings
# =============================================================================

print("Configuring settings...")

settings = openmc.Settings()
settings.batches = args.batches                # Total batches
settings.inactive = args.inactive              # Inactive batches (for source convergence)
settings.particles = args.particles            # Particles per batch
# Temperature treatment: nearest available temperature in cross section data
settings.temperature = {"method": "nearest", "tolerance": 100.0}

# --- Initial source ---
# Uniform box source spanning the core region
# This provides a reasonable initial fission source distribution
bounds = [
    -core_half_x, -core_half_y, -HALF_HEIGHT,   # Lower-left corner
    core_half_x, core_half_y, HALF_HEIGHT        # Upper-right corner
]
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(bounds[:3], bounds[3:]),
    constraints={"fissionable": True}
)

# --- Output settings ---
settings.output = {"tallies": True}


# =============================================================================
# Tallies
# =============================================================================

print("Setting up tallies...")

tallies = openmc.Tallies()

# --- Mesh tally for flux/fission rate distribution ---
# Create a regular mesh over the core region for visualization
if args.small_tallies:
    # Coarser mesh for faster runtime
    mesh_nx, mesh_ny = 50, 60
else:
    # Finer mesh for better resolution of plate structure
    mesh_nx, mesh_ny = 200, 240

mesh = openmc.RegularMesh(name="Core mesh")
mesh.dimension = [mesh_nx, mesh_ny, 1]
mesh.lower_left = [-core_half_x, -core_half_y, -HALF_HEIGHT]
mesh.upper_right = [core_half_x, core_half_y, HALF_HEIGHT]

# Fission rate tally on the mesh (for power distribution)
mesh_filter = openmc.MeshFilter(mesh)
fission_tally = openmc.Tally(name="Fission rate distribution")
fission_tally.filters = [mesh_filter]
fission_tally.scores = ["fission"]
tallies.append(fission_tally)

# Flux tally on the mesh
flux_tally = openmc.Tally(name="Flux distribution")
flux_tally.filters = [mesh_filter]
flux_tally.scores = ["flux"]
tallies.append(flux_tally)


# =============================================================================
# Export model
# =============================================================================

print("Exporting model to XML...")

model = openmc.Model(geometry=geometry, materials=materials,
                     settings=settings, tallies=tallies)
model.export_to_model_xml()

print(f"\nModel exported successfully!")
print(f"  Particles per batch: {args.particles}")
print(f"  Active batches:      {args.batches - args.inactive}")
print(f"  Inactive batches:    {args.inactive}")
print(f"  Core: {num_sfe} SFE ({NUM_PLATES_SFE} plates each) + "
      f"{num_cfe} CFE ({NUM_PLATES_CFE} plates each)")
print(f"  Total fuel plates:   {num_sfe * NUM_PLATES_SFE + num_cfe * NUM_PLATES_CFE}")
print(f"  Reflectors: Be ({BE_REFLECTOR_THICKNESS} cm) + "
      f"Graphite ({GRAPHITE_REFLECTOR_THICKNESS} cm)")
print(f"  Model type: 2D (XY) with reflective Z boundaries")
