"""
BEAVRS Full-Core PWR Benchmark Model
=====================================

MIT Benchmark for Evaluation And Validation of Reactor Simulations (BEAVRS)
Revision 2.0.1

Full-core model of a 4-loop Westinghouse PWR with 193 fuel assemblies.
Cycle 1 initial loading at Hot Zero Power (HZP) conditions with All Rods Out
(ARO) and 975 ppm soluble boron.

The core uses Westinghouse 17x17 fuel assembly design with three enrichment
zones in a low-leakage loading pattern:
  - 1.6 wt% U-235 (56 assemblies, center)
  - 2.4 wt% U-235 (64 assemblies, middle ring)
  - 3.1 wt% U-235 (73 assemblies, outer ring)

The geometry includes:
  - Exact pin dimensions from BEAVRS specification
  - 19x19 core lattice with stepped baffle (12 baffle universe types)
  - Core barrel, neutron shield region, and RPV
  - Burnable absorber (BPRA) rods per the exact loading map
  - All control rods fully withdrawn (ARO condition)

2D model (reflective Z boundaries) for initial validation.
3D axial detail to be added in a subsequent increment.

Reference k-eff at HZP, ARO, 975 ppm boron: ~0.99982 (measured critical).

References:
    N. Horelik et al., "Benchmark for Evaluation and Validation of Reactor
    Simulations (BEAVRS)," MIT Computational Reactor Physics Group, Rev. 2.0.1.
"""

import argparse

import numpy as np
import openmc
from openmc.data import atomic_mass, atomic_weight


# =============================================================================
# Cross-section data path
# =============================================================================
CROSS_SECTIONS_PATH = '/data/endfb-viii.0-hdf5/cross_sections.xml'

# =============================================================================
# BEAVRS Constants (Rev 2.0.1, from constants.py)
# =============================================================================

# Temperature (HZP isothermal, 560 F -> Kelvin)
HZP_TEMP = 566.483

# Pin cell dimensions [cm]
PELLET_OR = 0.39218
CLAD_IR = 0.40005
CLAD_OR = 0.45720
PIN_PITCH = 1.25984

# Guide tube dimensions [cm]
GT_IR = 0.56134
GT_OR = 0.60198
GT_DASH_IR = 0.50419  # Dashpot region (lower)
GT_DASH_OR = 0.54610

# Instrument tube dimensions [cm]
IT_IR = 0.43688
IT_OR = 0.48387

# Burnable absorber rod dimensions [cm]
BA_IR1 = 0.21400   # Inner air tube IR
BA_IR2 = 0.23051   # Inner SS tube OR
BA_IR3 = 0.24130   # Inner gap / glass IR
BA_IR4 = 0.42672   # Glass OR
BA_IR5 = 0.43688   # Outer gap / outer clad IR
BA_IR6 = 0.48387   # Outer clad OR

# Grid spacer dimensions [cm]
GRID_SIDE_TB = 1.22030   # Top/bottom grids
GRID_SIDE_I = 1.22098    # Intermediate grids

# Lattice parameters [cm]
LATTICE_PITCH = 21.50364
GRIDSTRAP_SIDE = 21.49595

# Core structural dimensions [cm]
CORE_BARREL_IR = 187.9600
CORE_BARREL_OR = 193.6750
BAFFLE_WIDTH = 2.2225
BAFFLE_WATER_GAP = 0.1627
NEUTRON_SHIELD_IR = 194.840
NEUTRON_SHIELD_OR = 201.630
LINER_IR = 219.150
RPV_IR = 219.710
RPV_OR = 241.300

# Axial dimensions [cm] (from ML033530020 and BEAVRS rod_axials.ods)
Z_LOWEST = 0.0
Z_SUPPORT_PLATE_BOT = 20.0
Z_SUPPORT_PLATE_TOP = 25.0
Z_LOWER_NOZZLE_BOT = 25.0
Z_LOWER_NOZZLE_TOP = 35.0

# Fuel rod axials
Z_FUEL_ROD_BOT = 35.0
Z_LOWER_FITTING_TOP = 36.748
Z_ACTIVE_FUEL_BOT = 36.748
Z_ACTIVE_FUEL_TOP = 402.508
Z_PLENUM_BOT = 402.508
Z_PLENUM_TOP = 417.164
Z_UPPER_FITTING_TOP = 419.704
Z_FUEL_ROD_TOP = 419.704

# Guide tube axials
Z_GT_DASHPOT_BOT = 35.0
Z_GT_DASHPOT_TOP = 39.958
Z_GT_UPPER_BOT = 39.958
Z_GT_UPPER_TOP = 423.049

# BPRA axials
Z_BPRA_ROD_BOT = 38.66
Z_BPRA_LOWER_FITTING_TOP = 40.558
Z_BPRA_ACTIVE_BOT = 40.558
Z_BPRA_ACTIVE_TOP = 401.238
Z_BPRA_PLENUM_TOP = 421.532
Z_BPRA_ROD_TOP = 431.876

# Instrument tube axials
Z_INSTR_BOT = 35.0
Z_INSTR_TOP = 423.049

# Upper core structures
Z_UPPER_NOZZLE_BOT = 423.049
Z_UPPER_NOZZLE_TOP = 431.876
Z_HIGHEST = 460.0

# Plenum spring OR
PLENUM_SPRING_OR = 0.06459

# Grid spacer positions [cm] (8 grids)
GRID_POSITIONS = [
    (37.1621, 40.5200),   # Grid 1 (bottom, Inconel)
    (98.0250, 103.740),   # Grid 2
    (150.222, 155.937),   # Grid 3
    (202.419, 208.134),   # Grid 4
    (254.616, 260.331),   # Grid 5
    (306.813, 312.528),   # Grid 6
    (359.010, 364.725),   # Grid 7
    (411.806, 415.164),   # Grid 8 (top, Inconel)
]

# Water density at HZP [g/cc]
H2O_DENSITY = 0.73986
BORON_PPM = 975

# Fuel enrichments and densities (BEAVRS spec, exact values)
FUEL_SPECS = {
    '1.6': {'enrichment': 0.0161006, 'density': 10.31341},
    '2.4': {'enrichment': 0.0239993, 'density': 10.29748},
    '3.1': {'enrichment': 0.0310221, 'density': 10.30166},
}

# =============================================================================
# Guide tube and instrument tube positions in 17x17 lattice (0-indexed)
# =============================================================================
GUIDE_TUBE_POSITIONS = [
    (2, 5),   (2, 8),   (2, 11),
    (3, 3),                       (3, 13),
    (5, 2),   (5, 5),   (5, 8),   (5, 11),  (5, 14),
    (8, 2),   (8, 5),             (8, 11),  (8, 14),
    (11, 2),  (11, 5),  (11, 8),  (11, 11), (11, 14),
    (13, 3),                      (13, 13),
    (14, 5),  (14, 8),  (14, 11),
]

INSTRUMENT_TUBE_POSITION = (8, 8)

# =============================================================================
# BPRA pin position maps (which of 24 guide tube positions get BA rods)
# Standard Westinghouse configurations from BEAVRS spec
# =============================================================================
# Positions are (row, col) in 17x17 lattice, subset of GUIDE_TUBE_POSITIONS

# 24 guide tube positions for reference (ordered)
_GT_SET = set(GUIDE_TUBE_POSITIONS)

# 20 BA rods: all 24 guide tubes except 4 face positions
BPRA_20 = [p for p in GUIDE_TUBE_POSITIONS if p not in
           [(2, 8), (8, 2), (8, 14), (14, 8)]]

# 16 BA rods: remove 8 outermost positions
BPRA_16 = [p for p in GUIDE_TUBE_POSITIONS if p not in
           [(2, 5), (2, 8), (2, 11), (5, 2), (5, 14),
            (8, 2), (8, 14), (14, 5), (14, 8), (14, 11),
            (11, 2), (11, 14), (3, 3), (3, 13), (13, 3), (13, 13)]]
# That's: (5,5),(5,8),(5,11),(8,5),(8,11),(11,5),(11,8),(11,11) = 8 positions
# Plus the 8 in the next ring
# Actually, standard 16-BA: inner 16 of the 24 positions
BPRA_16 = [
    (2, 5),            (2, 11),
    (3, 3),                       (3, 13),
    (5, 2),   (5, 5),   (5, 11),  (5, 14),
    (8, 2),                        (8, 14),
    (11, 2),  (11, 5),  (11, 11), (11, 14),
    (13, 3),                      (13, 13),
]

# 15 BA rods with orientation variants (one face position empty)
BPRA_15NW = [p for p in BPRA_16 if p != (2, 5)]     # Remove NW-ish position
BPRA_15NE = [p for p in BPRA_16 if p != (2, 11)]    # Remove NE-ish position
BPRA_15SW = [p for p in BPRA_16 if p != (13, 3)]    # Remove SW-ish position
BPRA_15SE = [p for p in BPRA_16 if p != (13, 13)]   # Remove SE-ish position

# 12 BA rods: alternating pattern
BPRA_12 = [
              (2, 5),             (2, 11),
    (5, 2),            (5, 8),             (5, 14),
              (8, 5),             (8, 11),
    (11, 2),           (11, 8),            (11, 14),
              (14, 5),            (14, 11),
]

# 6 BA rods with orientation variants (one face of the assembly)
BPRA_6N = [(2, 5), (2, 8), (2, 11), (3, 3), (3, 13), (5, 8)]
BPRA_6S = [(14, 5), (14, 8), (14, 11), (13, 3), (13, 13), (11, 8)]
BPRA_6E = [(5, 14), (8, 14), (11, 14), (3, 13), (13, 13), (8, 11)]
BPRA_6W = [(5, 2), (8, 2), (11, 2), (3, 3), (13, 3), (8, 5)]

BPRA_CONFIGS = {
    '20': BPRA_20, '16': BPRA_16,
    '15NW': BPRA_15NW, '15NE': BPRA_15NE,
    '15SW': BPRA_15SW, '15SE': BPRA_15SE,
    '12': BPRA_12,
    '6N': BPRA_6N, '6S': BPRA_6S, '6E': BPRA_6E, '6W': BPRA_6W,
}

# =============================================================================
# BPRA loading map (which assemblies get BAs, and which configuration)
# Keys are (row, col) in 19x19 core lattice (fuel positions only)
# Values are BPRA configuration strings
# =============================================================================
# Transcribed from BEAVRS core.py ba_positions, converted to 19x19 grid coords
# Row in 19x19 = BEAVRS_row_number + 1, Col mapping: R=2,P=3,N=4,...,A=16
BPRA_MAP = {
    # Row 2 (number 1)
    (2, 7): '6S', (2, 9): '6S', (2, 11): '6S',
    # Row 3 (number 2)
    (3, 6): '16', (3, 8): '20', (3, 10): '20', (3, 12): '16',
    # Row 4 (number 3)
    (4, 4): '15SE', (4, 5): '16', (4, 7): '16', (4, 9): '16', (4, 11): '16',
    (4, 13): '16', (4, 14): '15SW',
    # Row 5 (number 4)
    (5, 4): '16', (5, 6): '16', (5, 8): '12', (5, 10): '12', (5, 12): '16',
    (5, 14): '16',
    # Row 6 (number 5)
    (6, 3): '16', (6, 5): '16', (6, 7): '12', (6, 9): '12', (6, 11): '12',
    (6, 13): '16', (6, 15): '16',
    # Row 7 (number 6)
    (7, 2): '6E', (7, 4): '16', (7, 6): '12', (7, 8): '12', (7, 10): '12',
    (7, 12): '12', (7, 14): '16', (7, 16): '6W',
    # Row 8 (number 7)
    (8, 3): '20', (8, 5): '12', (8, 7): '12', (8, 9): '16', (8, 11): '12',
    (8, 13): '12', (8, 15): '20',
    # Row 9 (number 8, center)
    (9, 2): '6E', (9, 4): '16', (9, 6): '12', (9, 8): '16', (9, 10): '16',
    (9, 12): '12', (9, 14): '16', (9, 16): '6W',
    # Row 10 (number 9)
    (10, 3): '20', (10, 5): '12', (10, 7): '12', (10, 9): '16',
    (10, 11): '12', (10, 13): '12', (10, 15): '20',
    # Row 11 (number 10)
    (11, 2): '6E', (11, 4): '16', (11, 6): '12', (11, 8): '12',
    (11, 10): '12', (11, 12): '12', (11, 14): '16', (11, 16): '6W',
    # Row 12 (number 11)
    (12, 3): '16', (12, 5): '16', (12, 7): '12', (12, 9): '12',
    (12, 11): '12', (12, 13): '16', (12, 15): '16',
    # Row 13 (number 12)
    (13, 4): '16', (13, 6): '16', (13, 8): '12', (13, 10): '12',
    (13, 12): '16', (13, 14): '16',
    # Row 14 (number 13)
    (14, 4): '15NE', (14, 5): '16', (14, 7): '16', (14, 9): '16',
    (14, 11): '16', (14, 13): '16', (14, 14): '15NW',
    # Row 15 (number 14)
    (15, 6): '16', (15, 8): '20', (15, 10): '20', (15, 12): '16',
    # Row 16 (number 15)
    (16, 7): '6N', (16, 9): '6N', (16, 11): '6N',
}

# =============================================================================
# 19x19 Core Layout Map
# =============================================================================
# Cell type codes:
#   'W'   = water (outside core)
#   'F1'  = fuel 1.6%
#   'F2'  = fuel 2.4%
#   'F3'  = fuel 3.1%
#   'ES'  = edge baffle, steel faces south (-y)  [template 'bafn_', at north edge]
#   'EN'  = edge baffle, steel faces north (+y)  [template 'bafs_', at south edge]
#   'EE'  = edge baffle, steel faces east (+x)   [template 'bafw_', at west edge]
#   'EW'  = edge baffle, steel faces west (-x)   [template 'bafe_', at east edge]
#   'CSE' = corner baffle, steel faces SE (-y,+x) [template 'bfcnw']
#   'CSW' = corner baffle, steel faces SW (-y,-x) [template 'bfcne']
#   'CNE' = corner baffle, steel faces NE (+y,+x) [template 'bfcsw']
#   'CNW' = corner baffle, steel faces NW (+y,-x) [template 'bfcse']
#   'TSE' = tip baffle, steel faces SE            [template 'bafnw']
#   'TSW' = tip baffle, steel faces SW            [template 'bafne']
#   'TNE' = tip baffle, steel faces NE            [template 'bafsw']
#   'TNW' = tip baffle, steel faces NW            [template 'bafse']
#
# Row 0 = north (highest y), Row 18 = south (lowest y)
# Col 0 = west (lowest x), Col 18 = east (highest x)

CORE_MAP = [
  # Col: 0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
  ['W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W' ],  # 0
  ['W',  'W',  'W',  'W',  'W',  'TSE','ES', 'ES', 'ES', 'ES', 'ES', 'ES', 'ES', 'TSW','W',  'W',  'W',  'W',  'W' ],  # 1
  ['W',  'W',  'W',  'TSE','ES', 'CSE','F3', 'F3', 'F3', 'F3', 'F3', 'F3', 'F3', 'CSW','ES', 'TSW','W',  'W',  'W' ],  # 2
  ['W',  'W',  'TSE','CSE','F3', 'F3', 'F3', 'F1', 'F3', 'F1', 'F3', 'F1', 'F3', 'F3', 'F3', 'CSW','TSW','W',  'W' ],  # 3
  ['W',  'W',  'EE', 'F3', 'F3', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F3', 'F3', 'EW', 'W',  'W' ],  # 4
  ['W',  'TSE','CSE','F3', 'F2', 'F2', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F2', 'F2', 'F3', 'CSW','TSW','W' ],  # 5
  ['W',  'EE', 'F3', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'F3', 'EW', 'W' ],  # 6
  ['W',  'EE', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'EW', 'W' ],  # 7
  ['W',  'EE', 'F3', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'F3', 'EW', 'W' ],  # 8
  ['W',  'EE', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'EW', 'W' ],  # 9 center
  ['W',  'EE', 'F3', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'F3', 'EW', 'W' ],  # 10
  ['W',  'EE', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'EW', 'W' ],  # 11
  ['W',  'EE', 'F3', 'F3', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F3', 'F3', 'EW', 'W' ],  # 12
  ['W',  'TNE','CNE','F3', 'F2', 'F2', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F2', 'F2', 'F3', 'CNW','TNW','W' ],  # 13
  ['W',  'W',  'EE', 'F3', 'F3', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F1', 'F2', 'F3', 'F3', 'EW', 'W',  'W' ],  # 14
  ['W',  'W',  'TNE','CNE','F3', 'F3', 'F3', 'F1', 'F3', 'F1', 'F3', 'F1', 'F3', 'F3', 'F3', 'CNW','TNW','W',  'W' ],  # 15
  ['W',  'W',  'W',  'TNE','EN', 'CNE','F3', 'F3', 'F3', 'F3', 'F3', 'F3', 'F3', 'CNW','EN', 'TNW','W',  'W',  'W' ],  # 16
  ['W',  'W',  'W',  'W',  'W',  'TNE','EN', 'EN', 'EN', 'EN', 'EN', 'EN', 'EN', 'TNW','W',  'W',  'W',  'W',  'W' ],  # 17
  ['W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W',  'W' ],  # 18
]


# =============================================================================
# Material creation
# =============================================================================

def create_materials():
    """Create all materials per BEAVRS Rev 2.0.1 specification."""

    mats = {}

    # --- UO2 fuels (3 enrichments with exact isotopics including U-234) ---
    MB10 = atomic_mass('B10')
    MB11 = atomic_mass('B11')
    MU234 = atomic_mass('U234')
    MU235 = atomic_mass('U235')
    MU238 = atomic_mass('U238')

    for label, spec in FUEL_SPECS.items():
        enr_25 = spec['enrichment']
        enr_24 = 0.008 * enr_25      # U-234 fraction (0.8% of U-235)
        enr_28 = 1.0 - (enr_24 + enr_25)
        MU = 1.0 / (enr_24/MU234 + enr_25/MU235 + enr_28/MU238)
        MUO2 = MU + 2.0 * atomic_weight('O')
        wU = MU / MUO2
        a_U = wU * MUO2 / MU
        a_O = (1.0 - wU) * MUO2 / atomic_weight('O')

        fuel = openmc.Material(name=f'Fuel {label}%')
        fuel.set_density('g/cm3', spec['density'])
        fuel.add_element('O', a_O, 'ao')
        fuel.add_element('U', a_U, 'ao', enrichment=enr_25 * 100)
        fuel.temperature = HZP_TEMP
        mats[f'fuel_{label}'] = fuel

    # --- Zircaloy-4 (BEAVRS spec) ---
    zirc4 = openmc.Material(name='Zircaloy-4')
    zirc4.set_density('g/cm3', 6.55)
    zirc4.add_element('O',  0.00125, 'wo')
    zirc4.add_element('Cr', 0.0010, 'wo')
    zirc4.add_element('Fe', 0.0021, 'wo')
    zirc4.add_element('Zr', 0.98115, 'wo')
    zirc4.add_element('Sn', 0.0145, 'wo')
    zirc4.temperature = HZP_TEMP
    mats['zirc4'] = zirc4

    # --- Helium fill gas ---
    helium = openmc.Material(name='Helium')
    helium.set_density('g/cm3', 0.0015981)
    helium.add_element('He', 1.0)
    helium.temperature = HZP_TEMP
    mats['helium'] = helium

    # --- Borated water (975 ppm, NIST density at 2250 psia, 560 F) ---
    water = openmc.model.borated_water(
        BORON_PPM, density=H2O_DENSITY, name='Borated Water')
    water.temperature = HZP_TEMP
    mats['water'] = water

    # --- SS304 (BEAVRS spec) ---
    ss304 = openmc.Material(name='SS304')
    ss304.set_density('g/cm3', 8.03)
    ss304.add_element('Si', 0.0060, 'wo')
    ss304.add_element('Cr', 0.1900, 'wo')
    ss304.add_element('Mn', 0.0200, 'wo')
    ss304.add_element('Fe', 0.6840, 'wo')
    ss304.add_element('Ni', 0.1000, 'wo')
    ss304.temperature = HZP_TEMP
    mats['ss304'] = ss304

    # --- Borosilicate glass (for BPRA rods, CASMO weight fractions) ---
    wO = 0.5481
    wAl = 0.0344
    wSi = 0.3787
    wB10 = 0.0071
    wB11 = 0.0317
    M_bsg = 1.0 / (wO/atomic_weight('O') + wAl/atomic_weight('Al') +
                    wSi/atomic_weight('Si') + wB10/MB10 + wB11/MB11)

    bsg = openmc.Material(name='Borosilicate Glass')
    bsg.set_density('g/cm3', 2.26)
    bsg.add_element('O', wO * M_bsg / atomic_weight('O'), 'ao')
    bsg.add_element('Si', wSi * M_bsg / atomic_weight('Si'), 'ao')
    bsg.add_element('Al', wAl * M_bsg / atomic_weight('Al'), 'ao')
    bsg.add_nuclide('B10', wB10 * M_bsg / MB10, 'ao')
    bsg.add_nuclide('B11', wB11 * M_bsg / MB11, 'ao')
    bsg.temperature = HZP_TEMP
    mats['bsg'] = bsg

    # --- Air (for BPRA inner tubes) ---
    air = openmc.Material(name='Air')
    air.set_density('g/cm3', 0.000616)
    air.add_element('O', 0.2095, 'ao')
    air.add_element('N', 0.7809, 'ao')
    air.add_element('Ar', 0.00933, 'ao')
    air.add_nuclide('C12', 0.00027 * 0.9893, 'ao')
    air.add_nuclide('C13', 0.00027 * 0.0107, 'ao')
    air.temperature = HZP_TEMP
    mats['air'] = air

    # --- Inconel 718 (for grid spacers, future use) ---
    inconel = openmc.Material(name='Inconel 718')
    inconel.set_density('g/cm3', 8.2)
    inconel.add_element('Si', 0.0035, 'wo')
    inconel.add_element('Cr', 0.1896, 'wo')
    inconel.add_element('Mn', 0.0087, 'wo')
    inconel.add_element('Fe', 0.2863, 'wo')
    inconel.add_element('Ni', 0.5119, 'wo')
    inconel.temperature = HZP_TEMP
    mats['inconel'] = inconel

    # --- Carbon steel (for RPV) ---
    cs = openmc.Material(name='Carbon Steel')
    cs.set_density('g/cm3', 7.8)
    cs.add_nuclide('C12', 0.00270 * 0.9893, 'wo')
    cs.add_nuclide('C13', 0.00270 * 0.0107, 'wo')
    cs.add_element('Mn', 0.00750, 'wo')
    cs.add_element('P',  0.00025, 'wo')
    cs.add_element('S',  0.00025, 'wo')
    cs.add_element('Si', 0.00400, 'wo')
    cs.add_element('Ni', 0.00750, 'wo')
    cs.add_element('Cr', 0.00350, 'wo')
    cs.add_element('Mo', 0.00625, 'wo')
    cs.add_element('V',  0.00050, 'wo')
    cs.add_element('Nb', 0.00010, 'wo')
    cs.add_element('Cu', 0.00200, 'wo')
    cs.add_element('Ca', 0.00015, 'wo')
    cs.add_element('B',  0.00003, 'wo')
    cs.add_element('Ti', 0.00015, 'wo')
    cs.add_element('Al', 0.00025, 'wo')
    cs.add_element('Fe', 0.96487, 'wo')
    cs.temperature = HZP_TEMP
    mats['carbon_steel'] = cs

    return mats


# =============================================================================
# Pin universe constructors
# =============================================================================

def create_fuel_pin(fuel_mat, zirc4, helium, water, zplanes):
    """Create a 3D fuel pin universe with axial structure.

    Axial regions (bottom to top):
      - Below rod:     water
      - Lower fitting: solid Zirc4 plug (35.0 - 36.748 cm)
      - Active fuel:   UO2 pellet + He gap + Zirc4 clad (36.748 - 402.508 cm)
      - Gas plenum:    He + spring + clad (402.508 - 417.164 cm)
      - Upper fitting: solid Zirc4 plug (417.164 - 419.704 cm)
      - Above rod:     water
    """
    # Radial surfaces
    pellet = openmc.ZCylinder(r=PELLET_OR)
    clad_ir = openmc.ZCylinder(r=CLAD_IR)
    clad_or = openmc.ZCylinder(r=CLAD_OR)
    spring = openmc.ZCylinder(r=PLENUM_SPRING_OR)

    # Axial regions
    z_rod_bot = zplanes['fuel_rod_bot']
    z_fit_top = zplanes['lower_fitting_top']
    z_fuel_top = zplanes['active_fuel_top']
    z_plen_top = zplanes['plenum_top']
    z_ufit_top = zplanes['upper_fitting_top']

    u = openmc.Universe(name=f'Fuel Pin ({fuel_mat.name})')
    cells = []

    # Below fuel rod: water
    cells.append(openmc.Cell(fill=water, region=-z_rod_bot))

    # Lower fitting: solid Zirc4 inside cladding OD
    z_lf = +z_rod_bot & -z_fit_top
    cells.append(openmc.Cell(fill=zirc4, region=-clad_or & z_lf))
    cells.append(openmc.Cell(fill=water, region=+clad_or & z_lf))

    # Active fuel region
    z_af = +z_fit_top & -z_fuel_top
    cells.append(openmc.Cell(fill=fuel_mat, region=-pellet & z_af))
    cells.append(openmc.Cell(fill=helium, region=+pellet & -clad_ir & z_af))
    cells.append(openmc.Cell(fill=zirc4, region=+clad_ir & -clad_or & z_af))
    cells.append(openmc.Cell(fill=water, region=+clad_or & z_af))

    # Gas plenum: spring + helium + clad
    z_pl = +z_fuel_top & -z_plen_top
    cells.append(openmc.Cell(fill=zirc4, region=-spring & z_pl))
    cells.append(openmc.Cell(fill=helium, region=+spring & -clad_ir & z_pl))
    cells.append(openmc.Cell(fill=zirc4, region=+clad_ir & -clad_or & z_pl))
    cells.append(openmc.Cell(fill=water, region=+clad_or & z_pl))

    # Upper fitting: solid Zirc4
    z_uf = +z_plen_top & -z_ufit_top
    cells.append(openmc.Cell(fill=zirc4, region=-clad_or & z_uf))
    cells.append(openmc.Cell(fill=water, region=+clad_or & z_uf))

    # Above fuel rod: water
    cells.append(openmc.Cell(fill=water, region=+z_ufit_top))

    u.add_cells(cells)
    return u


def create_guide_tube(zirc4, water, zplanes):
    """Create a 3D guide tube universe with dashpot region.

    Axial regions:
      - Below GT:    water (z < 35.0)
      - Dashpot:     smaller IR/OR Zirc4 tube (35.0 - 39.958 cm)
      - Upper:       standard IR/OR Zirc4 tube (39.958 - 423.049 cm)
      - Above GT:    water (z > 423.049)
    """
    # Dashpot radii (smaller, thicker wall)
    gt_dash_ir = openmc.ZCylinder(r=GT_DASH_IR)
    gt_dash_or = openmc.ZCylinder(r=GT_DASH_OR)
    # Upper region radii
    gt_ir = openmc.ZCylinder(r=GT_IR)
    gt_or = openmc.ZCylinder(r=GT_OR)

    z_bot = zplanes['gt_bot']
    z_dash_top = zplanes['gt_dash_top']
    z_top = zplanes['gt_top']

    u = openmc.Universe(name='Guide Tube')
    cells = []

    # Below GT: water
    cells.append(openmc.Cell(fill=water, region=-z_bot))

    # Dashpot region (smaller tube)
    z_d = +z_bot & -z_dash_top
    cells.append(openmc.Cell(fill=water, region=-gt_dash_ir & z_d))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_dash_ir & -gt_dash_or & z_d))
    cells.append(openmc.Cell(fill=water, region=+gt_dash_or & z_d))

    # Upper region (standard tube)
    z_u = +z_dash_top & -z_top
    cells.append(openmc.Cell(fill=water, region=-gt_ir & z_u))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_ir & -gt_or & z_u))
    cells.append(openmc.Cell(fill=water, region=+gt_or & z_u))

    # Above GT: water
    cells.append(openmc.Cell(fill=water, region=+z_top))

    u.add_cells(cells)
    return u


def create_instrument_tube(zirc4, water, zplanes):
    """Create a 3D instrument tube universe (35.0 - 423.049 cm)."""
    it_ir = openmc.ZCylinder(r=IT_IR)
    it_or = openmc.ZCylinder(r=IT_OR)

    z_bot = zplanes['instr_bot']
    z_top = zplanes['instr_top']

    u = openmc.Universe(name='Instrument Tube')
    cells = []

    # Below tube: water
    cells.append(openmc.Cell(fill=water, region=-z_bot))

    # Tube region
    z_t = +z_bot & -z_top
    cells.append(openmc.Cell(fill=water, region=-it_ir & z_t))
    cells.append(openmc.Cell(fill=zirc4, region=+it_ir & -it_or & z_t))
    cells.append(openmc.Cell(fill=water, region=+it_or & z_t))

    # Above tube: water
    cells.append(openmc.Cell(fill=water, region=+z_top))

    u.add_cells(cells)
    return u


def create_bpra_pin(ss304, air, bsg, water, zirc4, zplanes):
    """Create a 3D BPRA pin universe.

    Axial regions:
      - Below rod:      water + guide tube (z < 38.66)
      - Lower fitting:  SS304 plug (38.66 - 40.558)
      - Active absorber: 6-zone radial structure (40.558 - 401.238)
      - Plenum:         He/air + SS clad (401.238 - 421.532)
      - Above rod:      water (> 421.532)

    The BPRA sits inside the guide tube, so below/above the BPRA rod,
    the guide tube structure is present.
    """
    # BPRA radial surfaces
    s1 = openmc.ZCylinder(r=BA_IR1)
    s2 = openmc.ZCylinder(r=BA_IR2)
    s3 = openmc.ZCylinder(r=BA_IR3)
    s4 = openmc.ZCylinder(r=BA_IR4)
    s5 = openmc.ZCylinder(r=BA_IR5)
    s6 = openmc.ZCylinder(r=BA_IR6)

    # Guide tube surfaces (BPRA replaces the guide tube interior)
    gt_ir = openmc.ZCylinder(r=GT_IR)
    gt_or = openmc.ZCylinder(r=GT_OR)
    gt_dash_ir = openmc.ZCylinder(r=GT_DASH_IR)
    gt_dash_or = openmc.ZCylinder(r=GT_DASH_OR)

    z_gt_bot = zplanes['gt_bot']
    z_gt_dash_top = zplanes['gt_dash_top']
    z_bpra_bot = zplanes['bpra_bot']
    z_bpra_fit_top = zplanes['bpra_fitting_top']
    z_bpra_active_top = zplanes['bpra_active_top']
    z_bpra_plenum_top = zplanes['bpra_plenum_top']
    z_bpra_top = zplanes['bpra_top']
    z_gt_top = zplanes['gt_top']

    u = openmc.Universe(name='BPRA Pin')
    cells = []

    # Below guide tube: water
    cells.append(openmc.Cell(fill=water, region=-z_gt_bot))

    # Guide tube dashpot region below BPRA rod (35.0 - 38.66)
    z_d = +z_gt_bot & -z_bpra_bot
    cells.append(openmc.Cell(fill=water, region=-gt_dash_ir & z_d))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_dash_ir & -gt_dash_or & z_d))
    cells.append(openmc.Cell(fill=water, region=+gt_dash_or & z_d))

    # BPRA lower fitting (38.66 - 40.558): SS plug inside GT
    z_lf = +z_bpra_bot & -z_bpra_fit_top
    cells.append(openmc.Cell(fill=ss304, region=-s6 & z_lf))
    cells.append(openmc.Cell(fill=water, region=+s6 & -gt_dash_ir & z_lf))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_dash_ir & -gt_dash_or & z_lf))
    cells.append(openmc.Cell(fill=water, region=+gt_dash_or & z_lf))

    # Active absorber region (40.558 - 401.238): 6-zone inside GT upper
    z_a = +z_bpra_fit_top & -z_bpra_active_top
    cells.append(openmc.Cell(fill=air, region=-s1 & z_a))
    cells.append(openmc.Cell(fill=ss304, region=+s1 & -s2 & z_a))
    cells.append(openmc.Cell(fill=air, region=+s2 & -s3 & z_a))
    cells.append(openmc.Cell(fill=bsg, region=+s3 & -s4 & z_a))
    cells.append(openmc.Cell(fill=air, region=+s4 & -s5 & z_a))
    cells.append(openmc.Cell(fill=ss304, region=+s5 & -s6 & z_a))
    cells.append(openmc.Cell(fill=water, region=+s6 & -gt_ir & z_a))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_ir & -gt_or & z_a))
    cells.append(openmc.Cell(fill=water, region=+gt_or & z_a))

    # BPRA plenum (401.238 - 421.532): air/He + SS clad inside GT
    z_p = +z_bpra_active_top & -z_bpra_plenum_top
    cells.append(openmc.Cell(fill=air, region=-s6 & z_p))
    cells.append(openmc.Cell(fill=water, region=+s6 & -gt_ir & z_p))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_ir & -gt_or & z_p))
    cells.append(openmc.Cell(fill=water, region=+gt_or & z_p))

    # Above BPRA rod to GT top (421.532 - 423.049): water in GT
    z_above = +z_bpra_plenum_top & -z_gt_top
    cells.append(openmc.Cell(fill=water, region=-gt_ir & z_above))
    cells.append(openmc.Cell(fill=zirc4, region=+gt_ir & -gt_or & z_above))
    cells.append(openmc.Cell(fill=water, region=+gt_or & z_above))

    # Above guide tube: water
    cells.append(openmc.Cell(fill=water, region=+z_gt_top))

    u.add_cells(cells)
    return u


# =============================================================================
# Assembly constructor
# =============================================================================

def create_assembly(enrichment_label, materials, zplanes, bpra_config=None):
    """
    Create a 17x17 fuel assembly universe with 3D axial detail.

    Parameters
    ----------
    enrichment_label : str
        '1.6', '2.4', or '3.1'
    materials : dict
        Materials dictionary
    zplanes : dict
        Dictionary of shared ZPlane surfaces
    bpra_config : str or None
        BPRA configuration key (e.g. '20', '16', '12', '6N'), or None
    """
    fuel_mat = materials[f'fuel_{enrichment_label}']
    zirc4 = materials['zirc4']
    helium = materials['helium']
    water = materials['water']

    fuel_pin = create_fuel_pin(fuel_mat, zirc4, helium, water, zplanes)
    guide_tube = create_guide_tube(zirc4, water, zplanes)
    instrument_tube = create_instrument_tube(zirc4, water, zplanes)

    # Create BPRA pin if needed
    bpra_pin = None
    bpra_positions = set()
    if bpra_config is not None:
        bpra_pin = create_bpra_pin(
            materials['ss304'], materials['air'], materials['bsg'],
            water, zirc4, zplanes)
        bpra_positions = set(BPRA_CONFIGS[bpra_config])

    # Build 17x17 pin array
    pin_array = np.full((17, 17), fuel_pin, dtype=openmc.Universe)

    for pos in GUIDE_TUBE_POSITIONS:
        if pos in bpra_positions:
            pin_array[pos[0], pos[1]] = bpra_pin
        else:
            pin_array[pos[0], pos[1]] = guide_tube

    pin_array[INSTRUMENT_TUBE_POSITION[0],
              INSTRUMENT_TUBE_POSITION[1]] = instrument_tube

    # Create lattice
    assembly = openmc.RectLattice(
        name=f'Assembly {enrichment_label}% BA={bpra_config}')
    assembly.pitch = (PIN_PITCH, PIN_PITCH)
    assembly.lower_left = (-17 * PIN_PITCH / 2.0, -17 * PIN_PITCH / 2.0)
    assembly.universes = pin_array
    assembly.outer = openmc.Universe()
    assembly.outer.add_cell(openmc.Cell(fill=water))

    # Wrap in universe
    u = openmc.Universe(
        name=f'Assem {enrichment_label}% BA={bpra_config}')
    u.add_cell(openmc.Cell(fill=assembly))
    return u


# =============================================================================
# Baffle universe constructors
# =============================================================================

def _create_edge_baffle(name, axis, sign, ss304, water):
    """
    Create an edge baffle universe with a steel plate near one cell edge.

    Parameters
    ----------
    axis : str
        'x' or 'y'
    sign : int
        +1 = steel near positive edge, -1 = steel near negative edge
    """
    hp = LATTICE_PITCH / 2
    gap_pos = sign * (hp - BAFFLE_WATER_GAP)
    steel_pos = sign * (hp - BAFFLE_WATER_GAP - BAFFLE_WIDTH)

    if axis == 'y':
        gap_surf = openmc.YPlane(y0=gap_pos)
        steel_surf = openmc.YPlane(y0=steel_pos)
    else:
        gap_surf = openmc.XPlane(x0=gap_pos)
        steel_surf = openmc.XPlane(x0=steel_pos)

    u = openmc.Universe(name=name)
    if sign > 0:
        # Steel near +axis edge
        u.add_cells([
            openmc.Cell(fill=water, region=+gap_surf),           # Water gap
            openmc.Cell(fill=ss304, region=-gap_surf & +steel_surf),  # Steel
            openmc.Cell(fill=water, region=-steel_surf),         # Outer water
        ])
    else:
        # Steel near -axis edge
        u.add_cells([
            openmc.Cell(fill=water, region=-gap_surf),           # Water gap
            openmc.Cell(fill=ss304, region=+gap_surf & -steel_surf),  # Steel
            openmc.Cell(fill=water, region=+steel_surf),         # Outer water
        ])
    return u


def _create_corner_baffle(name, y_sign, x_sign, ss304, water):
    """
    Create a corner baffle universe with L-shaped steel.

    The steel faces toward (y_sign, x_sign) direction.
    Water gap strips are along the outer edges (same direction as steel).
    Outer water fills the opposite corner.
    """
    hp = LATTICE_PITCH / 2

    # Y-direction surfaces (for horizontal plate)
    y_gap_pos = y_sign * (hp - BAFFLE_WATER_GAP)
    y_steel_pos = y_sign * (hp - BAFFLE_WATER_GAP - BAFFLE_WIDTH)

    # X-direction surfaces (for vertical plate)
    x_gap_pos = x_sign * (hp - BAFFLE_WATER_GAP)
    x_steel_pos = x_sign * (hp - BAFFLE_WATER_GAP - BAFFLE_WIDTH)

    y_gap = openmc.YPlane(y0=y_gap_pos)
    y_steel = openmc.YPlane(y0=y_steel_pos)
    x_gap = openmc.XPlane(x0=x_gap_pos)
    x_steel = openmc.XPlane(x0=x_steel_pos)

    u = openmc.Universe(name=name)

    # Define half-space selectors based on signs
    # "toward steel" = in the direction of sign
    # "away from steel" = opposite direction
    if y_sign > 0:
        y_beyond_gap = +y_gap      # y > gap (thin strip at +y edge)
        y_in_steel = -y_gap & +y_steel  # between gap and steel surfaces
        y_past_steel = -y_steel    # past the steel toward center
    else:
        y_beyond_gap = -y_gap
        y_in_steel = +y_gap & -y_steel
        y_past_steel = +y_steel

    if x_sign > 0:
        x_beyond_gap = +x_gap
        x_in_steel = -x_gap & +x_steel
        x_past_steel = -x_steel
    else:
        x_beyond_gap = -x_gap
        x_in_steel = +x_gap & -x_steel
        x_past_steel = +x_steel

    u.add_cells([
        # Water gap strip along y-edge (entire row)
        openmc.Cell(fill=water, region=y_beyond_gap),
        # Water gap strip along x-edge (below y-gap)
        openmc.Cell(fill=water, region=~y_beyond_gap & x_beyond_gap),
        # Horizontal steel plate (between y surfaces, past x gap)
        openmc.Cell(fill=ss304, region=y_in_steel & ~x_beyond_gap),
        # Vertical steel plate (past y steel, between x surfaces)
        openmc.Cell(fill=ss304, region=y_past_steel & x_in_steel),
        # Outer water (past both steel plates)
        openmc.Cell(fill=water, region=y_past_steel & x_past_steel),
    ])
    return u


def _create_tip_baffle(name, y_sign, x_sign, ss304, water):
    """
    Create a tip baffle universe with L-shaped steel at a corner.

    Similar to corner baffle but water gap is only at the corner intersection,
    and steel extends to the cell edges. Three water regions fill the rest.
    """
    hp = LATTICE_PITCH / 2

    y_gap_pos = y_sign * (hp - BAFFLE_WATER_GAP)
    y_steel_pos = y_sign * (hp - BAFFLE_WATER_GAP - BAFFLE_WIDTH)
    x_gap_pos = x_sign * (hp - BAFFLE_WATER_GAP)
    x_steel_pos = x_sign * (hp - BAFFLE_WATER_GAP - BAFFLE_WIDTH)

    y_gap = openmc.YPlane(y0=y_gap_pos)
    y_steel = openmc.YPlane(y0=y_steel_pos)
    x_gap = openmc.XPlane(x0=x_gap_pos)
    x_steel = openmc.XPlane(x0=x_steel_pos)

    u = openmc.Universe(name=name)

    if y_sign > 0:
        y_beyond_gap = +y_gap
        y_in_gap_or_steel = -y_gap & +y_steel  # between y_gap and y_steel
        y_past_steel = -y_steel
    else:
        y_beyond_gap = -y_gap
        y_in_gap_or_steel = +y_gap & -y_steel
        y_past_steel = +y_steel

    if x_sign > 0:
        x_beyond_gap = +x_gap
        x_between = -x_gap & +x_steel
        x_past_steel = -x_steel
    else:
        x_beyond_gap = -x_gap
        x_between = +x_gap & -x_steel
        x_past_steel = +x_steel

    u.add_cells([
        # Water gap: corner where both gap surfaces are exceeded
        openmc.Cell(fill=water, region=y_beyond_gap & x_beyond_gap),
        # Steel: L-shape
        # Horizontal arm: between y surfaces, between x surfaces
        openmc.Cell(fill=ss304, region=y_in_gap_or_steel & x_between),
        # Vertical arm: between y_gap and y_steel, beyond x gap (but not
        # in the water gap corner). Plus area below y_steel between x surfaces.
        # Actually following BEAVRS reference exactly:
        # C2 (steel N arm): y > y_steel AND x < x_steel AND x > x_gap
        # C3 (steel W arm): y < y_gap AND y > y_steel AND x < x_gap
        # Let me redefine for the general case.
    ])
    # Clear and redo with exact BEAVRS reference logic
    u._cells.clear()

    # Following the reference baffle.py NWT pattern exactly, generalized:
    # The tip has steel in an L-shape with water gap only at the corner
    # intersection. Using the NW tip as reference and generalizing signs.

    # For NW tip (y_sign=+1, x_sign=-1):
    #   C1: y > y_gap AND x < x_gap  (water gap at +y,-x corner)
    #   C2: y > y_steel AND x_gap < x < x_steel  (horizontal steel strip)
    #   C3: y_steel < y < y_gap AND x < x_gap  (vertical steel strip)
    #   C4: y > y_steel AND x > x_steel  (water beyond both plates, +x side)
    #   C5: y < y_steel AND x < x_steel  (water beyond both plates, -y side)
    #   C6: y < y_steel AND x > x_steel  (water in opposite corner)

    # Generalized: "beyond_gap" = past gap in steel direction
    #              "between" = between gap and steel surfaces on that axis
    #              "past_steel" = past steel away from edge

    # For y_sign>0: beyond_gap = +y_gap, between = -y_gap & +y_steel,
    #               past_steel = -y_steel
    # For y_sign<0: beyond_gap = -y_gap, between = +y_gap & -y_steel,
    #               past_steel = +y_steel
    # Similar for x

    # Need: "in_steel_zone" = between steel surface and edge (includes gap)
    if y_sign > 0:
        y_in_steel_zone = +y_steel    # y > y_steel (gap + steel + beyond)
        y_in_steel_only = -y_gap & +y_steel
        y_past = -y_steel
    else:
        y_in_steel_zone = -y_steel
        y_in_steel_only = +y_gap & -y_steel
        y_past = +y_steel

    if x_sign > 0:
        x_in_steel_zone = +x_steel
        x_in_steel_only = -x_gap & +x_steel
        x_past = -x_steel
    else:
        x_in_steel_zone = -x_steel
        x_in_steel_only = +x_gap & -x_steel
        x_past = +x_steel

    u.add_cells([
        # C1: Water gap (corner intersection of both gap regions)
        openmc.Cell(fill=water, region=y_beyond_gap & x_beyond_gap),
        # C2: Steel arm along y-edge (between x surfaces, in y steel zone)
        openmc.Cell(fill=ss304, region=y_in_steel_zone & x_in_steel_only),
        # C3: Steel arm along x-edge (between y surfaces, beyond x gap)
        openmc.Cell(fill=ss304, region=y_in_steel_only & x_beyond_gap),
        # C4: Water past y-steel on +x side (or whichever is the "other" side)
        openmc.Cell(fill=water, region=y_in_steel_zone & x_past),
        # C5: Water past both plates (opposite corner from steel)
        openmc.Cell(fill=water, region=y_past & x_in_steel_zone),
        # C6: Water in remaining quadrant
        openmc.Cell(fill=water, region=y_past & x_past),
    ])
    return u


def create_baffle_universes(ss304, water):
    """Create all 12 baffle universe types."""
    baffles = {}

    # Edge baffles (steel faces toward core center)
    baffles['EN'] = _create_edge_baffle('Baffle Edge N', 'y', +1, ss304, water)
    baffles['ES'] = _create_edge_baffle('Baffle Edge S', 'y', -1, ss304, water)
    baffles['EE'] = _create_edge_baffle('Baffle Edge E', 'x', +1, ss304, water)
    baffles['EW'] = _create_edge_baffle('Baffle Edge W', 'x', -1, ss304, water)

    # Corner baffles (L-shaped steel facing toward core center)
    baffles['CNW'] = _create_corner_baffle('Baffle Corner NW', +1, -1, ss304, water)
    baffles['CNE'] = _create_corner_baffle('Baffle Corner NE', +1, +1, ss304, water)
    baffles['CSW'] = _create_corner_baffle('Baffle Corner SW', -1, -1, ss304, water)
    baffles['CSE'] = _create_corner_baffle('Baffle Corner SE', -1, +1, ss304, water)

    # Tip baffles (L-shaped steel at step transitions)
    baffles['TNW'] = _create_tip_baffle('Baffle Tip NW', +1, -1, ss304, water)
    baffles['TNE'] = _create_tip_baffle('Baffle Tip NE', +1, +1, ss304, water)
    baffles['TSW'] = _create_tip_baffle('Baffle Tip SW', -1, -1, ss304, water)
    baffles['TSE'] = _create_tip_baffle('Baffle Tip SE', -1, +1, ss304, water)

    return baffles


# =============================================================================
# Model builder
# =============================================================================

def build_model(particles=20000, batches=200, inactive=50):
    """Build the complete BEAVRS full-core 3D model."""

    model = openmc.model.Model()

    # =========================================================================
    # Materials
    # =========================================================================
    materials = create_materials()
    model.materials = openmc.Materials(materials.values())
    model.materials.cross_sections = CROSS_SECTIONS_PATH

    # =========================================================================
    # Shared Z-plane surfaces for 3D axial structure
    # =========================================================================
    zplanes = {
        'fuel_rod_bot':      openmc.ZPlane(z0=Z_FUEL_ROD_BOT),
        'lower_fitting_top': openmc.ZPlane(z0=Z_LOWER_FITTING_TOP),
        'active_fuel_top':   openmc.ZPlane(z0=Z_ACTIVE_FUEL_TOP),
        'plenum_top':        openmc.ZPlane(z0=Z_PLENUM_TOP),
        'upper_fitting_top': openmc.ZPlane(z0=Z_UPPER_FITTING_TOP),
        'gt_bot':            openmc.ZPlane(z0=Z_GT_DASHPOT_BOT),
        'gt_dash_top':       openmc.ZPlane(z0=Z_GT_DASHPOT_TOP),
        'gt_top':            openmc.ZPlane(z0=Z_GT_UPPER_TOP),
        'instr_bot':         openmc.ZPlane(z0=Z_INSTR_BOT),
        'instr_top':         openmc.ZPlane(z0=Z_INSTR_TOP),
        'bpra_bot':          openmc.ZPlane(z0=Z_BPRA_ROD_BOT),
        'bpra_fitting_top':  openmc.ZPlane(z0=Z_BPRA_LOWER_FITTING_TOP),
        'bpra_active_top':   openmc.ZPlane(z0=Z_BPRA_ACTIVE_TOP),
        'bpra_plenum_top':   openmc.ZPlane(z0=Z_BPRA_PLENUM_TOP),
        'bpra_top':          openmc.ZPlane(z0=Z_BPRA_ROD_TOP),
    }

    # =========================================================================
    # Assembly universes (keyed by enrichment + BPRA config)
    # =========================================================================
    assembly_cache = {}

    def get_assembly(row, col, enrichment_label):
        """Get or create assembly universe for given position."""
        bpra = BPRA_MAP.get((row, col))
        key = (enrichment_label, bpra)
        if key not in assembly_cache:
            assembly_cache[key] = create_assembly(
                enrichment_label, materials, zplanes, bpra)
        return assembly_cache[key]

    # Water universe for empty positions
    water_univ = openmc.Universe(name='Water')
    water_univ.add_cell(openmc.Cell(fill=materials['water']))

    # Baffle universes
    baffles = create_baffle_universes(materials['ss304'], materials['water'])

    # =========================================================================
    # Build 19x19 core lattice
    # =========================================================================
    enrichment_map = {'F1': '1.6', 'F2': '2.4', 'F3': '3.1'}
    n_fuel = 0

    core_array = np.empty((19, 19), dtype=openmc.Universe)
    for i in range(19):
        for j in range(19):
            code = CORE_MAP[i][j]
            if code == 'W':
                core_array[i, j] = water_univ
            elif code.startswith('F'):
                enr = enrichment_map[code]
                core_array[i, j] = get_assembly(i, j, enr)
                n_fuel += 1
            else:
                core_array[i, j] = baffles[code]

    print(f"Core loading: {n_fuel} fuel assemblies (expected 193)")
    print(f"Unique assembly types: {len(assembly_cache)}")
    assert n_fuel == 193, f"Expected 193 assemblies, got {n_fuel}"

    core_lattice = openmc.RectLattice(name='Core Lattice')
    core_lattice.pitch = (LATTICE_PITCH, LATTICE_PITCH)
    core_lattice.lower_left = (-19 * LATTICE_PITCH / 2.0,
                                -19 * LATTICE_PITCH / 2.0)
    core_lattice.universes = core_array
    core_lattice.outer = water_univ

    # =========================================================================
    # Core structures with 3D Z-boundaries
    # =========================================================================
    # Radial surfaces
    barrel_ir = openmc.ZCylinder(r=CORE_BARREL_IR, name='Barrel IR')
    barrel_or = openmc.ZCylinder(r=CORE_BARREL_OR, name='Barrel OR')
    rpv_ir = openmc.ZCylinder(r=RPV_IR, name='RPV IR')
    rpv_or = openmc.ZCylinder(r=RPV_OR, name='RPV OR',
                               boundary_type='vacuum')

    # Axial boundaries
    z_bot = openmc.ZPlane(z0=Z_LOWEST, boundary_type='vacuum')
    z_top = openmc.ZPlane(z0=Z_HIGHEST, boundary_type='vacuum')
    z_axial = +z_bot & -z_top

    # Core lattice cell (inside barrel, full axial extent)
    # The pin universes handle their own axial structure internally
    core_cell = openmc.Cell(name='Core', fill=core_lattice,
                            region=-barrel_ir & z_axial)

    # Core barrel (full height)
    barrel_cell = openmc.Cell(name='Core Barrel', fill=materials['ss304'],
                              region=+barrel_ir & -barrel_or & z_axial)

    # Downcomer water (between barrel and RPV)
    downcomer_cell = openmc.Cell(name='Downcomer', fill=materials['water'],
                                 region=+barrel_or & -rpv_ir & z_axial)

    # RPV (carbon steel)
    rpv_cell = openmc.Cell(name='RPV', fill=materials['carbon_steel'],
                           region=+rpv_ir & -rpv_or & z_axial)

    root = openmc.Universe(name='Root')
    root.add_cells([core_cell, barrel_cell, downcomer_cell, rpv_cell])
    model.geometry = openmc.Geometry(root)

    # =========================================================================
    # Settings
    # =========================================================================
    model.settings.batches = batches
    model.settings.inactive = inactive
    model.settings.particles = particles
    model.settings.run_mode = 'eigenvalue'

    # Source over active fuel region (3D box)
    src_hw = 15 * LATTICE_PITCH / 2.0
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(
            (-src_hw, -src_hw, Z_ACTIVE_FUEL_BOT),
            (src_hw, src_hw, Z_ACTIVE_FUEL_TOP)
        ),
        constraints={'fissionable': True}
    )

    model.settings.temperature = {
        'method': 'interpolation',
        'default': HZP_TEMP,
    }

    # =========================================================================
    # Tallies
    # =========================================================================
    # 3D assembly power mesh (19x19 radial, 1 axial bin for now)
    mesh = openmc.RegularMesh()
    mesh.dimension = [19, 19, 1]
    mesh.lower_left = [-19 * LATTICE_PITCH / 2, -19 * LATTICE_PITCH / 2,
                       Z_ACTIVE_FUEL_BOT]
    mesh.upper_right = [19 * LATTICE_PITCH / 2, 19 * LATTICE_PITCH / 2,
                        Z_ACTIVE_FUEL_TOP]

    power_tally = openmc.Tally(name='Assembly Power')
    power_tally.filters = [openmc.MeshFilter(mesh)]
    power_tally.scores = ['fission']

    model.tallies = openmc.Tallies([power_tally])

    return model


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='BEAVRS Full-Core PWR Benchmark (Rev 2.0.1)')
    parser.add_argument('--particles', type=int, default=20000)
    parser.add_argument('--batches', type=int, default=200)
    parser.add_argument('--inactive', type=int, default=50)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--export', action='store_true')
    args = parser.parse_args()

    if not args.run and not args.export:
        parser.print_help()
        print("\nSpecify --export to generate XML files or --run to execute.")
        return

    print("=" * 70)
    print("BEAVRS Full-Core PWR Benchmark (Rev 2.0.1)")
    print("=" * 70)
    print(f"  Particles: {args.particles}  Batches: {args.batches}"
          f"  Inactive: {args.inactive}")

    model = build_model(args.particles, args.batches, args.inactive)

    if args.export:
        model.export_to_xml()
        print("Model exported to XML.")

    if args.run:
        print("Starting eigenvalue calculation...")
        sp = model.run()
        print(f"Done. Statepoint: {sp}")


if __name__ == '__main__':
    main()
