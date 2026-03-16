#!/usr/bin/env python3
"""
OECD/NEA VVER-1000 MOX Core Computational Benchmark - Full Core Model
======================================================================

NEA/NSC/DOC(2005)17: "VVER-1000 MOX Core Computational Benchmark"
V. Boyarinov et al., OECD Nuclear Energy Agency, 2005.

This is a 2-D benchmark (infinite axial dimension, vacuum on the side surface)
of a VVER-1000 core with heterogeneous 30% MOX-fuel loading. The core has 163
hexagonal fuel assemblies with 60-degree rotational symmetry (28 unique
assemblies in the symmetry sector).

State S1 (Working State)
------------------------
  - Fuel temperature:      1027 K
  - Moderator in FA:       575 K, 1300 ppm boron, 0.7241 g/cm3
  - Water in reflector:    560 K, 1300 ppm boron, 0.7533 g/cm3
  - All control rods withdrawn (ARO)

Assembly Types
--------------
  UOX (Type 1): Graded uranium oxide with U-Gd burnable absorber rods
    - U_4.2: 4.2 wt% U-235 (interior fuel pins)
    - U_3.7: 3.7 wt% U-235 (peripheral fuel pins)
    - TVEG_5: 3.3 wt% U-235 + 5 wt% Gd2O3 (12 BA rods)
    - Burnup levels: 0, 15, 32, 40 MWd/kg

  MOX (Type 2): Graded, profiled MOX with U-Gd burnable absorber rods
    - PU_3.6: 3.62 wt% fissile Pu (outer fuel zone)
    - PU_2.7: 2.69 wt% fissile Pu (middle fuel zone)
    - PU_2.4: 2.42 wt% fissile Pu (inner fuel zone)
    - TVEG_4: 3.6 wt% U-235 + 4 wt% Gd2O3 (12 BA rods)
    - Burnup levels: 0, 17, 33 MWd/kg

Pin Geometry (cm)
-----------------
  Fuel pellet outer radius:  0.386  (no central hole, no gas gap)
  Cladding outer radius:     0.455
  Pin pitch (hex):           1.275
  Guide tube inner radius:   0.55
  Guide tube outer radius:   0.63
  Assembly pitch (hex F2F):  23.6

Reflector Structure (from center outward)
-----------------------------------------
  Water gap:      3 mm between outermost FA row and steel buffer
  Steel buffer:   R = 175.0 to 181.0 cm (with water holes, modeled homogenized)
  Steel barrel:   R = 181.0 to 206.8 cm
  Downcomer:      R = 206.8 to 226.75 cm (water)
  Steel vessel:   beyond R = 226.75 cm (vacuum boundary)

All material compositions use exact number densities from the benchmark
specification (Tables A.3-A.11).

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# FUEL COMPOSITIONS (atoms/barn*cm) from Tables A.3-A.9
# =============================================================================
# Each fuel type is a dict mapping burnup (MWd/kg) to {nuclide: number_density}.
# Only nuclides with nonzero density at that burnup are included.

FUEL_U42 = {
    0: {
        'U235': 9.0411E-04, 'U238': 2.0362E-02, 'O16': 4.2532E-02,
    },
    15: {
        'U235': 5.8139E-04, 'U236': 5.7700E-05, 'U238': 2.0161E-02,
        'Np237': 3.7658E-06, 'Pu238': 4.5135E-07, 'Pu239': 9.6584E-05,
        'Pu240': 1.7820E-05, 'Pu241': 7.8291E-06, 'Pu242': 8.5918E-07,
        'Am241': 9.6169E-08, 'O16': 4.2532E-02,
        'Sm149': 9.1807E-08, 'Sm151': 3.8524E-07,
        'Tc99': 1.9865E-05, 'Rh103': 1.1241E-05, 'Cs133': 2.1660E-05,
        'Nd143': 1.7145E-05, 'Nd145': 1.2231E-05, 'Pm147': 5.1135E-06,
        'Sm152': 2.0001E-06,
    },
    32: {
        'U235': 3.2990E-04, 'U236': 9.7452E-05, 'U238': 1.9905E-02,
        'Np237': 1.0351E-05, 'Pu238': 2.7546E-06, 'Pu239': 1.2788E-04,
        'Pu240': 4.4214E-05, 'Pu241': 2.4361E-05, 'Pu242': 6.9401E-06,
        'Am241': 5.7445E-07, 'O16': 4.2532E-02,
        'Sm149': 8.9565E-08, 'Sm151': 4.9458E-07,
        'Tc99': 4.0006E-05, 'Rh103': 2.2541E-05, 'Cs133': 4.3134E-05,
        'Nd143': 3.0721E-05, 'Nd145': 2.3856E-05, 'Pm147': 7.0579E-06,
        'Sm152': 4.0165E-06,
    },
    40: {
        'U235': 2.4314E-04, 'U236': 1.0890E-04, 'U238': 1.9773E-02,
        'Np237': 1.3577E-05, 'Pu238': 4.6336E-06, 'Pu239': 1.3134E-04,
        'Pu240': 5.4717E-05, 'Pu241': 3.0501E-05, 'Pu242': 1.1811E-05,
        'Am241': 8.2251E-07, 'O16': 4.2532E-02,
        'Sm149': 8.5421E-08, 'Sm151': 5.3181E-07,
        'Tc99': 4.8545E-05, 'Rh103': 2.6925E-05, 'Cs133': 5.2008E-05,
        'Nd143': 3.5092E-05, 'Nd145': 2.8579E-05, 'Pm147': 7.2585E-06,
        'Sm152': 4.7868E-06,
    },
}

FUEL_TVEG5 = {
    0: {
        'U235': 6.6163E-04, 'U238': 1.9143E-02, 'O16': 4.1938E-02,
        'Gd152': 3.2142E-06, 'Gd154': 3.4579E-05, 'Gd155': 2.3321E-04,
        'Gd156': 3.2053E-04, 'Gd157': 2.4346E-04, 'Gd158': 3.8403E-04,
        'Gd160': 3.3373E-04,
    },
    15: {
        'U235': 4.8382E-04, 'U236': 3.5164E-05, 'U238': 1.8968E-02,
        'Np237': 2.5863E-06, 'Pu238': 2.9960E-07, 'Pu239': 8.8781E-05,
        'Pu240': 1.5353E-05, 'Pu241': 6.2524E-06, 'Pu242': 6.0179E-07,
        'Am241': 6.9561E-08, 'O16': 4.1938E-02,
        'Gd152': 2.2242E-06, 'Gd154': 3.1600E-05, 'Gd155': 3.2385E-07,
        'Gd156': 5.4422E-04, 'Gd157': 1.8330E-07, 'Gd158': 6.3019E-04,
        'Gd160': 3.3229E-04,
        'Sm149': 8.0821E-08, 'Sm151': 3.1290E-07,
        'Tc99': 1.2031E-05, 'Rh103': 7.4977E-06, 'Cs133': 1.3175E-05,
        'Nd143': 1.0380E-05, 'Nd145': 7.3069E-06, 'Pm147': 3.3064E-06,
        'Sm152': 1.2690E-06,
    },
    32: {
        'U235': 2.6776E-04, 'U236': 7.0094E-05, 'U238': 1.8728E-02,
        'Np237': 7.7190E-06, 'Pu238': 2.0282E-06, 'Pu239': 1.1559E-04,
        'Pu240': 4.1294E-05, 'Pu241': 2.2283E-05, 'Pu242': 6.3481E-06,
        'Am241': 4.9976E-07, 'O16': 4.1938E-02,
        'Gd152': 1.1130E-06, 'Gd154': 2.7446E-05, 'Gd155': 1.5769E-07,
        'Gd156': 5.3058E-04, 'Gd157': 1.6774E-07, 'Gd158': 6.3494E-04,
        'Gd160': 3.3037E-04,
        'Sm149': 8.0584E-08, 'Sm151': 4.0538E-07,
        'Tc99': 3.0625E-05, 'Rh103': 1.8625E-05, 'Cs133': 3.3172E-05,
        'Nd143': 2.3428E-05, 'Nd145': 1.8017E-05, 'Pm147': 5.7735E-06,
        'Sm152': 3.2490E-06,
    },
    40: {
        'U235': 1.9449E-04, 'U236': 8.0032E-05, 'U238': 1.8603E-02,
        'Np237': 1.0347E-05, 'Pu238': 3.5111E-06, 'Pu239': 1.1801E-04,
        'Pu240': 5.1119E-05, 'Pu241': 2.8043E-05, 'Pu242': 1.1103E-05,
        'Am241': 7.2654E-07, 'O16': 4.1938E-02,
        'Gd152': 7.6549E-07, 'Gd154': 2.5465E-05, 'Gd155': 1.4181E-07,
        'Gd156': 5.2333E-04, 'Gd157': 1.5847E-07, 'Gd158': 6.3726E-04,
        'Gd160': 3.2938E-04,
        'Sm149': 7.7330E-08, 'Sm151': 4.4104E-07,
        'Tc99': 3.8586E-05, 'Rh103': 2.2940E-05, 'Cs133': 4.1526E-05,
        'Nd143': 2.7732E-05, 'Nd145': 2.2419E-05, 'Pm147': 6.1251E-06,
        'Sm152': 3.9837E-06,
    },
}

FUEL_U37 = {
    0: {
        'U235': 7.9649E-04, 'U238': 2.0469E-02, 'O16': 4.2530E-02,
    },
    15: {
        'U235': 4.8884E-04, 'U236': 5.4042E-05, 'U238': 2.0262E-02,
        'Np237': 3.7015E-06, 'Pu238': 4.6522E-07, 'Pu239': 9.5675E-05,
        'Pu240': 1.9149E-05, 'Pu241': 8.3590E-06, 'Pu242': 1.0147E-06,
        'Am241': 1.0325E-07, 'O16': 4.2530E-02,
        'Sm149': 8.3005E-08, 'Sm151': 3.4935E-07,
        'Tc99': 1.9485E-05, 'Rh103': 1.1189E-05, 'Cs133': 2.1245E-05,
        'Nd143': 1.6565E-05, 'Nd145': 1.1935E-05, 'Pm147': 4.9515E-06,
        'Sm152': 2.0058E-06,
    },
    32: {
        'U235': 2.6496E-04, 'U236': 8.8494E-05, 'U238': 2.0000E-02,
        'Np237': 9.9287E-06, 'Pu238': 2.7516E-06, 'Pu239': 1.2378E-04,
        'Pu240': 4.5926E-05, 'Pu241': 2.4695E-05, 'Pu242': 7.7274E-06,
        'Am241': 5.7739E-07, 'O16': 4.2530E-02,
        'Sm149': 8.2706E-08, 'Sm151': 4.5427E-07,
        'Tc99': 3.8798E-05, 'Rh103': 2.2235E-05, 'Cs133': 4.1820E-05,
        'Nd143': 2.9007E-05, 'Nd145': 2.2961E-05, 'Pm147': 6.7031E-06,
        'Sm152': 3.9584E-06,
    },
    40: {
        'U235': 1.9092E-04, 'U236': 9.7726E-05, 'U238': 1.9864E-02,
        'Np237': 1.2901E-05, 'Pu238': 4.5708E-06, 'Pu239': 1.2659E-04,
        'Pu240': 5.6139E-05, 'Pu241': 3.0494E-05, 'Pu242': 1.2925E-05,
        'Am241': 8.0933E-07, 'O16': 4.2530E-02,
        'Sm149': 7.9571E-08, 'Sm151': 4.9115E-07,
        'Tc99': 4.6961E-05, 'Rh103': 2.6478E-05, 'Cs133': 5.0289E-05,
        'Nd143': 3.2889E-05, 'Nd145': 2.7419E-05, 'Pm147': 6.8639E-06,
        'Sm152': 4.6994E-06,
    },
}

FUEL_PU36 = {
    0: {
        'U235': 4.3057E-05, 'U238': 2.0386E-02,
        'Pu238': 1.0841E-06, 'Pu239': 7.5661E-04, 'Pu240': 5.3794E-05,
        'Pu241': 9.5720E-06, 'Pu242': 3.5119E-06,
        'O16': 4.2506E-02,
    },
    17: {
        'U235': 3.0534E-05, 'U236': 2.5385E-06, 'U238': 2.0144E-02,
        'Np237': 2.4045E-06, 'Pu238': 1.3292E-06, 'Pu239': 4.7406E-04,
        'Pu240': 1.4795E-04, 'Pu241': 5.7132E-05, 'Pu242': 9.8236E-06,
        'Am241': 1.3594E-06, 'O16': 4.2506E-02,
        'Sm149': 1.4783E-07, 'Sm151': 7.7056E-07,
        'Tc99': 2.2348E-05, 'Rh103': 2.2129E-05, 'Cs133': 2.4904E-05,
        'Nd143': 1.5770E-05, 'Nd145': 1.1215E-05, 'Pm147': 5.0355E-06,
        'Sm152': 3.2199E-06,
    },
    33: {
        'U235': 2.0186E-05, 'U236': 4.2696E-06, 'U238': 1.9894E-02,
        'Np237': 4.2797E-06, 'Pu238': 2.7311E-06, 'Pu239': 2.9852E-04,
        'Pu240': 1.7846E-04, 'Pu241': 8.3282E-05, 'Pu242': 2.5860E-05,
        'Am241': 2.9942E-06, 'O16': 4.2506E-02,
        'Sm149': 1.2062E-07, 'Sm151': 7.6447E-07,
        'Tc99': 4.0862E-05, 'Rh103': 3.6232E-05, 'Cs133': 4.4872E-05,
        'Nd143': 2.7603E-05, 'Nd145': 2.0679E-05, 'Pm147': 6.6403E-06,
        'Sm152': 5.2658E-06,
    },
}

FUEL_TVEG4 = {
    0: {
        'U235': 7.3225E-04, 'U238': 1.9360E-02, 'O16': 4.2056E-02,
        'Gd152': 2.5815E-06, 'Gd154': 2.7772E-05, 'Gd155': 1.8730E-04,
        'Gd156': 2.5743E-04, 'Gd157': 1.9553E-04, 'Gd158': 3.0843E-04,
        'Gd160': 2.6804E-04,
    },
    17: {
        'U235': 5.4783E-04, 'U236': 3.9989E-05, 'U238': 1.9139E-02,
        'Np237': 3.7973E-06, 'Pu238': 4.9031E-07, 'Pu239': 1.2109E-04,
        'Pu240': 1.7377E-05, 'Pu241': 7.0208E-06, 'Pu242': 5.4476E-07,
        'Am241': 8.8029E-08, 'O16': 4.2055E-02,
        'Gd152': 1.8321E-06, 'Gd154': 2.5080E-05, 'Gd155': 1.0215E-06,
        'Gd156': 4.3334E-04, 'Gd157': 2.5976E-07, 'Gd158': 5.0718E-04,
        'Gd160': 2.6656E-04,
        'Sm149': 1.0735E-07, 'Sm151': 4.0519E-07,
        'Tc99': 1.2602E-05, 'Rh103': 7.9971E-06, 'Cs133': 1.3779E-05,
        'Nd143': 1.0973E-05, 'Nd145': 7.6416E-06, 'Pm147': 3.2794E-06,
        'Sm152': 1.2631E-06,
    },
    33: {
        'U235': 3.4998E-04, 'U236': 7.3537E-05, 'U238': 1.8901E-02,
        'Np237': 9.6340E-06, 'Pu238': 2.5298E-06, 'Pu239': 1.5184E-04,
        'Pu240': 4.4993E-05, 'Pu241': 2.4072E-05, 'Pu242': 5.0595E-06,
        'Am241': 5.5623E-07, 'O16': 4.2055E-02,
        'Gd152': 1.0821E-06, 'Gd154': 2.2136E-05, 'Gd155': 1.6393E-07,
        'Gd156': 4.2171E-04, 'Gd157': 2.0690E-07, 'Gd158': 5.1162E-04,
        'Gd160': 2.6501E-04,
        'Sm149': 1.0224E-07, 'Sm151': 5.0972E-07,
        'Tc99': 2.9713E-05, 'Rh103': 1.8557E-05, 'Cs133': 3.2177E-05,
        'Nd143': 2.3818E-05, 'Nd145': 1.7521E-05, 'Pm147': 5.5753E-06,
        'Sm152': 3.0576E-06,
    },
}

FUEL_PU27 = {
    0: {
        'U235': 4.3057E-05, 'U238': 2.0598E-02,
        'Pu238': 8.0774E-07, 'Pu239': 5.6222E-04, 'Pu240': 3.9987E-05,
        'Pu241': 7.1160E-06, 'Pu242': 2.6131E-06,
        'O16': 4.2508E-02,
    },
    17: {
        'U235': 2.8612E-05, 'U236': 2.7944E-06, 'U238': 2.0347E-02,
        'Np237': 2.4008E-06, 'Pu238': 1.0993E-06, 'Pu239': 3.3450E-04,
        'Pu240': 1.2215E-04, 'Pu241': 4.9442E-05, 'Pu242': 9.5258E-06,
        'Am241': 1.1240E-06, 'O16': 4.2508E-02,
        'Sm149': 1.1354E-07, 'Sm151': 5.8719E-07,
        'Tc99': 2.0699E-05, 'Rh103': 2.0188E-05, 'Cs133': 2.3035E-05,
        'Nd143': 1.4448E-05, 'Nd145': 1.0429E-05, 'Pm147': 4.6128E-06,
        'Sm152': 3.0318E-06,
    },
    33: {
        'U235': 1.7606E-05, 'U236': 4.5372E-06, 'U238': 2.0087E-02,
        'Np237': 4.2361E-06, 'Pu238': 2.5148E-06, 'Pu239': 2.1853E-04,
        'Pu240': 1.4136E-04, 'Pu241': 6.9024E-05, 'Pu242': 2.5765E-05,
        'Am241': 2.3584E-06, 'O16': 4.2508E-02,
        'Sm149': 9.8177E-08, 'Sm151': 6.0027E-07,
        'Tc99': 3.7403E-05, 'Rh103': 3.2327E-05, 'Cs133': 4.0983E-05,
        'Nd143': 2.4621E-05, 'Nd145': 1.9007E-05, 'Pm147': 5.9544E-06,
        'Sm152': 4.7872E-06,
    },
}

FUEL_PU24 = {
    0: {
        'U235': 4.3057E-05, 'U238': 2.0660E-02,
        'Pu238': 7.2271E-07, 'Pu239': 5.0579E-04, 'Pu240': 3.5961E-05,
        'Pu241': 6.4023E-06, 'Pu242': 2.3413E-06,
        'O16': 4.2508E-02,
    },
    17: {
        'U235': 2.7777E-05, 'U236': 2.9076E-06, 'U238': 2.0405E-02,
        'Np237': 2.3897E-06, 'Pu238': 1.0335E-06, 'Pu239': 2.9332E-04,
        'Pu240': 1.1497E-04, 'Pu241': 4.6985E-05, 'Pu242': 9.6174E-06,
        'Am241': 1.0464E-06, 'O16': 4.2508E-02,
        'Sm149': 1.0281E-07, 'Sm151': 5.3010E-07,
        'Tc99': 2.0326E-05, 'Rh103': 1.9698E-05, 'Cs133': 2.2608E-05,
        'Nd143': 1.4108E-05, 'Nd145': 1.0254E-05, 'Pm147': 4.5137E-06,
        'Sm152': 2.9973E-06,
    },
    33: {
        'U235': 1.6391E-05, 'U236': 4.6650E-06, 'U238': 2.0140E-02,
        'Np237': 4.1976E-06, 'Pu238': 2.4576E-06, 'Pu239': 1.9329E-04,
        'Pu240': 1.3072E-04, 'Pu241': 6.4083E-05, 'Pu242': 2.6319E-05,
        'Am241': 2.1281E-06, 'O16': 4.2508E-02,
        'Sm149': 9.0532E-08, 'Sm151': 5.4672E-07,
        'Tc99': 3.6699E-05, 'Rh103': 3.1334E-05, 'Cs133': 4.0174E-05,
        'Nd143': 2.3798E-05, 'Nd145': 1.8670E-05, 'Pm147': 5.7880E-06,
        'Sm152': 4.6882E-06,
    },
}

# =============================================================================
# STRUCTURAL MATERIAL COMPOSITIONS (Table A.10, atoms/barn*cm)
# =============================================================================

ZIRCALLOY = {  # Fuel cladding, central tube, guide tube
    'Zr90': 4.257E-02 * 0.5145,   # Natural Zr isotopic fractions
    'Zr91': 4.257E-02 * 0.1122,
    'Zr92': 4.257E-02 * 0.1715,
    'Zr94': 4.257E-02 * 0.1738,
    'Zr96': 4.257E-02 * 0.0280,
    'Nb93': 4.223E-04,
    'Hf180': 6.594E-06 * 0.3508,  # Dominant Hf isotope
    'Hf178': 6.594E-06 * 0.2728,
    'Hf177': 6.594E-06 * 0.1860,
    'Hf179': 6.594E-06 * 0.1362,
    'Hf176': 6.594E-06 * 0.0526,
}

STEEL = {  # Steel buffer, barrel, vessel
    'Fe56': 5.933E-02 * 0.9175,
    'Fe54': 5.933E-02 * 0.0585,
    'Fe57': 5.933E-02 * 0.0212,
    'Fe58': 5.933E-02 * 0.0028,
    'Cr52': 1.687E-02 * 0.8379,
    'Cr53': 1.687E-02 * 0.0950,
    'Cr50': 1.687E-02 * 0.0435,
    'Cr54': 1.687E-02 * 0.0236,
    'Ni58': 8.477E-03 * 0.6808,
    'Ni60': 8.477E-03 * 0.2622,
    'Ni62': 8.477E-03 * 0.0364,
    'Ni61': 8.477E-03 * 0.0114,
    'Ni64': 8.477E-03 * 0.0093,
    'Ti48': 9.904E-04 * 0.7372,
    'Ti46': 9.904E-04 * 0.0825,
    'Ti47': 9.904E-04 * 0.0744,
    'Ti49': 9.904E-04 * 0.0541,
    'Ti50': 9.904E-04 * 0.0518,
    'C12': 4.737E-04 * 0.9893,
    'C13': 4.737E-04 * 0.0107,
}

# =============================================================================
# MODERATOR COMPOSITIONS (Table A.11, atoms/barn*cm)
# =============================================================================

# M575B1.3: Moderator in FA, 575 K, 1300 ppm boron, 0.7241 g/cm3
MOD_575_B13 = {
    'H1': 4.8410E-02, 'O16': 2.4205E-02,
    'B10': 1.0381E-05, 'B11': 4.2049E-05,
}

# M560B1.3: Water in reflector, 560 K, 1300 ppm boron, 0.7533 g/cm3
MOD_560_B13 = {
    'H1': 5.0362E-02, 'O16': 2.5181E-02,
    'B10': 1.0800E-05, 'B11': 4.3744E-05,
}


# =============================================================================
# GUIDE TUBE POSITIONS WITHIN ASSEMBLY (OpenMC HexLattice indexing)
# =============================================================================
# 18 guide tube positions in standard VVER-1000 layout.
# Format: (ring_index_from_outer, position_in_ring) for 11-ring HexLattice.

GUIDE_TUBE_POSITIONS = [
    # Physical ring 7 (ring_index 3 from outer in 11-ring lattice)
    (3, 2), (3, 5), (3, 8), (3, 11), (3, 14), (3, 17),
    # Physical ring 6 (ring_index 4)
    (4, 3), (4, 9), (4, 15), (4, 21), (4, 27), (4, 33),
    # Physical ring 5 (ring_index 5)
    (5, 0), (5, 5), (5, 10), (5, 15), (5, 20), (5, 25),
]

# 12 Gd burnable absorber rod positions (approximate standard VVER-1000 layout)
# Placed in physical rings 4 and 8 with 6-fold symmetry (2 per sector)
GD_ROD_POSITIONS = [
    # Physical ring 4 (ring_index 6): 6 positions evenly spaced
    (6, 0), (6, 4), (6, 8), (6, 12), (6, 16), (6, 20),
    # Physical ring 8 (ring_index 2): 6 positions evenly spaced
    (2, 0), (2, 8), (2, 16), (2, 24), (2, 32), (2, 40),
]


# =============================================================================
# 60-DEGREE SECTOR CORE MAP (Figure A.1)
# =============================================================================
# 28 assemblies in the symmetry sector, numbered 1-28.
# Each entry: (assembly_number, type, burnup_MWd_per_kg)
# Type 1 = UOX, Type 2 = MOX

SECTOR_MAP = {
    # Ring 0 (center): 1 assembly
    0: [(1, 1, 40)],
    # Ring 1: 1 assembly in sector
    1: [(2, 1, 32)],
    # Ring 2: 2 assemblies in sector
    2: [(3, 2, 0), (4, 2, 17)],
    # Ring 3: 3 assemblies
    3: [(5, 2, 33), (6, 1, 32), (7, 2, 0)],
    # Ring 4: 4 assemblies
    4: [(8, 1, 32), (9, 1, 15), (10, 2, 33), (11, 2, 33)],
    # Ring 5: 5 assemblies
    5: [(12, 1, 32), (13, 1, 15), (14, 1, 15), (15, 1, 32), (16, 2, 0)],
    # Ring 6: 6 assemblies
    6: [(17, 1, 0), (18, 2, 17), (19, 1, 15), (20, 2, 17),
        (21, 1, 0), (22, 1, 40)],
    # Ring 7: 6 fuel assemblies + 1 reflector (None)
    7: [(23, 1, 40), (24, 1, 0), (25, 1, 0), (26, 1, 15),
        (27, 1, 0), (28, 1, 40), None],
}


def _expand_sector_to_full_core():
    """Expand 60-degree sector to full core using rotational symmetry.

    Returns a dict mapping ring -> list of (type, burnup) tuples for all
    positions in the full core. For reflector positions, returns None.
    """
    full_core = {}

    # Ring 0: center assembly (not rotated)
    assy = SECTOR_MAP[0][0]
    full_core[0] = [(assy[1], assy[2])]

    # Rings 1-7: each sector position maps to 6 positions in full core
    for ring in range(1, 8):
        sector_positions = SECTOR_MAP[ring]
        n_sector = len(sector_positions)
        n_full = 6 * ring
        full_ring = [None] * n_full

        for k, entry in enumerate(sector_positions):
            if entry is None:
                # Reflector position
                for s in range(6):
                    idx = k + s * n_sector
                    if idx < n_full:
                        full_ring[idx] = None
            else:
                _, assy_type, burnup = entry
                for s in range(6):
                    idx = k + s * n_sector
                    if idx < n_full:
                        full_ring[idx] = (assy_type, burnup)

        full_core[ring] = full_ring

    return full_core


def build_model(particles=15000, batches=200, inactive=50):
    """Build the full-core VVER-1000 MOX benchmark OpenMC model (State S1)."""

    model = openmc.Model()

    # =========================================================================
    # MATERIALS
    # =========================================================================

    def make_fuel_material(name, composition):
        """Create a fuel material from number density dict."""
        mat = openmc.Material(name=name)
        mat.set_density('sum')
        for nuclide, density in composition.items():
            mat.add_nuclide(nuclide, density)
        mat.temperature = 1027.0  # State S1 fuel temperature
        return mat

    def make_nonfuel_material(name, composition, temperature):
        """Create a structural/moderator material from number density dict."""
        mat = openmc.Material(name=name)
        mat.set_density('sum')
        for nuclide, density in composition.items():
            mat.add_nuclide(nuclide, density)
        mat.temperature = temperature
        return mat

    # --- Fuel materials: one per (fuel_type, burnup) combination ---
    # UOX assembly fuels
    fuel_u42 = {}
    fuel_tveg5 = {}
    fuel_u37 = {}
    for bu in [0, 15, 32, 40]:
        fuel_u42[bu] = make_fuel_material(f'U_4.2 BU={bu}', FUEL_U42[bu])
        fuel_tveg5[bu] = make_fuel_material(f'TVEG_5 BU={bu}', FUEL_TVEG5[bu])
        fuel_u37[bu] = make_fuel_material(f'U_3.7 BU={bu}', FUEL_U37[bu])

    # MOX assembly fuels
    fuel_pu36 = {}
    fuel_tveg4 = {}
    fuel_pu27 = {}
    fuel_pu24 = {}
    for bu in [0, 17, 33]:
        fuel_pu36[bu] = make_fuel_material(f'PU_3.6 BU={bu}', FUEL_PU36[bu])
        fuel_tveg4[bu] = make_fuel_material(f'TVEG_4 BU={bu}', FUEL_TVEG4[bu])
        fuel_pu27[bu] = make_fuel_material(f'PU_2.7 BU={bu}', FUEL_PU27[bu])
        fuel_pu24[bu] = make_fuel_material(f'PU_2.4 BU={bu}', FUEL_PU24[bu])

    # --- Structural materials ---
    zircalloy = make_nonfuel_material('Zr Alloy', ZIRCALLOY, 575.0)

    steel = make_nonfuel_material('Steel', STEEL, 560.0)

    # --- Moderator in FA (M575B1.3) ---
    mod_fa = make_nonfuel_material('Moderator FA', MOD_575_B13, 575.0)
    mod_fa.add_s_alpha_beta('c_H_in_H2O')

    # --- Water in reflector (M560B1.3) ---
    mod_refl = make_nonfuel_material('Water Reflector', MOD_560_B13, 560.0)
    mod_refl.add_s_alpha_beta('c_H_in_H2O')

    # Collect all materials
    all_mats = []
    for d in [fuel_u42, fuel_tveg5, fuel_u37,
              fuel_pu36, fuel_tveg4, fuel_pu27, fuel_pu24]:
        all_mats.extend(d.values())
    all_mats.extend([zircalloy, steel, mod_fa, mod_refl])

    model.materials = openmc.Materials(all_mats)
    model.materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # PIN CELL GEOMETRY
    # =========================================================================
    # Benchmark specification (Table A.1):
    #   Fuel cell: pellet R = 0.386 cm, clad OR = 0.455 cm
    #   No gas gap (clad IR = pellet OR), no central hole
    #   GT/CT: tube IR = 0.55 cm, tube OR = 0.63 cm

    fuel_or = 0.386
    clad_or = 0.455
    gt_ir = 0.55
    gt_or = 0.63
    pin_pitch = 1.275
    assembly_pitch = 23.6
    n_pin_rings = 11  # rings 0-10, total 331 positions

    fuel_surf = openmc.ZCylinder(r=fuel_or)
    clad_surf = openmc.ZCylinder(r=clad_or)
    gt_ir_surf = openmc.ZCylinder(r=gt_ir)
    gt_or_surf = openmc.ZCylinder(r=gt_or)

    # --- Fuel pin universe factory ---
    def make_fuel_pin(fuel_mat, name):
        """Fuel pin: solid pellet | clad | moderator (no gap, no hole)."""
        c_fuel = openmc.Cell(name=f'{name} pellet', fill=fuel_mat)
        c_fuel.region = -fuel_surf
        c_clad = openmc.Cell(name=f'{name} clad', fill=zircalloy)
        c_clad.region = +fuel_surf & -clad_surf
        c_mod = openmc.Cell(name=f'{name} mod', fill=mod_fa)
        c_mod.region = +clad_surf
        return openmc.Universe(name=name, cells=[c_fuel, c_clad, c_mod])

    # --- Guide tube / central tube universe ---
    gt_inner = openmc.Cell(name='GT inner', fill=mod_fa, region=-gt_ir_surf)
    gt_wall = openmc.Cell(name='GT wall', fill=zircalloy,
                          region=+gt_ir_surf & -gt_or_surf)
    gt_outer = openmc.Cell(name='GT outer', fill=mod_fa, region=+gt_or_surf)
    guide_tube_univ = openmc.Universe(
        name='Guide Tube', cells=[gt_inner, gt_wall, gt_outer])

    ct_inner = openmc.Cell(name='CT inner', fill=mod_fa, region=-gt_ir_surf)
    ct_wall = openmc.Cell(name='CT wall', fill=zircalloy,
                          region=+gt_ir_surf & -gt_or_surf)
    ct_outer = openmc.Cell(name='CT outer', fill=mod_fa, region=+gt_or_surf)
    instrument_tube_univ = openmc.Universe(
        name='Central Tube', cells=[ct_inner, ct_wall, ct_outer])

    # --- Moderator-only universe (for positions outside assembly hex) ---
    mod_cell_outer = openmc.Cell(name='Pin outer mod', fill=mod_fa)
    pin_outer_univ = openmc.Universe(name='Pin Outer', cells=[mod_cell_outer])

    # =========================================================================
    # ASSEMBLY LATTICE BUILDER
    # =========================================================================
    # UOX assembly zones:
    #   Rings 1-7: U_4.2 (type 2, interior fuel)
    #   Rings 8-10: U_3.7 (type 5, peripheral fuel)
    #   12 Gd rods: TVEG_5 (type 4) at GD_ROD_POSITIONS
    #   18 guide tubes at GUIDE_TUBE_POSITIONS
    #   Center: central tube
    #
    # MOX assembly zones (profiled: lower Pu inside, higher outside):
    #   Rings 1-3: PU_2.4 (type 6, lowest Pu, inner zone)
    #   Rings 4-7: PU_2.7 (type 5, middle Pu)
    #   Rings 8-10: PU_3.6 (type 2, highest Pu, outer zone)
    #   12 Gd rods: TVEG_4 (type 4) at GD_ROD_POSITIONS
    #   18 guide tubes at GUIDE_TUBE_POSITIONS
    #   Center: central tube

    def build_uox_assembly(burnup, name):
        """Build a UOX assembly at the given burnup."""
        pin_u42 = make_fuel_pin(fuel_u42[burnup], f'U_4.2@{burnup}')
        pin_u37 = make_fuel_pin(fuel_u37[burnup], f'U_3.7@{burnup}')
        pin_gd = make_fuel_pin(fuel_tveg5[burnup], f'TVEG_5@{burnup}')

        universes = []
        for ring_idx in range(n_pin_rings):
            physical_ring = n_pin_rings - 1 - ring_idx
            if physical_ring == 0:
                universes.append([instrument_tube_univ])
            else:
                # Zone assignment
                if physical_ring >= 8:
                    pin = pin_u37  # peripheral
                else:
                    pin = pin_u42  # interior
                n_pos = 6 * physical_ring
                universes.append([pin] * n_pos)

        # Place guide tubes
        for ri, pi in GUIDE_TUBE_POSITIONS:
            universes[ri][pi] = guide_tube_univ

        # Place Gd rods
        for ri, pi in GD_ROD_POSITIONS:
            universes[ri][pi] = pin_gd

        lattice = openmc.HexLattice(name=f'{name} pins')
        lattice.orientation = 'x'
        lattice.center = (0.0, 0.0)
        lattice.pitch = [pin_pitch]
        lattice.universes = universes
        lattice.outer = pin_outer_univ

        # Wrap in hex prism
        edge_length = assembly_pitch / np.sqrt(3.0)
        hex_prism = openmc.model.HexagonalPrism(
            edge_length=edge_length, origin=(0.0, 0.0), orientation='x')

        assy_cell = openmc.Cell(name=f'{name} lattice',
                                fill=lattice, region=-hex_prism)
        gap_cell = openmc.Cell(name=f'{name} gap',
                               fill=mod_fa, region=+hex_prism)
        return openmc.Universe(name=name, cells=[assy_cell, gap_cell])

    def build_mox_assembly(burnup, name):
        """Build a MOX assembly at the given burnup."""
        pin_pu36 = make_fuel_pin(fuel_pu36[burnup], f'PU_3.6@{burnup}')
        pin_pu27 = make_fuel_pin(fuel_pu27[burnup], f'PU_2.7@{burnup}')
        pin_pu24 = make_fuel_pin(fuel_pu24[burnup], f'PU_2.4@{burnup}')
        pin_gd = make_fuel_pin(fuel_tveg4[burnup], f'TVEG_4@{burnup}')

        universes = []
        for ring_idx in range(n_pin_rings):
            physical_ring = n_pin_rings - 1 - ring_idx
            if physical_ring == 0:
                universes.append([instrument_tube_univ])
            else:
                # Profiled MOX: inner=PU_2.4, middle=PU_2.7, outer=PU_3.6
                if physical_ring <= 3:
                    pin = pin_pu24
                elif physical_ring <= 7:
                    pin = pin_pu27
                else:
                    pin = pin_pu36
                n_pos = 6 * physical_ring
                universes.append([pin] * n_pos)

        # Place guide tubes
        for ri, pi in GUIDE_TUBE_POSITIONS:
            universes[ri][pi] = guide_tube_univ

        # Place Gd rods
        for ri, pi in GD_ROD_POSITIONS:
            universes[ri][pi] = pin_gd

        lattice = openmc.HexLattice(name=f'{name} pins')
        lattice.orientation = 'x'
        lattice.center = (0.0, 0.0)
        lattice.pitch = [pin_pitch]
        lattice.universes = universes
        lattice.outer = pin_outer_univ

        edge_length = assembly_pitch / np.sqrt(3.0)
        hex_prism = openmc.model.HexagonalPrism(
            edge_length=edge_length, origin=(0.0, 0.0), orientation='x')

        assy_cell = openmc.Cell(name=f'{name} lattice',
                                fill=lattice, region=-hex_prism)
        gap_cell = openmc.Cell(name=f'{name} gap',
                               fill=mod_fa, region=+hex_prism)
        return openmc.Universe(name=name, cells=[assy_cell, gap_cell])

    # Build all 7 unique assembly types
    assy_univs = {}

    # UOX at burnups 0, 15, 32, 40
    for bu in [0, 15, 32, 40]:
        key = (1, bu)
        assy_univs[key] = build_uox_assembly(bu, f'UOX BU={bu}')

    # MOX at burnups 0, 17, 33
    for bu in [0, 17, 33]:
        key = (2, bu)
        assy_univs[key] = build_mox_assembly(bu, f'MOX BU={bu}')

    # Reflector universe (water)
    refl_cell = openmc.Cell(name='reflector water', fill=mod_refl)
    reflector_univ = openmc.Universe(name='Reflector', cells=[refl_cell])

    # =========================================================================
    # CORE-LEVEL HEXAGONAL LATTICE
    # =========================================================================

    full_core = _expand_sector_to_full_core()
    n_core_rings = 8  # rings 0-7

    core_universes = []
    n_fuel = 0
    n_mox = 0

    for ring_idx in range(n_core_rings):
        physical_ring = n_core_rings - 1 - ring_idx
        ring_data = full_core[physical_ring]
        ring_univs = []

        for entry in ring_data:
            if entry is None:
                ring_univs.append(reflector_univ)
            else:
                assy_type, burnup = entry
                ring_univs.append(assy_univs[(assy_type, burnup)])
                n_fuel += 1
                if assy_type == 2:
                    n_mox += 1

        core_universes.append(ring_univs)

    core_lattice = openmc.HexLattice(name='VVER-1000 Core')
    core_lattice.orientation = 'x'
    core_lattice.center = (0.0, 0.0)
    core_lattice.pitch = [assembly_pitch]
    core_lattice.universes = core_universes
    core_lattice.outer = reflector_univ

    assert n_fuel == 163, f"Expected 163 fuel assemblies, got {n_fuel}"
    mox_pct = n_mox / n_fuel * 100.0

    # =========================================================================
    # REFLECTOR STRUCTURE AND BOUNDARIES
    # =========================================================================
    # From Figure A.2:
    #   Core edge → 175.0 cm: water gap (3 mm)
    #   175.0 → 181.0 cm: steel buffer (with water holes, modeled as steel)
    #   181.0 → 206.8 cm: steel barrel
    #   206.8 → 226.75 cm: downcomer (water)
    #   226.75 cm: outer vacuum boundary (steel vessel not explicitly modeled)

    core_cyl = openmc.ZCylinder(r=175.0, name='Core/buffer boundary')
    buffer_cyl = openmc.ZCylinder(r=181.0, name='Buffer/barrel boundary')
    barrel_cyl = openmc.ZCylinder(r=206.8, name='Barrel/downcomer boundary')
    vessel_cyl = openmc.ZCylinder(r=226.75, boundary_type='vacuum',
                                  name='Vessel OR (vacuum)')

    # Core region with hex lattice
    core_cell = openmc.Cell(name='Core', fill=core_lattice, region=-core_cyl)

    # Steel buffer (simplified: homogeneous steel, ignoring water holes)
    buffer_cell = openmc.Cell(name='Steel buffer', fill=steel,
                              region=+core_cyl & -buffer_cyl)

    # Steel barrel
    barrel_cell = openmc.Cell(name='Steel barrel', fill=steel,
                              region=+buffer_cyl & -barrel_cyl)

    # Downcomer (water)
    dc_cell = openmc.Cell(name='Downcomer', fill=mod_refl,
                          region=+barrel_cyl & -vessel_cyl)

    root = openmc.Universe(name='root')
    root.add_cells([core_cell, buffer_cell, barrel_cell, dc_cell])
    model.geometry = openmc.Geometry(root)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    # Initial source in core region
    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(
        [-170.0, -170.0, -1.0], [170.0, 170.0, 1.0])
    source.constraints = {'fissionable': True}
    settings.source = source

    # Temperature treatment
    settings.temperature = {'method': 'interpolation', 'default': 575.0}

    model.settings = settings

    return model


def main():
    parser = argparse.ArgumentParser(
        description='VVER-1000 MOX Core Benchmark (NEA/NSC/DOC(2005)17)')
    parser.add_argument('--particles', type=int, default=15000)
    parser.add_argument('--batches', type=int, default=200)
    parser.add_argument('--inactive', type=int, default=50)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()

    model = build_model(
        particles=args.particles,
        batches=args.batches,
        inactive=args.inactive,
    )

    # Count assemblies
    full_core = _expand_sector_to_full_core()
    n_total = 0
    n_mox = 0
    for ring_data in full_core.values():
        for entry in ring_data:
            if entry is not None:
                n_total += 1
                if entry[0] == 2:
                    n_mox += 1

    model.export_to_model_xml()
    print(f"VVER-1000 MOX Core Benchmark (NEA/NSC/DOC(2005)17)")
    print(f"  State S1: Working state, ARO, 1300 ppm boron")
    print(f"  Fuel T = 1027 K, Moderator T = 575 K")
    print(f"  Core: {n_total} assemblies ({n_mox} MOX = {n_mox/n_total*100:.1f}%)")
    print(f"  7 unique assembly types (4 UOX burnups + 3 MOX burnups)")
    print(f"  Assembly pitch: 23.6 cm, Pin pitch: 1.275 cm")
    print(f"  Particles: {args.particles}, Batches: {args.batches}, "
          f"Inactive: {args.inactive}")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
