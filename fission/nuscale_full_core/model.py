#!/usr/bin/env python3
"""
McSAFER NuScale-like SMR Full-Core Benchmark (37 Fuel Assemblies)
==================================================================

Fridman, E., Bilodid, Y., Valtavirta, V. (2023). "Definition of the
neutronics benchmark of the NuScale-like core." Nuclear Engineering and
Technology, 55(10), 3639-3647. DOI: 10.1016/j.net.2023.06.029

Bousquet, J. et al. (2025). "Neutronics Benchmark of the NuScale-like SMR
in the McSAFER Project." J. Nucl. Eng. 6(4), 44.

Number densities from RODARE dataset: DOI: 10.14278/rodare.2457

Core Layout (37 assemblies, 7x7 grid minus corners)
-----------------------------------------------------
  Row 1:             C01   B02   C01
  Row 2:       C02   B01   A01   B01   C02
  Row 3: C01   B01   A02   A01   A02   B01   C01
  Row 4: B02   A01   A01   C03   A01   A01   B02
  Row 5: C01   B01   A02   A01   A02   B01   C01
  Row 6:       C02   B01   A01   B01   C02
  Row 7:             C01   B02   C01

Assembly Types
--------------
  A01: 1.50 wt% U-235, 8 assemblies
  A02: 1.60 wt% U-235, 4 assemblies
  B01: 2.50 wt% U-235, 8 assemblies
  B02: 2.60 wt% U-235, 4 assemblies
  C01: 4.05 wt% U-235, 8 assemblies
  C02: 4.55 wt% U-235, 4 assemblies (16 Gd rods with 8 wt% Gd2O3)
  C03: 2.60 wt% U-235, 1 assembly (center)

Pin Dimensions (Fridman et al. / RODARE dataset)
--------------------------------------------------
  Fuel pellet OR:    0.4058 cm
  Cladding IR:       0.4140 cm
  Cladding OR:       0.4750 cm
  Guide tube IR:     0.5715 cm
  Guide tube OR:     0.6121 cm
  Instrument tube:   same as guide tube
  Pin pitch:         1.2598 cm
  Assembly pitch:    21.5036 cm
  Active fuel height: 200.0 cm

Operating Conditions (BOC, HFP benchmark state)
-------------------------------------------------
  Fuel temperature:     900 K
  Coolant temperature:  600 K (average)
  Boron concentration:  1000 ppm
  System pressure:      12.755 MPa
  Core power:           160 MWth

Usage:
    python model.py [--particles N] [--batches N] [--inactive N] [--run]
"""

import argparse
import numpy as np
import openmc


# =============================================================================
# PIN LAYOUT (17x17 Westinghouse-type)
# =============================================================================
GUIDE_TUBE_POSITIONS = {
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5),          (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
}

INSTRUMENT_TUBE_POSITION = (8, 8)

# 16 Gd rod positions for C02 assemblies (approximate symmetric pattern;
# exact positions from Framatome HTP fuel design not publicly available)
GD_ROD_POSITIONS = {
    (2, 2), (2, 14),
    (3, 6), (3, 10),
    (6, 3), (6, 13),
    (7, 7), (7, 9),
    (9, 7), (9, 9),
    (10, 3), (10, 13),
    (13, 6), (13, 10),
    (14, 2), (14, 14),
}

# =============================================================================
# CORE MAP (7x7, 0=water)
# =============================================================================
# Assembly type strings
CORE_MAP = [
    [None, None, 'C01', 'B02', 'C01', None, None],
    [None, 'C02', 'B01', 'A01', 'B01', 'C02', None],
    ['C01', 'B01', 'A02', 'A01', 'A02', 'B01', 'C01'],
    ['B02', 'A01', 'A01', 'C03', 'A01', 'A01', 'B02'],
    ['C01', 'B01', 'A02', 'A01', 'A02', 'B01', 'C01'],
    [None, 'C02', 'B01', 'A01', 'B01', 'C02', None],
    [None, None, 'C01', 'B02', 'C01', None, None],
]

_n_assy = sum(1 for row in CORE_MAP for v in row if v is not None)
assert _n_assy == 37, f"Expected 37, got {_n_assy}"

# =============================================================================
# FUEL COMPOSITIONS (number densities in atoms/barn-cm from RODARE dataset)
# =============================================================================
# Each fuel type has a name and nuclide number densities.
# C02 assemblies have two pin types: regular (C02) and Gd-bearing (C02_Gd).

FUEL_DATA = {
    'A01': {
        'name': 'UO2 1.50%',
        'nuclides': {
            'U235': 3.56456e-4, 'U238': 2.31116e-2, 'O16': 4.69361e-2,
        }
    },
    'A02': {
        'name': 'UO2 1.60%',
        'nuclides': {
            'U235': 3.80219e-4, 'U238': 2.30881e-2, 'O16': 4.69367e-2,
        }
    },
    'B01': {
        'name': 'UO2 2.50%',
        'nuclides': {
            'U235': 5.94084e-4, 'U238': 2.28766e-2, 'O16': 4.69414e-2,
        }
    },
    'B02': {
        'name': 'UO2 2.60%',
        'nuclides': {
            'U235': 6.17846e-4, 'U238': 2.28531e-2, 'O16': 4.69419e-2,
        }
    },
    'C01': {
        'name': 'UO2 4.05%',
        'nuclides': {
            'U235': 9.62391e-4, 'U238': 2.25123e-2, 'O16': 4.69495e-2,
        }
    },
    'C02': {
        'name': 'UO2 4.55%',
        'nuclides': {
            'U235': 1.0812e-3, 'U238': 2.23948e-2, 'O16': 4.69521e-2,
        }
    },
    'C03': {
        'name': 'UO2 2.60% (center)',
        'nuclides': {
            'U235': 6.17846e-4, 'U238': 2.28531e-2, 'O16': 4.69419e-2,
        }
    },
    'C02_Gd': {
        'name': 'UO2-Gd2O3 4.55% 8wt%Gd',
        'nuclides': {
            'U235': 9.61347e-4, 'U238': 1.99124e-2, 'O16': 4.58016e-2,
            'Gd152': 5.40541e-6, 'Gd154': 5.89189e-5,
            'Gd155': 4.00000e-4, 'Gd156': 5.53243e-4,
            'Gd157': 4.22973e-4, 'Gd158': 6.71352e-4,
            'Gd160': 5.90811e-4,
        }
    },
}


def build_model(particles=20000, batches=200, inactive=50):
    """Build the McSAFER NuScale-like SMR full-core benchmark model."""

    # =========================================================================
    # DIMENSIONS
    # =========================================================================
    pin_pitch = 1.2598       # cm
    assembly_pitch = 21.5036  # cm
    active_height = 200.0     # cm

    fuel_or = 0.4058    # fuel pellet outer radius
    clad_ir = 0.4140    # cladding inner radius
    clad_or = 0.4750    # cladding outer radius
    gt_ir = 0.5715      # guide tube inner radius
    gt_or = 0.6121      # guide tube outer radius

    core_barrel_ir = 93.98   # cm (ECP/RODARE value)
    core_barrel_or = 99.06   # cm
    rpv_ir = 122.68          # cm

    fuel_temp = 900.0   # K
    mod_temp = 600.0    # K, average coolant

    # =========================================================================
    # MATERIALS
    # =========================================================================

    # --- Fuel materials (one per enrichment + Gd variant) ---
    fuel_mats = {}
    for key, spec in FUEL_DATA.items():
        mat = openmc.Material(name=spec['name'])
        mat.set_density('sum')
        for nuclide, density in spec['nuclides'].items():
            mat.add_nuclide(nuclide, density)
        mat.temperature = fuel_temp
        fuel_mats[key] = mat

    # --- Helium gap ---
    he_gap = openmc.Material(name='Helium')
    he_gap.set_density('sum')
    he_gap.add_nuclide('He4', 2.4044e-4)
    he_gap.temperature = fuel_temp

    # --- Zircaloy-4 cladding (isotopic from RODARE) ---
    zr4 = openmc.Material(name='Zircaloy-4')
    zr4.set_density('sum')
    for nuc, nd in [
        ('Cr50', 3.2962e-6), ('Cr52', 6.3564e-5),
        ('Cr53', 7.2076e-6), ('Cr54', 1.7941e-6),
        ('Fe54', 8.6698e-6), ('Fe56', 1.361e-4),
        ('Fe57', 3.1431e-6), ('Fe58', 4.1829e-7),
        ('O16', 3.0744e-4),
        ('Sn112', 4.6735e-6), ('Sn114', 3.1799e-6),
        ('Sn115', 1.6381e-6), ('Sn116', 7.0055e-5),
        ('Sn117', 3.7003e-5), ('Sn118', 1.1669e-4),
        ('Sn119', 4.1387e-5), ('Sn120', 1.5697e-4),
        ('Sn122', 2.2308e-5), ('Sn124', 2.7897e-5),
        ('Zr90', 2.1828e-2), ('Zr91', 4.7601e-3),
        ('Zr92', 7.2759e-3), ('Zr94', 7.3734e-3),
        ('Zr96', 1.1879e-3),
    ]:
        zr4.add_nuclide(nuc, nd)
    zr4.temperature = mod_temp

    # --- Borated water coolant (1000 ppm boron) ---
    water = openmc.Material(name='Borated Water 1000ppm')
    water.set_density('sum')
    water.add_nuclide('H1', 5.02932e-2)
    water.add_nuclide('O16', 2.51573e-2)
    water.add_nuclide('B10', 8.33778e-6)
    water.add_nuclide('B11', 3.35608e-5)
    water.add_s_alpha_beta('c_H_in_H2O')
    water.temperature = mod_temp

    # --- SS304L core barrel (isotopic from RODARE) ---
    ss304l = openmc.Material(name='SS304L')
    ss304l.set_density('sum')
    for nuc, nd in [
        ('Cr50', 7.6778e-4), ('Cr52', 1.4806e-2),
        ('Cr53', 1.6789e-3), ('Cr54', 4.1791e-4),
        ('Fe54', 3.462e-3), ('Fe56', 5.4345e-2),
        ('Fe57', 1.2551e-3), ('Fe58', 1.6703e-4),
        ('Mn55', 1.7604e-3),
        ('Ni58', 5.6089e-3), ('Ni60', 2.1605e-3),
        ('Ni61', 9.3917e-5), ('Ni62', 2.9945e-4),
        ('Ni64', 7.6261e-5),
        ('Si28', 9.5281e-4), ('Si29', 4.8381e-5),
        ('Si30', 3.1893e-5),
    ]:
        ss304l.add_nuclide(nuc, nd)
    ss304l.temperature = mod_temp

    all_materials = list(fuel_mats.values()) + [he_gap, zr4, water, ss304l]
    materials = openmc.Materials(all_materials)
    materials.cross_sections = '/data/endfb-viii.0-hdf5/cross_sections.xml'

    # =========================================================================
    # PIN SURFACES
    # =========================================================================
    fuel_or_s = openmc.ZCylinder(r=fuel_or)
    clad_ir_s = openmc.ZCylinder(r=clad_ir)
    clad_or_s = openmc.ZCylinder(r=clad_or)
    gt_ir_s = openmc.ZCylinder(r=gt_ir)
    gt_or_s = openmc.ZCylinder(r=gt_or)

    # =========================================================================
    # PIN UNIVERSES
    # =========================================================================

    def make_fuel_pin(fuel_mat, name):
        """Create fuel pin: pellet -> gap -> clad -> water."""
        cells = [
            openmc.Cell(name=f'{name} pellet', fill=fuel_mat,
                        region=-fuel_or_s),
            openmc.Cell(name=f'{name} gap', fill=he_gap,
                        region=+fuel_or_s & -clad_ir_s),
            openmc.Cell(name=f'{name} clad', fill=zr4,
                        region=+clad_ir_s & -clad_or_s),
            openmc.Cell(name=f'{name} water', fill=water,
                        region=+clad_or_s),
        ]
        return openmc.Universe(name=name, cells=cells)

    # Fuel pin universes for each material type
    fuel_pins = {}
    for key, mat in fuel_mats.items():
        fuel_pins[key] = make_fuel_pin(mat, f'Pin {key}')

    # Guide tube universe (water-filled, ARO condition)
    gt_univ = openmc.Universe(name='Guide Tube', cells=[
        openmc.Cell(name='GT water', fill=water, region=-gt_ir_s),
        openmc.Cell(name='GT wall', fill=zr4,
                    region=+gt_ir_s & -gt_or_s),
        openmc.Cell(name='GT outer', fill=water, region=+gt_or_s),
    ])

    # Instrument tube (same dimensions as GT per benchmark)
    it_univ = openmc.Universe(name='Instrument Tube', cells=[
        openmc.Cell(name='IT water', fill=water, region=-gt_ir_s),
        openmc.Cell(name='IT wall', fill=zr4,
                    region=+gt_ir_s & -gt_or_s),
        openmc.Cell(name='IT outer', fill=water, region=+gt_or_s),
    ])

    # =========================================================================
    # ASSEMBLY UNIVERSES (one per type)
    # =========================================================================
    # Outer universe for lattice (inter-assembly water)
    lat_outer = openmc.Universe(name='Lattice outer',
                                cells=[openmc.Cell(fill=water)])

    assembly_types = ['A01', 'A02', 'B01', 'B02', 'C01', 'C02', 'C03']
    assembly_universes = {}

    for assy_type in assembly_types:
        lattice = openmc.RectLattice(name=f'Assembly {assy_type}')
        lattice.pitch = (pin_pitch, pin_pitch)
        lattice.lower_left = (-17 * pin_pitch / 2.0, -17 * pin_pitch / 2.0)

        universes = []
        for row in range(17):
            row_list = []
            for col in range(17):
                pos = (row, col)
                if pos == INSTRUMENT_TUBE_POSITION:
                    row_list.append(it_univ)
                elif pos in GUIDE_TUBE_POSITIONS:
                    row_list.append(gt_univ)
                elif assy_type == 'C02' and pos in GD_ROD_POSITIONS:
                    row_list.append(fuel_pins['C02_Gd'])
                else:
                    row_list.append(fuel_pins[assy_type])
            universes.append(row_list)
        lattice.universes = universes
        lattice.outer = lat_outer

        assy_cell = openmc.Cell(name=f'Assy {assy_type}', fill=lattice)
        assembly_universes[assy_type] = openmc.Universe(
            name=f'Assembly {assy_type}', cells=[assy_cell])

    # Water universe for empty core positions
    water_univ = openmc.Universe(name='Water',
                                 cells=[openmc.Cell(fill=water)])

    # =========================================================================
    # CORE LATTICE (7x7)
    # =========================================================================
    core_lattice = openmc.RectLattice(name='NuScale Core 7x7')
    core_lattice.pitch = (assembly_pitch, assembly_pitch)
    core_lattice.lower_left = (-7 * assembly_pitch / 2.0,
                                -7 * assembly_pitch / 2.0)

    core_univs = []
    for row in CORE_MAP:
        row_list = []
        for val in row:
            if val is None:
                row_list.append(water_univ)
            else:
                row_list.append(assembly_universes[val])
        core_univs.append(row_list)
    core_lattice.universes = core_univs
    core_lattice.outer = water_univ

    # =========================================================================
    # ROOT GEOMETRY
    # =========================================================================
    barrel_ir_s = openmc.ZCylinder(r=core_barrel_ir)
    barrel_or_s = openmc.ZCylinder(r=core_barrel_or)
    rpv_s = openmc.ZCylinder(r=rpv_ir, boundary_type='vacuum')

    z_min = openmc.ZPlane(z0=-active_height / 2.0, boundary_type='reflective')
    z_max = openmc.ZPlane(z0=+active_height / 2.0, boundary_type='reflective')

    axial = +z_min & -z_max

    core_cell = openmc.Cell(name='Core', fill=core_lattice,
                            region=-barrel_ir_s & axial)
    barrel_cell = openmc.Cell(name='Core barrel', fill=ss304l,
                              region=+barrel_ir_s & -barrel_or_s & axial)
    downcomer_cell = openmc.Cell(name='Downcomer', fill=water,
                                 region=+barrel_or_s & -rpv_s & axial)

    root = openmc.Universe(name='Root',
                           cells=[core_cell, barrel_cell, downcomer_cell])
    geometry = openmc.Geometry(root)

    # =========================================================================
    # SETTINGS
    # =========================================================================
    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive

    settings.temperature = {'method': 'interpolation', 'default': mod_temp}

    hw = 3.5 * assembly_pitch
    source = openmc.IndependentSource(
        space=openmc.stats.Box([-hw, -hw, -active_height / 2.0],
                               [hw, hw, active_height / 2.0]))
    source.constraints = {'fissionable': True}
    settings.source = [source]

    # =========================================================================
    # TALLIES
    # =========================================================================
    tallies = openmc.Tallies()
    mesh = openmc.RegularMesh(name='Assembly power mesh')
    mesh.dimension = (7, 7)
    mesh.lower_left = (-hw, -hw)
    mesh.upper_right = (hw, hw)
    tally = openmc.Tally(name='Assembly Powers')
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ['fission']
    tallies.append(tally)

    # =========================================================================
    # MODEL
    # =========================================================================
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings, tallies=tallies)

    # Print summary
    counts = {}
    for row in CORE_MAP:
        for v in row:
            if v is not None:
                counts[v] = counts.get(v, 0) + 1
    print("McSAFER NuScale-like SMR Benchmark")
    print("=" * 50)
    for t in assembly_types:
        enr = {'A01': '1.50', 'A02': '1.60', 'B01': '2.50', 'B02': '2.60',
               'C01': '4.05', 'C02': '4.55+Gd', 'C03': '2.60'}[t]
        print(f"  {t} ({enr:>7s}%): {counts.get(t, 0):2d} assemblies")
    print(f"  Total: 37 assemblies")
    print(f"  Boron: 1000 ppm, Fuel T: {fuel_temp} K, Mod T: {mod_temp} K")
    print("=" * 50)

    return model


def main():
    parser = argparse.ArgumentParser(
        description='McSAFER NuScale-like SMR Full-Core Benchmark')
    parser.add_argument('--particles', type=int, default=20000)
    parser.add_argument('--batches', type=int, default=200)
    parser.add_argument('--inactive', type=int, default=50)
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()

    model = build_model(particles=args.particles, batches=args.batches,
                        inactive=args.inactive)
    model.export_to_xml()
    print("\nExported to XML files")

    if args.run:
        openmc.run()


if __name__ == '__main__':
    main()
