#!/usr/bin/env python3
"""
Generate 2D slice plots of the BEAVRS single fuel assembly geometry.

This script creates four plots:
  1. Full assembly XY cross-section colored by material
  2. Full assembly XY cross-section colored by cell
  3. Zoomed view showing fuel pin and guide tube detail
  4. Quarter-assembly zoom colored by material

The BEAVRS assembly has slightly different dimensions from the VERA benchmark
(smaller fuel pellet radius, different pitch), but the overall layout is the
same standard Westinghouse 17x17 pattern.

Prerequisites:
    Run model.py first to generate the XML input files.

Usage:
    python visualize.py
"""

import openmc
from model import build_model


def make_plots():
    """Create and export geometry visualization plots."""

    # Build model to get material references for custom coloring
    model = build_model()

    # Get material references by name for consistent coloring
    materials = {m.name: m for m in model.materials}
    fuel_mat = materials['UO2 Fuel 3.1%']
    gap_mat = materials['Helium Gap']
    clad_mat = materials['Zircaloy-4']
    mod_mat = materials['Borated Water']

    # BEAVRS assembly pitch
    assembly_pitch = 21.50364  # cm

    # Color scheme: consistent across all plots.
    # Red for fuel (hot!), yellow for gap, grey for structural metal,
    # blue for water.
    mat_colors = {
        fuel_mat: (180, 30, 30),       # dark red - UO2 fuel
        gap_mat:  (255, 255, 100),     # yellow - helium gap
        clad_mat: (160, 160, 160),     # grey - Zircaloy-4
        mod_mat:  (100, 150, 255),     # blue - borated water
    }

    # --- Plot 1: Full assembly, material coloring ---
    # Shows the entire 17x17 array at the BEAVRS assembly pitch.
    # Fuel pins appear as red dots with grey cladding rings. Guide tubes
    # and instrument tube appear as larger blue circles with grey walls.
    mat_plot = openmc.Plot(name='assembly_materials_xy')
    mat_plot.filename = 'plot_assembly_materials'
    mat_plot.origin = (0.0, 0.0, 0.0)
    mat_plot.width = (assembly_pitch, assembly_pitch)
    mat_plot.pixels = (1200, 1200)  # high resolution for 17x17 detail
    mat_plot.color_by = 'material'
    mat_plot.colors = mat_colors

    # --- Plot 2: Full assembly, cell coloring ---
    # Each cell type gets a different color. Useful for verifying that the
    # correct universe (fuel pin / guide tube / instrument tube) is placed
    # at each lattice position.
    cell_plot = openmc.Plot(name='assembly_cells_xy')
    cell_plot.filename = 'plot_assembly_cells'
    cell_plot.origin = (0.0, 0.0, 0.0)
    cell_plot.width = (assembly_pitch, assembly_pitch)
    cell_plot.pixels = (1200, 1200)
    cell_plot.color_by = 'cell'

    # --- Plot 3: Zoomed view of fuel pin and guide tube region ---
    # Zoom into the area around a guide tube and adjacent fuel pins.
    # Guide tube at position (5,5) in 0-indexed lattice coordinates:
    #   x offset from center = (5 - 8) * pitch = -3 * 1.25984 = -3.77952
    #   y offset from center = (8 - 5) * pitch = +3 * 1.25984 = +3.77952
    # (y is flipped because row 0 is top = highest y)
    pin_pitch = 1.25984
    zoom_x = (5 - 8) * pin_pitch  # near guide tube at col 5
    zoom_y = (8 - 5) * pin_pitch  # near guide tube at row 5
    zoom_plot = openmc.Plot(name='zoom_fuel_and_gt')
    zoom_plot.filename = 'plot_zoom_detail'
    zoom_plot.origin = (zoom_x, zoom_y, 0.0)
    zoom_plot.width = (4.0, 4.0)  # ~3 pin pitches across
    zoom_plot.pixels = (800, 800)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = mat_colors

    # --- Plot 4: Quarter assembly zoom ---
    # Shows the top-right quadrant of the assembly. Good for seeing the
    # guide tube pattern within one quadrant.
    quarter_plot = openmc.Plot(name='quarter_assembly')
    quarter_plot.filename = 'plot_quarter_assembly'
    quarter_plot.origin = (5.0, 5.0, 0.0)
    quarter_plot.width = (11.0, 11.0)
    quarter_plot.pixels = (900, 900)
    quarter_plot.color_by = 'material'
    quarter_plot.colors = mat_colors

    # Collect and export all plots
    plots = openmc.Plots([mat_plot, cell_plot, zoom_plot, quarter_plot])
    plots.export_to_xml()
    print("Exported plots.xml with 4 plot definitions.")

    # Run OpenMC in plotting mode to generate PNG files
    openmc.plot_geometry()
    print("Generated plot images:")
    print("  - plot_assembly_materials.png  (full assembly, material colors)")
    print("  - plot_assembly_cells.png      (full assembly, cell colors)")
    print("  - plot_zoom_detail.png         (zoomed fuel pin + guide tube)")
    print("  - plot_quarter_assembly.png    (quarter assembly view)")


if __name__ == '__main__':
    make_plots()
