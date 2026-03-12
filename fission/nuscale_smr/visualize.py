#!/usr/bin/env python3
"""
Generate 2D slice plots of the NuScale-like SMR fuel assembly geometry.

This script creates visualization plots for the 17x17 fuel assembly model:
  1. Full assembly XY cross-section colored by material
  2. Full assembly XY cross-section colored by cell
  3. Zoomed view showing fuel pin and guide tube detail

The zoomed view highlights the difference between fuel pins (with UO2 pellet,
gap, and cladding) and the larger-diameter guide tubes (water-filled when
control rods are withdrawn).

Prerequisites:
    Run model.py first to generate the XML input files, or this script
    will build the model directly.

Usage:
    python visualize.py
"""

import openmc
from model import build_assembly


def make_plots():
    """Create and export geometry visualization plots."""

    # Build the assembly model to get material references for custom coloring
    model = build_assembly()

    # Get material references by name for consistent color assignment
    materials = {m.name: m for m in model.materials}
    fuel_mat = materials['UO2 Fuel 3.1%']
    gap_mat = materials['Helium Gap']
    clad_mat = materials['Zircaloy-4']
    mod_mat = materials['Unborated Water']

    # Assembly pitch for setting plot widths
    assembly_pitch = 21.5036  # cm

    # Color scheme: consistent colors across all material plots
    #   Red   = UO2 fuel (hot, fissile material)
    #   Yellow = helium gap (thin, hard to see at assembly scale)
    #   Grey  = Zircaloy-4 cladding and tube walls
    #   Blue  = unborated water moderator
    mat_colors = {
        fuel_mat: (180, 30, 30),       # dark red - UO2 fuel
        gap_mat:  (255, 255, 100),     # yellow - helium gap
        clad_mat: (160, 160, 160),     # grey - Zircaloy-4
        mod_mat:  (100, 150, 255),     # blue - unborated water
    }

    # --- Plot 1: Full assembly, material coloring ---
    # Shows the entire 17x17 array with material-based coloring.
    # Fuel pins appear as red dots (UO2) surrounded by grey cladding rings.
    # Guide tubes and the instrument tube appear as larger blue circles
    # (water-filled) with grey tube walls. The inter-assembly water gap
    # is visible as a thin blue border around the lattice.
    mat_plot = openmc.Plot(name='assembly_materials_xy')
    mat_plot.filename = 'plot_assembly_materials'
    mat_plot.origin = (0.0, 0.0, 0.0)
    mat_plot.width = (assembly_pitch, assembly_pitch)
    mat_plot.pixels = (1200, 1200)  # high resolution for 17x17 detail
    mat_plot.color_by = 'material'
    mat_plot.colors = mat_colors

    # --- Plot 2: Full assembly, cell coloring ---
    # Each cell gets a distinct color, useful for verifying that the
    # correct universe type (fuel pin, guide tube, instrument tube) is
    # placed at each lattice position.
    cell_plot = openmc.Plot(name='assembly_cells_xy')
    cell_plot.filename = 'plot_assembly_cells'
    cell_plot.origin = (0.0, 0.0, 0.0)
    cell_plot.width = (assembly_pitch, assembly_pitch)
    cell_plot.pixels = (1200, 1200)
    cell_plot.color_by = 'cell'

    # --- Plot 3: Zoomed view of fuel pin and guide tube region ---
    # Zoom into an area near a guide tube to show the structural detail:
    # the concentric rings of the fuel pin (pellet, gap, clad) and the
    # larger guide tube (inner water, tube wall, outer water).
    #
    # We center on the guide tube at lattice position (5,5):
    #   x = (5 - 8) * 1.2598 = -3.7794 cm (offset from center)
    #   y = (8 - 5) * 1.2598 = +3.7794 cm
    zoom_x = -3.78  # near guide tube at column 5
    zoom_y = 3.78   # near guide tube at row 5
    zoom_plot = openmc.Plot(name='zoom_fuel_and_gt')
    zoom_plot.filename = 'plot_zoom_detail'
    zoom_plot.origin = (zoom_x, zoom_y, 0.0)
    zoom_plot.width = (5.0, 5.0)  # ~4 pin pitches across
    zoom_plot.pixels = (800, 800)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = mat_colors

    # Collect and export all plots
    plots = openmc.Plots([mat_plot, cell_plot, zoom_plot])
    plots.export_to_xml()
    print("Exported plots.xml with 3 plot definitions.")

    # Run OpenMC in plotting mode to generate PNG image files
    openmc.plot_geometry()
    print("Generated plot images:")
    print("  - plot_assembly_materials.png  (full assembly, material colors)")
    print("  - plot_assembly_cells.png      (full assembly, cell colors)")
    print("  - plot_zoom_detail.png         (zoomed fuel pin + guide tube)")


if __name__ == '__main__':
    make_plots()
