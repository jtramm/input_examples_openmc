#!/usr/bin/env python3
"""
Generate 2D slice plots of the VERA Problem 1A pin cell geometry.

This script creates three plots:
  1. Material-colored XY cross-section at z=0
  2. Cell-colored XY cross-section at z=0
  3. A zoomed-in material plot of the fuel/gap/clad region

Prerequisites:
    Run model.py first to generate the XML input files.

Usage:
    python visualize.py
"""

import openmc
from model import build_model


def make_plots():
    """Create and export geometry visualization plots."""

    # Build model to get material/cell references
    model = build_model()

    # Get material references from the model
    materials = {m.name: m for m in model.materials}
    fuel_mat = materials['UO2 Fuel 3.1%']
    gap_mat = materials['Helium Gap']
    clad_mat = materials['Zircaloy-4']
    mod_mat = materials['Borated Water']

    # --- Plot 1: Material coloring ---
    mat_plot = openmc.Plot(name='materials_xy')
    mat_plot.filename = 'plot_materials'
    mat_plot.origin = (0.0, 0.0, 0.0)
    mat_plot.width = (1.26, 1.26)
    mat_plot.pixels = (600, 600)
    mat_plot.color_by = 'material'
    mat_plot.colors = {
        fuel_mat: (180, 30, 30),       # dark red
        gap_mat:  (255, 255, 100),     # yellow
        clad_mat: (160, 160, 160),     # grey
        mod_mat:  (100, 150, 255),     # blue
    }

    # --- Plot 2: Cell coloring ---
    cell_plot = openmc.Plot(name='cells_xy')
    cell_plot.filename = 'plot_cells'
    cell_plot.origin = (0.0, 0.0, 0.0)
    cell_plot.width = (1.26, 1.26)
    cell_plot.pixels = (600, 600)
    cell_plot.color_by = 'cell'

    # --- Plot 3: Zoomed-in material plot ---
    zoom_plot = openmc.Plot(name='zoom_fuel_region')
    zoom_plot.filename = 'plot_zoom'
    zoom_plot.origin = (0.0, 0.0, 0.0)
    zoom_plot.width = (1.1, 1.1)
    zoom_plot.pixels = (800, 800)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = {
        fuel_mat: (180, 30, 30),
        gap_mat:  (255, 255, 100),
        clad_mat: (160, 160, 160),
        mod_mat:  (100, 150, 255),
    }

    # Collect and export
    plots = openmc.Plots([mat_plot, cell_plot, zoom_plot])
    plots.export_to_xml()
    print("Exported plots.xml")

    # Run OpenMC in plotting mode to generate PNG files
    openmc.plot_geometry()
    print("Generated plot images.")


if __name__ == '__main__':
    make_plots()
