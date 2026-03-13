#!/usr/bin/env python3
"""
Generate 2D slice plots of the NuScale-like SMR full-core geometry.

This script creates visualization plots for the 37-assembly full core:
  1. Full core XY cross-section colored by material
  2. Full core XY cross-section colored by cell
  3. XZ axial slice through the core center
  4. Zoomed view showing a single assembly with pin-level detail

The core barrel and downcomer water are visible as concentric rings
surrounding the fuel assembly lattice.

Prerequisites:
    Run model.py first to generate the XML input files, or this script
    will build the model directly.

Usage:
    python visualize.py
"""

import openmc
from model import build_model


def make_plots():
    """Create and export geometry visualization plots."""

    # Build the full-core model
    model = build_model()

    # Get material references for custom coloring
    materials = {m.name: m for m in model.materials}
    fuel_16 = materials['UO2 Fuel 1.6%']
    fuel_24 = materials['UO2 Fuel 2.4%']
    fuel_31 = materials['UO2 Fuel 3.1%']
    gap_mat = materials['Helium Gap']
    clad_mat = materials['Zircaloy-4']
    mod_mat = materials['Unborated Water']
    ss_mat = materials['SS304']

    # Color scheme:
    #   Different reds for different enrichment zones (darker = higher enrichment)
    #   Yellow = helium gap
    #   Grey = Zircaloy-4 cladding and tube walls
    #   Blue = unborated water moderator
    #   Dark grey = SS304 core barrel
    mat_colors = {
        fuel_16: (255, 120, 120),      # light red - 1.6% fuel (center)
        fuel_24: (220, 60, 60),        # medium red - 2.4% fuel (middle)
        fuel_31: (160, 20, 20),        # dark red - 3.1% fuel (outer)
        gap_mat:  (255, 255, 100),     # yellow - helium gap
        clad_mat: (180, 180, 180),     # light grey - Zircaloy-4
        mod_mat:  (100, 150, 255),     # blue - unborated water
        ss_mat:   (100, 100, 100),     # dark grey - SS304
    }

    # Core dimensions for plot sizing
    assembly_pitch = 21.50   # cm
    core_barrel_or = 105.0   # cm
    rpv_ir = 115.0           # cm
    active_height = 200.0    # cm

    # --- Plot 1: Full core XY, material coloring ---
    # Shows all 37 assemblies with the three enrichment zones visible
    # as different shades of red. The core barrel (dark grey) and
    # downcomer water (blue) surround the assembly lattice.
    core_width = 2.0 * rpv_ir + 10.0  # slightly larger than RPV
    core_mat_plot = openmc.Plot(name='core_materials_xy')
    core_mat_plot.filename = 'plot_core_materials'
    core_mat_plot.origin = (0.0, 0.0, 0.0)
    core_mat_plot.width = (core_width, core_width)
    core_mat_plot.pixels = (2000, 2000)
    core_mat_plot.color_by = 'material'
    core_mat_plot.colors = mat_colors

    # --- Plot 2: Full core XY, cell coloring ---
    # Each cell gets a distinct auto-assigned color, useful for verifying
    # that the correct assembly type is placed at each core position.
    core_cell_plot = openmc.Plot(name='core_cells_xy')
    core_cell_plot.filename = 'plot_core_cells'
    core_cell_plot.origin = (0.0, 0.0, 0.0)
    core_cell_plot.width = (core_width, core_width)
    core_cell_plot.pixels = (2000, 2000)
    core_cell_plot.color_by = 'cell'

    # --- Plot 3: XZ axial slice through the core center ---
    # Shows the active fuel height (200 cm) with reflective axial BCs.
    # The radial extent shows the core barrel and downcomer.
    xz_plot = openmc.Plot(name='core_xz_slice')
    xz_plot.filename = 'plot_core_xz'
    xz_plot.basis = 'xz'
    xz_plot.origin = (0.0, 0.0, 0.0)
    xz_plot.width = (core_width, active_height + 20.0)
    xz_plot.pixels = (2000, 1800)
    xz_plot.color_by = 'material'
    xz_plot.colors = mat_colors

    # --- Plot 4: Zoomed view of one assembly ---
    # Center on the (0,0) assembly position (which is at the core center)
    # to show individual fuel pins, guide tubes, and the instrument tube.
    zoom_plot = openmc.Plot(name='assembly_zoom_xy')
    zoom_plot.filename = 'plot_assembly_zoom'
    zoom_plot.origin = (0.0, 0.0, 0.0)
    zoom_plot.width = (assembly_pitch, assembly_pitch)
    zoom_plot.pixels = (1200, 1200)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = mat_colors

    # Collect and export all plots
    plots = openmc.Plots([core_mat_plot, core_cell_plot, xz_plot, zoom_plot])
    plots.export_to_xml()
    print("Exported plots.xml with 4 plot definitions.")

    # Run OpenMC in plotting mode
    openmc.plot_geometry()
    print("Generated plot images:")
    print("  - plot_core_materials.png   (full core XY, material colors)")
    print("  - plot_core_cells.png       (full core XY, cell colors)")
    print("  - plot_core_xz.png          (XZ axial slice)")
    print("  - plot_assembly_zoom.png    (single assembly detail)")


if __name__ == '__main__':
    make_plots()
