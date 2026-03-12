#!/usr/bin/env python3
"""
Generate 2D slice plots of the VVER-1000 hexagonal fuel assembly geometry.

This script creates three plots:
  1. Full assembly XY cross-section showing hexagonal pin arrangement
  2. Zoomed view of the central region showing pin detail
  3. Zoomed view of fuel pin, guide tube, and central tube detail

The hexagonal lattice structure is the defining feature of VVER reactors
(as opposed to the square lattice used in Western PWR designs). The plots
clearly show the triangular arrangement of pins within each hexagonal ring.

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

    # Get material references by name for consistent color mapping
    materials = {m.name: m for m in model.materials}
    fuel_mat = materials['UO2 Fuel 3.0%']
    zirc_mat = materials['Zr Alloy Cladding']
    water_mat = materials['Borated Water']
    b4c_mat = materials['B4C Absorber']
    steel_mat = materials['Stainless Steel']

    assembly_pitch = 23.6  # cm, flat-to-flat

    # Color scheme: physically intuitive colors
    # Red for fuel, grey for cladding, blue for water, black for absorber
    mat_colors = {
        fuel_mat:  (200, 40, 40),       # dark red - UO2 fuel pellets
        zirc_mat:  (180, 180, 180),     # grey - Zr alloy cladding/tubes
        water_mat: (100, 150, 255),     # blue - borated water moderator
        b4c_mat:   (30, 30, 30),        # near-black - B4C absorber
        steel_mat: (120, 100, 80),      # dark brown - stainless steel
    }

    # --- Plot 1: Full assembly, material coloring ---
    # This shows the entire hexagonal assembly with all 331 lattice positions.
    # The hexagonal boundary and the triangular arrangement of pins within
    # each ring are clearly visible. The 312 fuel pins appear as red dots,
    # the 18 guide tubes appear as darker circles, and the central tube
    # is at the assembly center.
    full_plot = openmc.Plot(name='assembly_full_xy')
    full_plot.filename = 'plot_assembly_full'
    full_plot.origin = (0.0, 0.0, 0.0)
    full_plot.width = (assembly_pitch * 1.05, assembly_pitch * 1.05)
    full_plot.pixels = (1500, 1500)
    full_plot.color_by = 'material'
    full_plot.colors = mat_colors

    # --- Plot 2: Zoomed center showing central tube and inner rings ---
    # This zoomed view shows the central instrumentation tube (Zr alloy)
    # surrounded by the first few rings of fuel pins, and some of the
    # absorber guide tubes. The different tube geometries are clearly
    # distinguishable: fuel pins (small red), guide tubes (larger with
    # dark absorber center), and the central tube (medium, water-filled).
    zoom_center = openmc.Plot(name='assembly_zoom_center')
    zoom_center.filename = 'plot_assembly_zoom_center'
    zoom_center.origin = (0.0, 0.0, 0.0)
    zoom_center.width = (8.0, 8.0)  # show ~6 rings
    zoom_center.pixels = (1200, 1200)
    zoom_center.color_by = 'material'
    zoom_center.colors = mat_colors

    # --- Plot 3: Tight zoom showing individual pin detail ---
    # This view is zoomed in enough to clearly see the concentric ring
    # structure of individual pins: fuel pellet -> cladding -> water gap
    # for fuel pins, and absorber -> clad -> water -> guide tube for
    # guide tubes.
    zoom_pin = openmc.Plot(name='assembly_zoom_pin')
    zoom_pin.filename = 'plot_assembly_zoom_pin'
    zoom_pin.origin = (0.0, 3.5, 0.0)  # offset to show a guide tube area
    zoom_pin.width = (4.0, 4.0)
    zoom_pin.pixels = (1000, 1000)
    zoom_pin.color_by = 'material'
    zoom_pin.colors = mat_colors

    # Collect all plots and export
    plots = openmc.Plots([full_plot, zoom_center, zoom_pin])
    plots.export_to_xml()

    print("Plot XML exported. Run 'openmc --plot' to generate images.")
    print("Plots defined:")
    print(f"  1. {full_plot.filename} - Full hexagonal assembly "
          f"({full_plot.pixels[0]}x{full_plot.pixels[1]} px)")
    print(f"  2. {zoom_center.filename} - Zoomed center region "
          f"({zoom_center.pixels[0]}x{zoom_center.pixels[1]} px)")
    print(f"  3. {zoom_pin.filename} - Zoomed pin detail "
          f"({zoom_pin.pixels[0]}x{zoom_pin.pixels[1]} px)")


if __name__ == '__main__':
    make_plots()
