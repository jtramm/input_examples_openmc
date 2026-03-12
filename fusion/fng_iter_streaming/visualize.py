#!/usr/bin/env python3
"""
FNG-ITER Streaming Experiment -- Geometry Visualization
========================================================

Produces two cross-section plots of the FNG-ITER streaming assembly:

  1. Full XZ view showing the complete layer structure, channel, cavity,
     and coil simulator at the assembly midplane (y=0).

  2. Zoomed XZ view focusing on the channel entrance, streaming duct,
     and detector cavity region.

The beam/channel axis is along z (horizontal in plots), with x vertical.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    fng_streaming_full.png    -- full assembly cross-section
    fng_streaming_zoomed.png  -- zoomed channel/cavity region
"""

import openmc


def main():
    # =========================================================================
    # Load model geometry and materials from XML files
    # =========================================================================
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    # Build a material-colour dictionary for clarity:
    # SS316:        steel-grey
    # Copper:       copper-brown
    # Perspex:      light blue (water-equivalent)
    # Polyethylene: light green
    # Air:          white/light grey
    mat_colors = {}
    for mat in materials:
        if "SS316" in mat.name:
            mat_colors[mat] = (140, 140, 160)      # steel grey
        elif "Copper" in mat.name:
            mat_colors[mat] = (184, 115, 51)        # copper brown
        elif "Perspex" in mat.name:
            mat_colors[mat] = (135, 206, 235)       # light blue
        elif "Polyethylene" in mat.name:
            mat_colors[mat] = (144, 238, 144)       # light green
        elif "Air" in mat.name:
            mat_colors[mat] = (240, 240, 240)       # near-white

    # =========================================================================
    # Plot 1: Full assembly cross-section (XZ at y=0)
    # =========================================================================
    # This plot shows the entire assembly from source region to rear shield.
    # z-axis (beam) is horizontal, x-axis is vertical.
    plot_full = openmc.Plot()
    plot_full.filename = "fng_streaming_full"
    plot_full.basis = "xz"                          # x vertical, z horizontal
    plot_full.origin = (0.0, 0.0, 47.0)            # centre of view
    plot_full.width = (110.0, 110.0)                # field of view in cm
    plot_full.pixels = (1100, 1100)                 # resolution
    plot_full.color_by = "material"
    plot_full.colors = mat_colors

    # =========================================================================
    # Plot 2: Zoomed view of channel and cavity region
    # =========================================================================
    # Focus on z = -2 to 48 cm (channel + cavity + first plates behind)
    # and x = -5 to 5 cm (channel and cavity cross-section)
    plot_zoom = openmc.Plot()
    plot_zoom.filename = "fng_streaming_zoomed"
    plot_zoom.basis = "xz"
    plot_zoom.origin = (0.0, 0.0, 23.0)            # centred on channel midpoint
    plot_zoom.width = (12.0, 52.0)                  # narrow x view, full channel z
    plot_zoom.pixels = (600, 2600)                  # high z-resolution
    plot_zoom.color_by = "material"
    plot_zoom.colors = mat_colors

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_full, plot_zoom])
    plots.export_to_xml()
    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  fng_streaming_full.png   -- full assembly XZ cross-section at y=0")
    print("  fng_streaming_zoomed.png -- zoomed channel/cavity region")
    print()
    print("Colour legend:")
    print("  Copper-brown:  Copper (first wall, coil simulator)")
    print("  Steel-grey:    SS316 (channel wall, plates, cavity walls)")
    print("  Light blue:    Perspex (PMMA moderator blocks)")
    print("  Light green:   Polyethylene (rear shield)")
    print("  Near-white:    Air (channel bore, cavity, gaps)")


if __name__ == "__main__":
    main()
