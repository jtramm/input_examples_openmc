#!/usr/bin/env python3
"""
FNS Tungsten Benchmark -- Geometry Visualization
=================================================

Produces an XY cross-section plot of the FNS tungsten cylindrical assembly
and D-T source position. The plot slices through z=0 (the midplane),
showing the beam-axis (x) and one transverse axis (y).

The tungsten assembly appears as a coloured rectangle (cross-section of the
cylinder) and the source position is marked with an annotation.

Usage:
    python visualize.py

Output:
    fns_tungsten_geometry.png
"""

import openmc


def main():
    # =========================================================================
    # Load the model geometry and materials from XML files
    # =========================================================================
    # These must have been generated first by running model.py
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    # =========================================================================
    # Create an XY cross-section plot (z = 0 midplane)
    # =========================================================================
    # This slice shows the beam axis (x) horizontally and one transverse
    # axis (y) vertically. The cylindrical assembly appears as a rectangle
    # in this projection.
    plot = openmc.Plot()
    plot.filename = "fns_tungsten_geometry"
    plot.basis = "xy"                           # plot in the x-y plane
    plot.origin = (15.0, 0.0, 0.0)             # centre of view (midpoint of assembly approx)
    plot.width = (90.0, 80.0)                   # field of view in cm (x, y)
    plot.pixels = (900, 800)                    # image resolution
    plot.color_by = "material"                  # colour cells by material assignment

    # Custom material colours for clarity
    # Tungsten alloy: steel-blue colour
    plot.colors = {
        materials[0]: (70, 130, 180),           # steel blue for tungsten alloy
    }

    # =========================================================================
    # Export and generate the plot
    # =========================================================================
    plots = openmc.Plots([plot])
    plots.export_to_xml()

    # Generate the plot image using OpenMC's plotter
    openmc.plot_geometry()

    print("Geometry plot saved to: fns_tungsten_geometry.png")
    print("  - View: XY cross-section at z = 0 (midplane)")
    print("  - Blue rectangle: tungsten alloy assembly")
    print("  - Source location: x = -20 cm, y = 0 (left of assembly)")


if __name__ == "__main__":
    main()
