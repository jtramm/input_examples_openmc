#!/usr/bin/env python3
"""
FNG Copper Benchmark -- Geometry Visualization
================================================

Produces an XZ cross-section plot of the FNG copper slab assembly and D-T
source position. The plot slices through y=0 (the midplane), showing the
beam axis (x) horizontally and the vertical axis (z) vertically.

The seven copper plates appear as coloured rectangles and the source
position is indicated by the void region to the left of the assembly.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    fng_copper_geometry.png
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
    # Create an XZ cross-section plot (y = 0 midplane)
    # =========================================================================
    # This slice shows the beam axis (x) horizontally and the vertical axis
    # (z) vertically. The rectangular copper assembly appears as a filled
    # rectangle. The source is in the void region to the left.
    #
    # Assembly: x = 0 to ~69.9 cm, z = -30 to +30 cm
    # Source:   x = -5.3 cm
    # We centre the view to show both source and full assembly.
    plot = openmc.Plot()
    plot.filename = "fng_copper_geometry"
    plot.basis = "xz"                             # plot in the x-z plane
    plot.origin = (32.0, 0.0, 0.0)               # centre of view (midpoint of assembly)
    plot.width = (90.0, 80.0)                     # field of view in cm (x, z)
    plot.pixels = (900, 800)                      # image resolution
    plot.color_by = "material"                    # colour cells by material assignment

    # Custom material colours for clarity:
    # Copper: distinctive orange-brown colour reminiscent of real copper
    plot.colors = {
        materials[0]: (184, 115, 51),             # copper-brown for OFHC copper
    }

    # =========================================================================
    # Export and generate the plot
    # =========================================================================
    plots = openmc.Plots([plot])
    plots.export_to_xml()

    # Generate the plot image using OpenMC's plotter
    openmc.plot_geometry()

    print("Geometry plot saved to: fng_copper_geometry.png")
    print("  - View: XZ cross-section at y = 0 (midplane)")
    print("  - Copper-brown rectangles: seven OFHC copper plates")
    print("  - Source location: x = -5.3 cm (in void region, left of assembly)")


if __name__ == "__main__":
    main()
