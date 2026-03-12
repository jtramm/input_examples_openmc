#!/usr/bin/env python3
"""
FNG-ITER Bulk Shield Mock-up -- Geometry Visualization
========================================================

Produces an XZ cross-section plot of the FNG-ITER bulk shield mock-up,
showing the multi-layer slab geometry with distinct colours for each
material (copper, SS316, Perspex).

The plot slices through y=0 (the midplane), showing the beam axis (x)
horizontally and the vertical axis (z) vertically. The 17 layers of the
mock-up are visible as alternating coloured bands.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    fng_iter_shield_geometry.png
"""

import openmc


def main():
    # =========================================================================
    # Load the model geometry and materials from XML files
    # =========================================================================
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    # Identify materials by name for colour assignment
    mat_by_name = {m.name: m for m in materials}
    copper = mat_by_name["Copper"]
    ss316 = mat_by_name["SS316"]
    perspex = mat_by_name["Perspex (PMMA)"]

    # =========================================================================
    # Create an XZ cross-section plot (y = 0 midplane)
    # =========================================================================
    # This slice shows the beam axis (x) horizontally and the vertical axis
    # (z) vertically. The layered mock-up appears as alternating coloured
    # bands. The source position is in the void region to the left.
    #
    # Assembly: x = 0 to ~94.3 cm, z = -30 to +30 cm
    # Source:   x = -5.3 cm
    plot = openmc.Plot()
    plot.filename = "fng_iter_shield_geometry"
    plot.basis = "xz"                               # plot in the x-z plane
    plot.origin = (45.0, 0.0, 0.0)                 # centre of view
    plot.width = (115.0, 80.0)                       # field of view (x, z)
    plot.pixels = (1150, 800)                        # image resolution
    plot.color_by = "material"                       # colour by material type

    # Custom material colours chosen for visual clarity:
    #   Copper:  orange-brown (like real copper metal)
    #   SS316:   steel-blue/grey (like real stainless steel)
    #   Perspex: pale yellow/cream (like real clear PMMA)
    plot.colors = {
        copper:  (184, 115, 51),    # copper-brown
        ss316:   (140, 160, 180),   # steel-blue-grey
        perspex: (255, 240, 180),   # pale cream/yellow
    }

    # =========================================================================
    # Export and generate the plot
    # =========================================================================
    plots = openmc.Plots([plot])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plot saved to: fng_iter_shield_geometry.png")
    print("  - View: XZ cross-section at y = 0 (midplane)")
    print("  - Copper-brown bands: Cu layers (first wall, TF coil sections)")
    print("  - Steel-grey bands: SS316 layers (structure, VV, TF coil)")
    print("  - Pale yellow bands: Perspex layers (moderator)")
    print("  - Source location: x = -5.3 cm (in void region, left of assembly)")


if __name__ == "__main__":
    main()
