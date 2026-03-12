#!/usr/bin/env python3
"""
FNG HCPB Tritium Breeder Module Mock-up -- Geometry Visualization
===================================================================

Produces cross-section plots of the FNG HCPB assembly showing the layered
structure of beryllium (neutron multiplier), Li2CO3 breeder layers, and
stainless steel walls.

Two views are generated:
  1. XY cross-section (z=0 midplane): shows the beam axis (y) horizontally
     and the lateral extent (x) vertically. This is the primary view showing
     the layer structure.
  2. XZ cross-section (y=10 cm): shows a slice through the central beryllium
     zone between the two breeder sections.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    fng_hcpb_geometry_xy.png   (XY cross-section, beam axis view)
    fng_hcpb_geometry_xz.png   (XZ cross-section, transverse view)
"""

import openmc


def main():
    # =========================================================================
    # Load the model geometry and materials from XML files
    # =========================================================================
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    # =========================================================================
    # Material colour assignments
    # =========================================================================
    # Colours chosen to visually distinguish the different materials:
    #   - Beryllium: grey (metallic appearance)
    #   - Li2CO3 front: green (breeder material, front sections)
    #   - Li2CO3 rear: light green (breeder material, rear cassette)
    #   - AISI-303 SS: blue (structural steel, main box)
    #   - AISI-316 SS: dark blue (structural steel, rear cassette)
    mat_colors = {}
    for mat in materials:
        if "Beryllium" in mat.name:
            mat_colors[mat] = (180, 180, 180)       # grey
        elif "front" in mat.name:
            mat_colors[mat] = (50, 180, 50)          # green
        elif "rear" in mat.name:
            mat_colors[mat] = (120, 220, 120)        # light green
        elif "303" in mat.name:
            mat_colors[mat] = (70, 70, 200)          # blue
        elif "316" in mat.name:
            mat_colors[mat] = (40, 40, 140)          # dark blue

    # =========================================================================
    # Plot 1: XY cross-section at z=0 (beam axis view)
    # =========================================================================
    # This is the most informative view, showing the full layer structure
    # along the beam axis (y): SS walls, beryllium zones, breeder double-layers.
    # The source position is in the void region to the left (negative y).
    #
    # Assembly: y = 0 to 43.8 cm, x = -15.5 to +15.5 cm
    plot_xy = openmc.Plot()
    plot_xy.filename = "fng_hcpb_geometry_xy"
    plot_xy.basis = "xy"                              # plot in x-y plane
    plot_xy.origin = (0.0, 20.0, 0.0)                # centre of view
    plot_xy.width = (45.0, 55.0)                      # field of view (x, y) in cm
    plot_xy.pixels = (900, 1100)                      # image resolution
    plot_xy.color_by = "material"
    plot_xy.colors = mat_colors

    # =========================================================================
    # Plot 2: XZ cross-section at y=10 cm (transverse view through Be zone 2)
    # =========================================================================
    # Shows the rectangular cross-section of the assembly at a depth of 10 cm,
    # which is in the middle of the central beryllium zone (between the two
    # breeder sections). This view confirms the lateral extent of the assembly.
    plot_xz = openmc.Plot()
    plot_xz.filename = "fng_hcpb_geometry_xz"
    plot_xz.basis = "xz"                              # plot in x-z plane
    plot_xz.origin = (0.0, 10.0, 0.0)                 # slice at y=10 cm (Be zone 2)
    plot_xz.width = (45.0, 45.0)                      # field of view (x, z) in cm
    plot_xz.pixels = (900, 900)                        # image resolution
    plot_xz.color_by = "material"
    plot_xz.colors = mat_colors

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_xy, plot_xz])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  1. fng_hcpb_geometry_xy.png -- XY cross-section (beam axis view)")
    print("     Shows layer structure: SS walls, Be zones, Li2CO3 breeder layers")
    print("  2. fng_hcpb_geometry_xz.png -- XZ cross-section at y=10 cm")
    print("     Shows assembly cross-section through central beryllium zone")
    print()
    print("Material colour key:")
    print("  Grey:        Beryllium (neutron multiplier)")
    print("  Green:       Li2CO3 breeder (front sections, 1.123 g/cm3)")
    print("  Light green: Li2CO3 breeder (rear cassette, 0.9413 g/cm3)")
    print("  Blue:        AISI-303 stainless steel (main box)")
    print("  Dark blue:   AISI-316 stainless steel (rear cassette)")


if __name__ == "__main__":
    main()
