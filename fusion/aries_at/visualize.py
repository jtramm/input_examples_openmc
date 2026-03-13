#!/usr/bin/env python3
"""
ARIES-AT Advanced Tokamak -- Geometry Visualization
======================================================

Produces cross-section plots of the ARIES-AT toroidal sector geometry showing
the radial build from plasma through the TF coil.

Two views are generated:
  1. RZ cross-section (XZ plane at y=0): shows the poloidal cross-section
     of the tokamak, with the radial build visible as concentric rings
     around the plasma. This is the most informative view for a toroidal
     geometry.
  2. XY midplane (z=0): shows the toroidal sector from above, with the
     22.5-degree wedge shape visible.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    aries_at_geometry_rz.png   (RZ poloidal cross-section)
    aries_at_geometry_xy.png   (XY toroidal midplane)
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
    # Colours chosen to match conventional fusion engineering colour schemes:
    #   - SiC: tan/brown (ceramic)
    #   - Pb-17Li: yellow (liquid metal)
    #   - SiC/LiPb blanket: orange (blended material)
    #   - WC/H2O shield: dark grey (dense shielding)
    #   - F82H steel: blue (structural steel)
    #   - TF coil: green (superconductor)
    mat_colors = {}
    for mat in materials:
        if mat.name == "SiC composite":
            mat_colors[mat] = (180, 140, 80)        # tan (ceramic)
        elif "Pb-17Li" in mat.name:
            mat_colors[mat] = (220, 200, 50)         # yellow (liquid metal)
        elif "blanket" in mat.name.lower():
            mat_colors[mat] = (230, 150, 50)         # orange (blanket)
        elif "Back wall" in mat.name:
            mat_colors[mat] = (200, 130, 40)         # darker orange
        elif "WC" in mat.name or "shield" in mat.name.lower():
            mat_colors[mat] = (100, 100, 100)        # dark grey (shield)
        elif "F82H" in mat.name:
            mat_colors[mat] = (70, 70, 200)          # blue (steel)
        elif "TF coil" in mat.name:
            mat_colors[mat] = (50, 180, 80)          # green (coil)

    # =========================================================================
    # Plot 1: RZ cross-section (poloidal view, XZ plane at y=0)
    # =========================================================================
    # This shows the tokamak in the R-Z plane, where R is the major radius
    # direction (x-axis) and Z is the vertical direction. The concentric
    # torus layers appear as nested annular rings around the plasma.
    #
    # View centred on the outboard midplane (R = R0 = 520 cm, Z = 0).
    # Width covers the full radial extent from inboard to outboard.
    plot_rz = openmc.Plot()
    plot_rz.filename = "aries_at_geometry_rz"
    plot_rz.basis = "xz"                            # R-Z plane (y=0)
    plot_rz.origin = (520.0, 0.0, 0.0)              # centred on plasma axis
    plot_rz.width = (700.0, 700.0)                   # covers full torus cross-section
    plot_rz.pixels = (1200, 1200)                    # high resolution
    plot_rz.color_by = "material"
    plot_rz.colors = mat_colors

    # =========================================================================
    # Plot 2: XY midplane (toroidal view, z=0)
    # =========================================================================
    # This shows the tokamak from above at the midplane. The 22.5-degree
    # sector is visible as a wedge shape with the reflective boundary
    # planes on either side. The radial layers appear as concentric arcs.
    plot_xy = openmc.Plot()
    plot_xy.filename = "aries_at_geometry_xy"
    plot_xy.basis = "xy"                            # toroidal midplane
    plot_xy.origin = (400.0, 150.0, 0.0)            # centred to show the sector
    plot_xy.width = (700.0, 500.0)                   # covers sector width
    plot_xy.pixels = (1400, 1000)                    # wide aspect ratio
    plot_xy.color_by = "material"
    plot_xy.colors = mat_colors

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_rz, plot_xy])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  1. aries_at_geometry_rz.png -- RZ poloidal cross-section")
    print("     Shows radial build: plasma -> FW -> blanket -> shield -> VV -> TF")
    print("  2. aries_at_geometry_xy.png -- XY toroidal midplane")
    print("     Shows 22.5-degree sector with reflective boundaries")
    print()
    print("Material colour key:")
    print("  Tan:         SiC composite (first wall)")
    print("  Orange:      SiC/LiPb blanket (breeding zone)")
    print("  Dark orange: Back wall / manifold")
    print("  Dark grey:   WC/H2O shield")
    print("  Blue:        F82H ferritic steel (vacuum vessel)")
    print("  Green:       TF coil winding pack (SS316/Cu)")


if __name__ == "__main__":
    main()
