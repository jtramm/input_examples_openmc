#!/usr/bin/env python3
"""
Simplified ITER 40-Degree Sector -- Geometry Visualization
=============================================================

Produces two geometry plots of the ITER sector model:

  1. **RZ cross-section** (xz plane at y=0): Shows the major radial build
     through the midplane, from the plasma core outward through the first
     wall, blanket, vacuum vessel, and TF coils. This is the standard view
     for tokamak shielding analysis, showing the attenuation path that
     14.1 MeV neutrons must traverse. The toroidal shells appear as nested
     ellipses (circles for this simplified model).

  2. **XY midplane** (xy plane at z=0): Shows the 40-degree sector from
     above, with the two reflective sector boundary planes visible as
     straight edges. The concentric toroidal shells appear as arcs. This
     view confirms the sector geometry and reflective boundary conditions.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    iter_sector_rz.png    -- RZ cross-section at y=0
    iter_sector_xy.png    -- XY midplane at z=0
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
    beryllium = mat_by_name["Beryllium armor"]
    copper = mat_by_name["Copper heat sink"]
    ss316ln = mat_by_name["SS316LN"]
    blanket_fill = mat_by_name["Blanket fill (SS316LN+H2O)"]
    tf_winding = mat_by_name["TF winding pack (homogenized)"]

    # Material colours chosen for visual clarity and physical association
    color_map = {
        beryllium:    (180, 180, 180),   # light grey (Be metal)
        copper:       (184, 115, 51),    # copper-brown
        ss316ln:      (100, 130, 170),   # steel-blue
        blanket_fill: (140, 180, 140),   # muted green (steel+water mix)
        tf_winding:   (160, 100, 160),   # purple (superconductor)
    }

    # =========================================================================
    # Plot 1: RZ cross-section (XZ plane at y=0)
    # =========================================================================
    # This is the most informative view for a tokamak shielding model.
    # The x-axis shows major radius (R), the z-axis shows vertical position.
    # The nested toroidal shells appear as concentric circles (or ellipses
    # for elongated models) centred at R = R0 = 620 cm.
    #
    # View window: R from 200 to 1050 cm, Z from -450 to +450 cm
    # This captures the full radial build from the inner bore to the outer
    # TF coil case, with some margin.
    plot_rz = openmc.Plot()
    plot_rz.filename = "iter_sector_rz"
    plot_rz.basis = "xz"
    plot_rz.origin = (620.0, 0.0, 0.0)     # centred at magnetic axis
    plot_rz.width = (850.0, 900.0)          # R and Z extent
    plot_rz.pixels = (1700, 1800)
    plot_rz.color_by = "material"
    plot_rz.colors = color_map

    # =========================================================================
    # Plot 2: XY midplane (XY plane at z=0)
    # =========================================================================
    # The top-down view showing the 40-degree sector. The sector boundaries
    # appear as straight radial lines from the origin. The concentric
    # toroidal shells appear as arcs of circles.
    #
    # The view is centred at the origin to show the full sector geometry.
    plot_xy = openmc.Plot()
    plot_xy.filename = "iter_sector_xy"
    plot_xy.basis = "xy"
    plot_xy.origin = (500.0, 250.0, 0.0)   # offset to show sector clearly
    plot_xy.width = (1200.0, 1200.0)
    plot_xy.pixels = (1600, 1600)
    plot_xy.color_by = "material"
    plot_xy.colors = color_map

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_rz, plot_xy])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  iter_sector_rz.png  -- RZ cross-section (xz plane at y=0)")
    print("    Nested toroidal shells: grey=Be, brown=Cu, blue=SS316LN,")
    print("    green=blanket/VV fill, purple=TF winding pack")
    print()
    print("  iter_sector_xy.png  -- XY midplane (xy plane at z=0)")
    print("    40-degree sector with reflective boundaries")
    print("    Concentric arcs show radial build components")


if __name__ == "__main__":
    main()
