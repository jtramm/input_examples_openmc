#!/usr/bin/env python3
"""
CFETR 22.5-Degree Toroidal Sector -- Geometry Visualization
==============================================================

Produces cross-section plots of the CFETR sector showing the concentric
toroidal shell structure: plasma, first wall, breeding blanket, back support,
vacuum vessel, thermal shield, and TF coil.

Two views are generated:
  1. XZ cross-section (y=0 midplane): shows the R-Z view of the tokamak
     sector, which is the standard way to visualize a tokamak cross-section.
     The toroidal shells appear as nested annuli cut by the midplane.
  2. XY cross-section (z=0 midplane): shows the top-down view of the
     22.5-degree sector, with the two reflective boundary planes visible.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    cfetr_geometry_rz.png   (R-Z cross-section, standard tokamak view)
    cfetr_geometry_xy.png   (top-down view showing sector boundaries)
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
    # Colours chosen to visually distinguish each radial layer:
    #   - CLF-1 steel:   steel grey (first wall, structural)
    #   - WCCB blanket:  green (breeding region, the key component)
    #   - Back support:  tan/brown (support structure)
    #   - SS316:         blue (vacuum vessel walls)
    #   - VV fill:       light blue (VV shielding)
    #   - SS304:         purple (thermal shield)
    #   - TF coil:       orange (magnet casing)
    mat_colors = {}
    for mat in materials:
        name = mat.name.lower()
        if "clf-1" in name and "back" not in name:
            mat_colors[mat] = (160, 160, 170)      # steel grey
        elif "wccb" in name or "homogenized" in name:
            mat_colors[mat] = (50, 180, 50)         # green
        elif "back support" in name:
            mat_colors[mat] = (180, 140, 80)        # tan
        elif "vv fill" in name:
            mat_colors[mat] = (130, 170, 220)       # light blue
        elif "ss316" in name:
            mat_colors[mat] = (70, 70, 200)         # blue
        elif "ss304" in name or "thermal" in name:
            mat_colors[mat] = (150, 80, 180)        # purple
        elif "tf coil" in name:
            mat_colors[mat] = (220, 150, 50)        # orange
        elif "li4sio4" in name:
            mat_colors[mat] = (80, 200, 80)         # bright green
        elif "beryllium" in name:
            mat_colors[mat] = (200, 200, 200)       # light grey
        elif "water" in name:
            mat_colors[mat] = (100, 150, 255)       # water blue
        elif "copper" in name:
            mat_colors[mat] = (200, 120, 50)        # copper

    # =========================================================================
    # Plot 1: R-Z cross-section (XZ plane at y=0)
    # =========================================================================
    # This is the standard tokamak cross-section view showing the poloidal
    # plane. The torus cross-sections appear as nested circles (since we
    # use circular cross-sections with b=c).
    #
    # The major radius R0=560 cm means the torus tube centres are at
    # x = +/- 560 cm. We focus on the outboard side (positive x).
    # The view spans from the inboard side (x ~ 560-252 = 308 cm) to the
    # outboard side (x ~ 560+252 = 812 cm).
    plot_rz = openmc.Plot()
    plot_rz.filename = "cfetr_geometry_rz"
    plot_rz.basis = "xz"
    plot_rz.origin = (560.0, 0.0, 0.0)            # centred on magnetic axis
    plot_rz.width = (600.0, 600.0)                  # wide enough to show full cross-section
    plot_rz.pixels = (1200, 1200)
    plot_rz.color_by = "material"
    plot_rz.colors = mat_colors

    # =========================================================================
    # Plot 2: Top-down view (XY plane at z=0)
    # =========================================================================
    # Shows the 22.5-degree sector from above. The sector boundaries
    # (reflective planes) are visible as straight edges. The toroidal
    # shells appear as concentric arcs.
    plot_xy = openmc.Plot()
    plot_xy.filename = "cfetr_geometry_xy"
    plot_xy.basis = "xy"
    plot_xy.origin = (400.0, 180.0, 0.0)           # offset to show sector well
    plot_xy.width = (900.0, 900.0)
    plot_xy.pixels = (1200, 1200)
    plot_xy.color_by = "material"
    plot_xy.colors = mat_colors

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_rz, plot_xy])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  1. cfetr_geometry_rz.png -- R-Z cross-section (standard tokamak view)")
    print("     Shows nested toroidal shells: FW, blanket, BSS, VV, shield, TF")
    print("  2. cfetr_geometry_xy.png -- Top-down view (XY midplane)")
    print("     Shows 22.5-degree sector with reflective boundaries")
    print()
    print("Material colour key:")
    print("  Steel grey:   CLF-1 RAFM steel (first wall)")
    print("  Green:        WCCB homogenized blanket (breeder)")
    print("  Tan:          Back support structure (CLF-1 + H2O)")
    print("  Blue:         SS316LN (VV walls)")
    print("  Light blue:   VV fill (SS316 + H2O)")
    print("  Purple:       SS304 (thermal shield)")
    print("  Orange:       TF coil (SS316 + Cu composite)")


if __name__ == "__main__":
    main()
