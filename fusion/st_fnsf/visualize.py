#!/usr/bin/env python3
"""
ST-FNSF Benchmark -- Geometry Visualization
=============================================

Produces an RZ (xz) cross-section plot of the ST-FNSF toroidal sector,
slicing through the midplane at y=0. This view clearly shows the
asymmetric inboard/outboard structure that defines spherical tokamak
geometry:

  - The compact center column (Cu post + WC shield) on the inboard side
  - The plasma region spanning from R=70 to R=270 cm
  - The full blanket/shield/VV stack on the outboard side only
  - The absence of breeding blanket on the inboard side

The plot highlights the key neutronics challenge: the inboard side has
only 20 cm of shielding between the plasma and the copper center post.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files.

Output:
    st_fnsf_geometry_rz.png   -- RZ cross-section at y=0
    st_fnsf_geometry_xy.png   -- XY cross-section at z=0 (plan view)
"""

import openmc


def main():
    # =========================================================================
    # Load the model geometry and materials from XML files
    # =========================================================================
    materials = openmc.Materials.from_xml("materials.xml")
    geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

    # =========================================================================
    # Define material colours
    # =========================================================================
    # Use distinctive colours that reflect the physical nature of each material
    mat_colors = {}
    for mat in materials:
        if "Copper" in mat.name:
            mat_colors[mat] = (184, 115, 51)     # copper-brown
        elif "inboard shield" in mat.name.lower():
            mat_colors[mat] = (100, 100, 100)     # dark grey (WC)
        elif "F82H" in mat.name:
            mat_colors[mat] = (180, 180, 200)     # light steel blue
        elif "DCLL" in mat.name:
            mat_colors[mat] = (50, 120, 200)      # blue (LiPb)
        elif "outboard shield" in mat.name.lower():
            mat_colors[mat] = (120, 120, 120)     # medium grey
        elif "SS316" in mat.name:
            mat_colors[mat] = (200, 200, 180)     # pale yellow-grey

    # =========================================================================
    # Plot 1: RZ cross-section (xz plane at y=0)
    # =========================================================================
    # This is the canonical view for tokamak neutronics. The x-axis shows
    # the radial direction (R) and the z-axis shows the vertical direction.
    # The asymmetry between inboard and outboard is immediately visible.
    #
    # View window: R from -50 to 450 cm, Z from -350 to 350 cm
    # (showing center column at left, outboard components at right)
    plot_rz = openmc.Plot()
    plot_rz.filename = "st_fnsf_geometry_rz"
    plot_rz.basis = "xz"
    plot_rz.origin = (200.0, 0.0, 0.0)       # centered on tokamak midplane
    plot_rz.width = (500.0, 700.0)            # R range and Z range in cm
    plot_rz.pixels = (1000, 1400)             # high resolution
    plot_rz.color_by = "material"
    plot_rz.colors = mat_colors

    # =========================================================================
    # Plot 2: XY cross-section (plan view at z=0)
    # =========================================================================
    # Shows the 30-degree toroidal sector from above. The sector boundaries
    # (reflective planes) are visible, and the concentric toroidal shells
    # appear as arcs.
    plot_xy = openmc.Plot()
    plot_xy.filename = "st_fnsf_geometry_xy"
    plot_xy.basis = "xy"
    plot_xy.origin = (200.0, 60.0, 0.0)      # centered on sector midpoint
    plot_xy.width = (500.0, 250.0)            # x and y range in cm
    plot_xy.pixels = (1000, 500)
    plot_xy.color_by = "material"
    plot_xy.colors = mat_colors

    # =========================================================================
    # Export and generate plots
    # =========================================================================
    plots = openmc.Plots([plot_rz, plot_xy])
    plots.export_to_xml()

    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  st_fnsf_geometry_rz.png -- RZ cross-section at y=0 (midplane)")
    print("    - Left: compact center column (Cu post + WC shield)")
    print("    - Center: plasma (void)")
    print("    - Right: outboard blanket/shield/VV stack")
    print("    - Note: NO blanket on inboard side (ST constraint)")
    print()
    print("  st_fnsf_geometry_xy.png -- XY cross-section at z=0 (plan view)")
    print("    - 30-degree sector with reflective boundaries")
    print("    - Concentric arcs: FW, blanket, shield, VV")


if __name__ == "__main__":
    main()
