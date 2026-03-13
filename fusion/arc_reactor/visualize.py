#!/usr/bin/env python3
"""
ARC Reactor Fusion Neutronics Benchmark -- Geometry Visualization
==================================================================

Produces cross-section plots of the ARC reactor toroidal sector geometry.
Two views are generated:

  1. RZ (poloidal) cross-section: a vertical slice through the tokamak
     at the sector midplane, showing the concentric annular layers
     (plasma, first wall, blanket, vacuum vessel, shield, coil).

  2. XY (toroidal) cross-section: a horizontal slice at the midplane
     (z = 0), showing the 20-degree sector with the toroidal curvature
     of each layer visible.

Usage:
    python visualize.py

Prerequisite:
    Run model.py first to generate the XML files (or at least export
    the model).

Output:
    arc_reactor_rz.png   -- poloidal (RZ) cross-section
    arc_reactor_xy.png   -- toroidal (XY) midplane cross-section
"""

import openmc


def main():
    # =========================================================================
    # Load the model geometry and materials from the model XML
    # =========================================================================
    # The model.py script exports a model_xml file that contains all
    # geometry, materials, settings, and tallies.
    model = openmc.Model.from_model_xml("model.xml")
    materials = model.materials
    geometry = model.geometry

    # Build a material colour map for visual clarity.
    # Assign distinctive colours to each material so the layers are
    # easy to distinguish in the plots.
    color_map = {}
    for mat in materials:
        if "FLiBe" in mat.name:
            # Green for molten salt blanket (evocative of fluoride salts)
            color_map[mat] = (100, 200, 100)
        elif "Inconel" in mat.name:
            # Silver-grey for the nickel superalloy
            color_map[mat] = (180, 180, 190)
        elif "Borated" in mat.name:
            # Dark blue-grey for the borated steel shield
            color_map[mat] = (80, 80, 120)
        elif "Coil" in mat.name:
            # Copper-orange for the TF coil composite
            color_map[mat] = (200, 140, 60)

    # =========================================================================
    # Plot 1: RZ (poloidal) cross-section
    # =========================================================================
    # This view slices through the tokamak at a fixed toroidal angle,
    # showing the vertical (Z) vs radial (R) structure. We use the
    # xz-plane (y = 0), which is one of the sector boundary planes.
    #
    # The plot is centred on the magnetic axis (x = R_MAJOR = 330 cm,
    # z = 0) and wide enough to show the full radial build from the
    # inboard side (R_MAJOR - R_COIL_OUTER) to the outboard side
    # (R_MAJOR + R_COIL_OUTER).
    rz_plot = openmc.Plot()
    rz_plot.filename = "arc_reactor_rz"
    rz_plot.basis = "xz"
    rz_plot.origin = (330.0, 0.0, 0.0)     # centre on magnetic axis
    rz_plot.width = (600.0, 550.0)          # wide enough for full cross-section
    rz_plot.pixels = (1200, 1100)           # high resolution
    rz_plot.color_by = "material"
    rz_plot.colors = color_map

    # =========================================================================
    # Plot 2: XY (toroidal) midplane cross-section
    # =========================================================================
    # This view slices horizontally at z = 0, showing the top-down view
    # of the 20-degree sector. The toroidal curvature of each layer is
    # visible as concentric circular arcs.
    #
    # The sector spans from the xz-plane (phi = 0) to phi = 20 degrees.
    # We centre the view on the sector midpoint.
    xy_plot = openmc.Plot()
    xy_plot.filename = "arc_reactor_xy"
    xy_plot.basis = "xy"
    xy_plot.origin = (310.0, 60.0, 0.0)    # offset to centre the sector view
    xy_plot.width = (650.0, 350.0)          # wide enough for the sector arc
    xy_plot.pixels = (1300, 700)
    xy_plot.color_by = "material"
    xy_plot.colors = color_map

    # =========================================================================
    # Export and generate the plots
    # =========================================================================
    plots = openmc.Plots([rz_plot, xy_plot])
    plots.export_to_xml()

    # Generate the plot images using OpenMC's built-in plotter.
    # This calls the C++ plotting routines via the openmc executable.
    openmc.plot_geometry()

    print("Geometry plots saved:")
    print("  arc_reactor_rz.png  -- Poloidal (RZ) cross-section at y = 0")
    print("    Shows: concentric torus layers (plasma void, FW, blanket,")
    print("           vacuum vessel, shield, TF coil)")
    print("  arc_reactor_xy.png  -- Toroidal (XY) midplane cross-section")
    print("    Shows: 20-degree sector with reflective boundaries,")
    print("           toroidal curvature of each layer")
    print()
    print("Material colour key:")
    print("  Green        -- FLiBe (Li2BeF4) molten salt blanket")
    print("  Silver-grey  -- Inconel 718 (first wall and vacuum vessel)")
    print("  Dark blue    -- Borated steel neutron shield")
    print("  Copper-orange -- TF coil composite (Cu/SS/HTS)")
    print("  Black/void   -- Plasma region and external void")


if __name__ == "__main__":
    main()
