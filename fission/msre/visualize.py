#!/usr/bin/env python3
"""
===============================================================================
MSRE Visualization Script
===============================================================================

Generates cross-section plots of the MSRE OpenMC model:
  1. XY cross section (radial view at core midplane)
  2. XZ cross section (axial view through the center)

Prerequisites:
  Run model.py first to generate the model XML files:
    python model.py --model homogeneous
  or
    python model.py --model channel

Usage:
  python visualize.py

The script detects whether the homogeneous or channel model was built
by checking the geometry extents, and adjusts the plot dimensions
accordingly.
"""

import openmc


def create_plots():
    """Create XY and XZ cross-section plots of the MSRE model.

    Cross-section plots slice through the geometry at a specified point
    and color each pixel according to the material (or cell) at that
    location. They are essential for verifying that the geometry was
    constructed correctly before running a full Monte Carlo simulation.
    """
    model = openmc.Model.from_model_xml()

    # Detect model type by looking at the bounding box.
    # The homogeneous model has a vessel radius of ~76 cm,
    # while the channel model is only ~2.5 cm half-pitch.
    all_cells = model.geometry.get_all_cells()
    is_channel = len(all_cells) <= 3  # channel model has 2-3 cells

    if is_channel:
        # Channel unit cell model
        # The unit cell is 5.08 cm x 5.08 cm x 162.56 cm
        xy_width = 6.0   # slightly larger than the 5.08 cm pitch
        xz_width = 6.0
        xz_height = 170.0

        print("Detected CHANNEL unit cell model")
        print(f"  XY plot: {xy_width} x {xy_width} cm at z=0")
        print(f"  XZ plot: {xz_width} x {xz_height} cm at y=0")
    else:
        # Homogeneous full-core model
        # The vessel outer radius is ~76.2 cm, total height ~244 cm
        xy_width = 180.0
        xz_width = 180.0
        xz_height = 280.0

        print("Detected HOMOGENEOUS full-core model")
        print(f"  XY plot: {xy_width} x {xy_width} cm at z=0")
        print(f"  XZ plot: {xz_width} x {xz_height} cm at y=0")

    plots = openmc.Plots()

    # --- XY cross section at core midplane (z = 0) ---
    # This shows the radial layout: for the homogeneous model you see
    # concentric rings (core, downcomer, vessel wall). For the channel
    # model you see the fuel channel in the graphite block.
    xy_plot = openmc.Plot(name="XY Cross Section (z=0)")
    xy_plot.basis = "xy"
    xy_plot.origin = (0.0, 0.0, 0.0)
    xy_plot.width = (xy_width, xy_width)
    xy_plot.pixels = (800, 800)
    xy_plot.color_by = "material"
    plots.append(xy_plot)

    # --- XZ cross section through center (y = 0) ---
    # This shows the axial layout: for the homogeneous model you see
    # the core, upper/lower plenums, and vessel wall. For the channel
    # model you see the long fuel channel running vertically.
    xz_plot = openmc.Plot(name="XZ Cross Section (y=0)")
    xz_plot.basis = "xz"
    xz_plot.origin = (0.0, 0.0, 0.0)
    xz_plot.width = (xz_width, xz_height)
    xz_plot.pixels = (600, 1000)
    xz_plot.color_by = "material"
    plots.append(xz_plot)

    model.plots = plots
    model.plots.export_to_xml()

    # Generate the plot images
    openmc.plot_geometry()

    print()
    print("Plot files generated:")
    print("  1_xy_cross_section_z0.png  - Radial cross section at midplane")
    print("  2_xz_cross_section_y0.png  - Axial cross section through center")
    print()
    print("Open these PNG files to visually verify the geometry.")


if __name__ == "__main__":
    create_plots()
