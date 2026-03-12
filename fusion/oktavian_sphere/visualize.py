#!/usr/bin/env python3
"""
OKTAVIAN Iron Sphere -- Geometry Visualization
===============================================

Generates cross-section plots of the iron sphere geometry using OpenMC's
built-in plotting capabilities.  Two plots are produced:

  1. XY cross-section through the origin  (z = 0 plane)
  2. XZ cross-section through the origin  (y = 0 plane)

Materials are colored to distinguish iron from void.

Prerequisites
-------------
Run model.py first to generate the XML input files.
"""

import openmc


def main():
    # ------------------------------------------------------------------
    # Plot 1: XY cross-section (looking down the z-axis)
    # ------------------------------------------------------------------
    xy_plot = openmc.Plot()
    xy_plot.filename = "plot_xy"
    xy_plot.basis = "xy"                  # slice in the xy-plane
    xy_plot.origin = (0.0, 0.0, 0.0)     # centered at the origin
    xy_plot.width = (500.0, 500.0)        # 500 cm x 500 cm view
    xy_plot.pixels = (800, 800)           # image resolution
    xy_plot.color_by = "material"         # color cells by material

    # ------------------------------------------------------------------
    # Plot 2: XZ cross-section (looking down the y-axis)
    # ------------------------------------------------------------------
    xz_plot = openmc.Plot()
    xz_plot.filename = "plot_xz"
    xz_plot.basis = "xz"                  # slice in the xz-plane
    xz_plot.origin = (0.0, 0.0, 0.0)
    xz_plot.width = (500.0, 500.0)
    xz_plot.pixels = (800, 800)
    xz_plot.color_by = "material"

    # ------------------------------------------------------------------
    # Generate the plots
    # ------------------------------------------------------------------
    plots = openmc.Plots([xy_plot, xz_plot])
    plots.export_to_xml()

    # Run OpenMC in plot mode (reads geometry.xml + materials.xml + plots.xml)
    openmc.plot_geometry()

    print("Plots generated:")
    print("  plot_xy.png  -- XY cross-section through origin")
    print("  plot_xz.png  -- XZ cross-section through origin")


if __name__ == "__main__":
    main()
