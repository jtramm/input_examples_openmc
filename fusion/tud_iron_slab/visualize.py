#!/usr/bin/env python3
"""
TUD Iron Slab Benchmark -- Geometry Visualization
===================================================

Produces XZ cross-section plots for each of the three TUD iron slab
configurations (A0, A1, A2). The plots slice through y=0 (the midplane)
for A0, and through the gap centre for A1/A2, showing how the vertical
gap appears in the slab cross-section.

For configurations with a gap (A1, A2), an additional plot is produced
slicing through the gap centre in the XY plane to show the gap position
relative to the beam axis.

Usage:
    python visualize.py

Output:
    tud_iron_slab_A0_xz.png   -- A0 solid slab, XZ midplane
    tud_iron_slab_A1_xz.png   -- A1 gap at 10 cm, XZ slice through gap
    tud_iron_slab_A2_xz.png   -- A2 gap at 20 cm, XZ slice through gap
    tud_iron_slab_A1_xy.png   -- A1 gap at 10 cm, XY slice showing gap
    tud_iron_slab_A2_xy.png   -- A2 gap at 20 cm, XY slice showing gap

Note:
    model.py must be run for each configuration before visualizing.
    This script generates all three configurations sequentially.
"""

import argparse
import subprocess
import sys
import openmc


def build_model_for_config(config_name):
    """Run model.py with the specified configuration to generate XML files.

    Parameters
    ----------
    config_name : str
        One of 'A0', 'A1', 'A2'.
    """
    print(f"\n--- Generating XML files for configuration {config_name} ---")
    result = subprocess.run(
        [sys.executable, "model.py", "--config", config_name, "--particles", "1000"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: model.py failed for {config_name}:")
        print(result.stderr)
        sys.exit(1)
    print(f"  XML files generated for {config_name}")


def create_xz_plot(config_name, materials, geometry):
    """Create an XZ cross-section plot for the given configuration.

    The XZ plane (y=0) shows the slab cross-section along the beam axis (x)
    and the vertical direction (z). For the solid A0 configuration, this
    shows a uniform iron rectangle. For A1/A2, the gap is not visible in
    this view since the gap is offset in y.

    Parameters
    ----------
    config_name : str
        Configuration name for the filename.
    materials : openmc.Materials
        Materials loaded from XML.
    geometry : openmc.Geometry
        Geometry loaded from XML.
    """
    plot = openmc.Plot()
    plot.filename = f"tud_iron_slab_{config_name}_xz"
    plot.basis = "xz"                              # XZ plane (beam axis vs vertical)
    plot.origin = (15.0, 0.0, 0.0)                # centre view on slab midpoint
    plot.width = (80.0, 120.0)                     # field of view: 80 cm along x, 120 cm along z
    plot.pixels = (800, 1200)                       # image resolution
    plot.color_by = "material"                      # colour by material assignment

    # Iron shown as a steel-grey colour
    plot.colors = {
        materials[0]: (160, 160, 170),              # grey for iron
    }

    return plot


def create_xy_plot(config_name, materials, geometry):
    """Create an XY cross-section plot for gap configurations.

    The XY plane (z=0) shows the slab cross-section along the beam axis (x)
    and the horizontal transverse direction (y). For A1/A2, this view clearly
    shows the vertical gap as a void channel through the iron slab.

    Parameters
    ----------
    config_name : str
        Configuration name for the filename.
    materials : openmc.Materials
        Materials loaded from XML.
    geometry : openmc.Geometry
        Geometry loaded from XML.
    """
    plot = openmc.Plot()
    plot.filename = f"tud_iron_slab_{config_name}_xy"
    plot.basis = "xy"                              # XY plane (beam axis vs horizontal)
    plot.origin = (15.0, 0.0, 0.0)                # centre view on slab midpoint
    plot.width = (80.0, 120.0)                     # field of view: 80 cm along x, 120 cm along y
    plot.pixels = (800, 1200)                       # image resolution
    plot.color_by = "material"                      # colour by material assignment

    # Iron shown as a steel-grey colour
    plot.colors = {
        materials[0]: (160, 160, 170),              # grey for iron
    }

    return plot


def main():
    # =========================================================================
    # Configuration definitions
    # =========================================================================
    # List of all three configurations to visualize
    configs = ["A0", "A1", "A2"]

    # Configurations with gaps (need an additional XY plot to show the gap)
    gap_configs = ["A1", "A2"]

    for config_name in configs:
        # -----------------------------------------------------------------
        # Step 1: Generate XML files by running model.py for this config
        # -----------------------------------------------------------------
        build_model_for_config(config_name)

        # -----------------------------------------------------------------
        # Step 2: Load geometry and materials from the generated XML
        # -----------------------------------------------------------------
        materials = openmc.Materials.from_xml("materials.xml")
        geometry = openmc.Geometry.from_xml("geometry.xml", materials=materials)

        # -----------------------------------------------------------------
        # Step 3: Create XZ cross-section plot (all configurations)
        # -----------------------------------------------------------------
        plots_list = []

        xz_plot = create_xz_plot(config_name, materials, geometry)
        plots_list.append(xz_plot)

        # -----------------------------------------------------------------
        # Step 4: For gap configs, also create an XY plot showing the gap
        # -----------------------------------------------------------------
        if config_name in gap_configs:
            xy_plot = create_xy_plot(config_name, materials, geometry)
            plots_list.append(xy_plot)

        # -----------------------------------------------------------------
        # Step 5: Export plots and generate images
        # -----------------------------------------------------------------
        plots = openmc.Plots(plots_list)
        plots.export_to_xml()

        # Generate the plot images using OpenMC's geometry plotter
        openmc.plot_geometry()

        # Print summary of generated files
        print(f"\nPlots generated for {config_name}:")
        print(f"  - tud_iron_slab_{config_name}_xz.png (XZ cross-section)")
        if config_name in gap_configs:
            print(f"  - tud_iron_slab_{config_name}_xy.png (XY cross-section, showing gap)")

    print("\nAll visualization plots complete.")


if __name__ == "__main__":
    main()
