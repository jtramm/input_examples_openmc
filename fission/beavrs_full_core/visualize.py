"""
BEAVRS Full-Core Visualization Script
======================================

Generates slice plots of the BEAVRS full-core PWR model for visual
verification of the geometry. Produces:

1. XY slice (radial cross-section) colored by material -- shows all 193
   fuel assemblies, the core baffle, barrel, and water reflector.
2. XY slice colored by cell -- useful for verifying lattice structure.
3. Zoomed XY slice of a single assembly to verify pin layout.

The plots are exported as PNG images using OpenMC's built-in plotting
capability, which ray-traces through the geometry and colors each pixel
according to the material or cell at that point.

Usage:
    python visualize.py              # Generate and display plots
    python visualize.py --no-show    # Generate plots without displaying
"""

import argparse
import os
import sys

import openmc

# Import the model builder
from model import build_model, ASSEMBLY_PITCH, BARREL_OR


def create_plots():
    """
    Create OpenMC plot objects for the BEAVRS full-core model.

    Returns a list of Plot objects that can be added to a Plots collection
    and rendered.
    """

    plots = []

    # -------------------------------------------------------------------------
    # Plot 1: Full-core radial cross-section (XY plane), colored by material
    # -------------------------------------------------------------------------
    # This is the "classic" reactor cross-section view. It shows the circular
    # arrangement of square fuel assemblies within the cylindrical barrel.
    # Materials are color-coded so you can distinguish fuel enrichment zones,
    # cladding, coolant, and structural steel.

    full_core_mat = openmc.SlicePlot()
    full_core_mat.basis = 'xy'
    full_core_mat.origin = (0.0, 0.0, 0.0)
    # Width covers from barrel edge to barrel edge with some margin
    plot_width = 2.0 * BARREL_OR + 10.0  # ~400 cm total
    full_core_mat.width = (plot_width, plot_width)
    # High resolution for clear pin-level detail
    full_core_mat.pixels = (4000, 4000)
    full_core_mat.color_by = 'material'
    full_core_mat.filename = 'beavrs_full_core_xy_material'

    # Custom colors for better visual distinction
    # (OpenMC will assign default colors if we don't specify)
    plots.append(full_core_mat)

    # -------------------------------------------------------------------------
    # Plot 2: Full-core XY slice colored by cell
    # -------------------------------------------------------------------------
    # Coloring by cell helps verify that each assembly position has the
    # correct lattice structure and that guide tube / instrument tube
    # positions are correct.

    full_core_cell = openmc.SlicePlot()
    full_core_cell.basis = 'xy'
    full_core_cell.origin = (0.0, 0.0, 0.0)
    full_core_cell.width = (plot_width, plot_width)
    full_core_cell.pixels = (4000, 4000)
    full_core_cell.color_by = 'cell'
    full_core_cell.filename = 'beavrs_full_core_xy_cell'
    plots.append(full_core_cell)

    # -------------------------------------------------------------------------
    # Plot 3: Zoomed view of a single assembly (center assembly)
    # -------------------------------------------------------------------------
    # This shows the 17x17 pin layout in detail. You should be able to see:
    # - 264 fuel pins (small circles with fuel, gap, and clad rings)
    # - 24 guide tubes (larger circles, water-filled)
    # - 1 instrument tube (center, slightly different diameter)

    single_assembly = openmc.SlicePlot()
    single_assembly.basis = 'xy'
    single_assembly.origin = (0.0, 0.0, 0.0)  # Center assembly
    single_assembly.width = (ASSEMBLY_PITCH, ASSEMBLY_PITCH)
    single_assembly.pixels = (1000, 1000)
    single_assembly.color_by = 'material'
    single_assembly.filename = 'beavrs_single_assembly_xy'
    plots.append(single_assembly)

    # -------------------------------------------------------------------------
    # Plot 4: Zoomed view of a single fuel pin
    # -------------------------------------------------------------------------
    # Shows the concentric ring structure: fuel pellet, helium gap, cladding,
    # and surrounding moderator.

    single_pin = openmc.SlicePlot()
    single_pin.basis = 'xy'
    single_pin.origin = (0.0, 0.0, 0.0)  # Center pin of center assembly
    single_pin.width = (1.5, 1.5)  # Slightly larger than pin pitch
    single_pin.pixels = (500, 500)
    single_pin.color_by = 'material'
    single_pin.filename = 'beavrs_single_pin_xy'
    plots.append(single_pin)

    return plots


def main():
    """Generate and optionally display BEAVRS geometry plots."""

    parser = argparse.ArgumentParser(
        description='Generate visualization plots of the BEAVRS full-core model'
    )
    parser.add_argument('--no-show', action='store_true',
                        help='Generate plot files without attempting to display them')
    parser.add_argument('--low-res', action='store_true',
                        help='Use lower resolution (faster rendering)')
    args = parser.parse_args()

    # Build the model (needed for geometry definition)
    print("Building BEAVRS model for visualization...")
    model = build_model()

    # Create the plots
    plots = create_plots()

    # Optionally reduce resolution for faster rendering
    if args.low_res:
        for plot in plots:
            plot.pixels = (min(plot.pixels[0], 1000),
                          min(plot.pixels[1], 1000))
        print("Using low-resolution mode for faster rendering.")

    # Add plots to the model
    model.plots = openmc.Plots(plots)

    # Export model to XML (required for plotting)
    model.export_to_xml()

    # Run the plotter
    print("\nGenerating geometry plots...")
    print("This may take a few minutes for high-resolution full-core plots.")
    openmc.plot_geometry()

    # List generated files
    print("\nGenerated plot files:")
    for plot in plots:
        filename = f'{plot.filename}.png'
        if os.path.exists(filename):
            size_mb = os.path.getsize(filename) / (1024 * 1024)
            print(f"  {filename} ({size_mb:.1f} MB)")
        else:
            # OpenMC may use .ppm format
            ppm_name = f'{plot.filename}.ppm'
            if os.path.exists(ppm_name):
                size_mb = os.path.getsize(ppm_name) / (1024 * 1024)
                print(f"  {ppm_name} ({size_mb:.1f} MB)")
            else:
                print(f"  {filename} (not found)")

    print("\nVisualization complete.")


if __name__ == '__main__':
    main()
