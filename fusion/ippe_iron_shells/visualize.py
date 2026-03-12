#!/usr/bin/env python3
"""
================================================================================
IPPE Iron Spherical Shell Transmission Benchmark - Geometry Visualization
================================================================================

Generates XY cross-section plots for all five IPPE iron shell configurations.
Each plot shows the three regions:
  - Inner void (vacuum cavity where the source is located)
  - Iron shell wall (the transmission medium)
  - Outer void (vacuum outside the shell)

The plots are saved as PNG files and can be used for geometry verification.

Usage:
    python visualize.py             # Plot all 5 shells
    python visualize.py --shell 3   # Plot only shell 3
================================================================================
"""

import argparse
import openmc

# Import the model builder from model.py
from model import build_model, SHELL_DATA


def plot_shell(shell_number: int):
    """
    Create an XY cross-section plot for the specified shell configuration.

    Parameters
    ----------
    shell_number : int
        Shell configuration number (1-5).
    """
    # Build the model (with minimal particles since we only need geometry)
    model = build_model(shell_number, particles=100)

    # Retrieve shell dimensions for setting the plot bounds
    shell = SHELL_DATA[shell_number]
    r_inner = shell["inner_radius"]
    r_outer = r_inner + shell["thickness"]

    # Set the plot width to show the full shell plus some margin
    # We want to see the shell clearly with some void on each side
    plot_width = 2.0 * (r_outer + 10.0)  # 10 cm margin beyond outer surface

    # Create an XY cross-section plot centered at the origin
    plot = openmc.Plot()
    plot.filename = f"shell_{shell_number}_xy"
    plot.basis = "xy"                      # XY plane (z=0)
    plot.origin = (0.0, 0.0, 0.0)         # Center at origin
    plot.width = (plot_width, plot_width)  # Plot dimensions in cm
    plot.pixels = (800, 800)               # Resolution
    plot.color_by = "material"             # Color cells by material

    # Assign colors: void regions are white, iron is steel-blue
    # (Void cells have no material, so they default to white/background)
    plot.background = (255, 255, 255)  # White background for void

    # Add plot to the model and generate the image
    plots = openmc.Plots([plot])
    model.plots = plots

    # Export and generate
    model.export_to_xml()
    openmc.plot_geometry()

    print(f"Shell {shell_number}: Plot saved as '{plot.filename}.png'")
    print(f"  Inner radius: {r_inner:.1f} cm")
    print(f"  Outer radius: {r_outer:.1f} cm")
    print(f"  Plot width:   {plot_width:.1f} cm")


def main():
    """Generate geometry visualization plots."""

    parser = argparse.ArgumentParser(
        description="Visualize IPPE iron shell geometry (XY cross-sections)."
    )
    parser.add_argument(
        "--shell", type=int, default=None, choices=[1, 2, 3, 4, 5],
        help="Shell to plot (1-5). Default: plot all 5 shells."
    )
    args = parser.parse_args()

    if args.shell is not None:
        # Plot a single shell
        plot_shell(args.shell)
    else:
        # Plot all five shells
        print("Generating XY cross-section plots for all 5 shell configurations...\n")
        for shell_num in range(1, 6):
            plot_shell(shell_num)
            print()

    print("Done. Geometry plots generated successfully.")


if __name__ == "__main__":
    main()
