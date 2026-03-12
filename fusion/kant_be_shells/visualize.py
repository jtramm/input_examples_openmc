#!/usr/bin/env python3
"""
KANT Beryllium Spherical Shell Benchmark - Geometry Visualization
=================================================================

Produces XY and XZ cross-section plots for all three shell configurations.
"""

import openmc
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt
import numpy as np


def build_model(shell_num):
    """Build geometry for a given shell configuration (1, 2, or 3).

    Returns (geometry, outer_radius) so the caller can set plot bounds.
    """
    INNER_RADIUS = 5.0
    SHELL_CONFIG = {
        1: {"thickness": 5.0,  "outer_radius": 10.0},
        2: {"thickness": 10.0, "outer_radius": 15.0},
        3: {"thickness": 17.0, "outer_radius": 22.0},
    }
    config = SHELL_CONFIG[shell_num]
    outer_radius = config["outer_radius"]
    thickness = config["thickness"]

    # Material
    beryllium = openmc.Material(name="Beryllium")
    beryllium.add_nuclide("Be9", 1.0)
    beryllium.set_density("g/cm3", 1.85)
    materials = openmc.Materials([beryllium])

    # Surfaces
    inner_sphere = openmc.Sphere(r=INNER_RADIUS)
    outer_sphere = openmc.Sphere(r=outer_radius)
    boundary_sphere = openmc.Sphere(r=outer_radius + 10.0, boundary_type="vacuum")

    # Cells
    inner_void = openmc.Cell(name="Central void")
    inner_void.region = -inner_sphere

    be_shell = openmc.Cell(name=f"Be shell ({thickness} cm)")
    be_shell.fill = beryllium
    be_shell.region = +inner_sphere & -outer_sphere

    outer_void = openmc.Cell(name="Outer void")
    outer_void.region = +outer_sphere & -boundary_sphere

    universe = openmc.Universe(cells=[inner_void, be_shell, outer_void])
    geometry = openmc.Geometry(universe)

    return geometry, materials, outer_radius


def main():
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    for col, shell_num in enumerate([1, 2, 3]):
        geometry, materials, outer_radius = build_model(shell_num)

        # Determine plot extent (slightly larger than vacuum boundary)
        extent = outer_radius + 12.0

        for row, basis in enumerate(["xy", "xz"]):
            ax = axes[row, col]

            # Create OpenMC plot
            plot = openmc.Plot()
            plot.basis = basis
            plot.origin = (0.0, 0.0, 0.0)
            plot.width = (2 * extent, 2 * extent)
            plot.pixels = (400, 400)
            plot.color_by = "material"

            # Use the model to generate the plot image
            model = openmc.Model(geometry=geometry, materials=materials)
            plot_img = model.plot_geometry(plot)

            ax.imshow(plot_img, extent=[-extent, extent, -extent, extent])
            ax.set_xlabel("X [cm]" if "x" in basis else "Y [cm]")
            ax.set_ylabel("Y [cm]" if basis == "xy" else "Z [cm]")

            thickness = {1: 5, 2: 10, 3: 17}[shell_num]
            ax.set_title(
                f"Shell {shell_num} ({thickness} cm) - {basis.upper()} plane"
            )

    plt.suptitle(
        "KANT Beryllium Spherical Shell Benchmark - Geometry Cross Sections",
        fontsize=14, fontweight="bold"
    )
    plt.tight_layout()
    plt.savefig("geometry_cross_sections.png", dpi=150, bbox_inches="tight")
    print("Saved geometry_cross_sections.png")


if __name__ == "__main__":
    main()
