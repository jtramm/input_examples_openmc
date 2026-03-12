#!/usr/bin/env python3
"""
=============================================================================
SNAP-10A/2 Space Reactor Benchmark - Geometry Visualization
=============================================================================

Generates cross-section plots of the SNAP-10A/2 SCA-4B Case 8490a model:
  1. XY cross-section at z=0 (midplane) showing the 37-position array
     with 28 fuel elements and 9 water-filled vacant positions
  2. XZ cross-section at y=0 showing the axial extent of the core,
     vessel, grid plates, water tanks, and control cap

These plots use the OpenMC plotting API to render material-colored
cross sections of the geometry.
=============================================================================
"""

import openmc
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for headless environments
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# =============================================================================
# Define material colors for consistent visualization
# =============================================================================
# These colors help distinguish the different materials in the reactor model.
# Colors are chosen to be visually distinct and roughly match conventions:
#   - Red/orange for fuel (fissile material)
#   - Gray for structural metals (steel, Hastelloy)
#   - Blue for water (moderator/reflector)
#   - Green for aluminum (grid plates, cover)
#   - Yellow for the hydrogen barrier coating
material_colors = {
    "U-ZrH fuel (SCA-4)": "orangered",
    "H-barrier coating (Sm2O3/Al2O3/SiO2)": "gold",
    "Type 1100 Aluminum": "limegreen",
    "Light water": "royalblue",
    "Hastelloy N": "silver",
    "SS316 stainless steel": "dimgray",
}

# Build a color dictionary keyed by material ID
# First, load the model to get material IDs
model = openmc.Model.from_xml()
color_by_id = {}
for mat in model.materials:
    if mat.name in material_colors:
        color_by_id[mat] = material_colors[mat.name]


# =============================================================================
# Plot 1: XY Cross-Section at Midplane (z=0)
# =============================================================================
# This shows the arrangement of 28 fuel elements on a triangular pitch
# inside the cylindrical core vessel. The 9 vacant positions (29-37) are
# filled with water since the core is flooded.

print("Generating XY cross-section plot at z=0 (midplane)...")

# Create an OpenMC plot for the XY plane
xy_plot = openmc.Plot()
xy_plot.id = 101                        # Fixed ID to avoid collisions
xy_plot.filename = "plot_xy_full"
xy_plot.basis = "xy"                    # XY plane
xy_plot.origin = (0.0, 0.0, 0.0)       # Centered at midplane
xy_plot.width = (65.0, 65.0)           # Width in cm (enough to see tanks)
xy_plot.pixels = (800, 800)             # Resolution
xy_plot.color_by = "material"           # Color by material type
xy_plot.colors = color_by_id

# Also create a zoomed-in view of just the core
xy_plot_zoom = openmc.Plot()
xy_plot_zoom.id = 102
xy_plot_zoom.filename = "plot_xy_core"
xy_plot_zoom.basis = "xy"
xy_plot_zoom.origin = (0.0, 0.0, 0.0)
xy_plot_zoom.width = (28.0, 28.0)      # Just the core vessel
xy_plot_zoom.pixels = (800, 800)
xy_plot_zoom.color_by = "material"
xy_plot_zoom.colors = color_by_id


# =============================================================================
# Plot 2: XZ Cross-Section (side view)
# =============================================================================
# This shows the axial structure: fuel elements, grid plates, vessel walls,
# water tanks (upper and lower), and the control cap tank below.

print("Generating XZ cross-section plot (side view)...")

xz_plot = openmc.Plot()
xz_plot.id = 103
xz_plot.filename = "plot_xz_side"
xz_plot.basis = "xz"                    # XZ plane (side view)
xz_plot.origin = (0.0, 0.0, 0.0)       # Centered at midplane
xz_plot.width = (70.0, 70.0)           # Wide enough to see all tanks
xz_plot.pixels = (800, 800)             # Resolution
xz_plot.color_by = "material"
xz_plot.colors = color_by_id


# =============================================================================
# Generate plots using OpenMC
# =============================================================================
plots = openmc.Plots([xy_plot, xy_plot_zoom, xz_plot])
plots.export_to_xml()

# Run the OpenMC plotter
openmc.plot_geometry()

# =============================================================================
# Convert to publication-quality figures using matplotlib
# =============================================================================
# OpenMC generates plot_1.png, plot_2.png, plot_3.png.
# We read these PNGs and add axes labels, titles, and legends.

from PIL import Image

plot_configs = [
    {
        "file_in": "plot_xy_full.png",
        "file_out": "snap10a_xy_full.png",
        "extent": [-32.5, 32.5, -32.5, 32.5],
        "xlabel": "X (cm)",
        "ylabel": "Y (cm)",
        "title": ("SNAP-10A/2 Case 8490a - XY Cross-Section at Midplane (z=0)\n"
                  "28 fuel elements, water-flooded core, full water reflection"),
    },
    {
        "file_in": "plot_2.png",
        "file_out": "snap10a_xy_core.png",
        "extent": [-14.0, 14.0, -14.0, 14.0],
        "xlabel": "X (cm)",
        "ylabel": "Y (cm)",
        "title": ("SNAP-10A/2 Case 8490a - Core Detail (XY at z=0)\n"
                  "37-position triangular pitch array, 28 fuel elements loaded"),
    },
    {
        "file_in": "plot_3.png",
        "file_out": "snap10a_xz_side.png",
        "extent": [-35.0, 35.0, -35.0, 35.0],
        "xlabel": "X (cm)",
        "ylabel": "Z (cm)",
        "title": ("SNAP-10A/2 Case 8490a - XZ Cross-Section (side view)\n"
                  "Core vessel, water tanks, control cap, grid plates"),
    },
]

# Legend patches for material colors
legend_patches = [mpatches.Patch(color=color, label=name)
                  for name, color in material_colors.items()]

for cfg in plot_configs:
    fig, ax = plt.subplots(figsize=(10, 10))
    img = np.array(Image.open(cfg["file_in"]))
    ax.imshow(img, extent=cfg["extent"])
    ax.set_xlabel(cfg["xlabel"], fontsize=12)
    ax.set_ylabel(cfg["ylabel"], fontsize=12)
    ax.set_title(cfg["title"], fontsize=11)
    ax.legend(handles=legend_patches, loc="upper right", fontsize=8,
              framealpha=0.9)
    fig.tight_layout()
    fig.savefig(cfg["file_out"], dpi=150)
    print(f"  Saved: {cfg['file_out']}")
    plt.close(fig)

print("\nAll visualization plots generated successfully.")
