#!/usr/bin/env python3
"""
Visualization script for the IAEA 10 MW MTR Research Reactor Benchmark.

Generates:
1. Full core XY cross-section showing all fuel elements and reflectors
2. Zoomed view of a single fuel element showing individual plate structure
3. Material color legend

Prerequisites:
    Run model.py first to generate the model XML files.

Usage:
    python visualize.py
"""

import openmc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# =============================================================================
# Load the model
# =============================================================================

print("Loading model...")
model = openmc.Model.from_model_xml()
geometry = model.geometry

# =============================================================================
# Material color map
# =============================================================================
# Define consistent colors for each material across all plots.
# Colors chosen for clear visual distinction of the plate structure.

# Find materials by name
mat_colors = {}
for mat in model.materials:
    if "fuel" in mat.name.lower():
        mat_colors[mat] = "red"          # Fuel meat: red (high visibility)
    elif "clad" in mat.name.lower() or "aluminum" in mat.name.lower():
        mat_colors[mat] = "gray"         # Aluminum cladding: gray
    elif "water" in mat.name.lower():
        mat_colors[mat] = "cornflowerblue"  # Water: blue
    elif "beryllium" in mat.name.lower():
        mat_colors[mat] = "gold"         # Be reflector: gold/yellow
    elif "graphite" in mat.name.lower():
        mat_colors[mat] = "dimgray"      # Graphite: dark gray


# =============================================================================
# Plot 1: Full core XY cross-section
# =============================================================================

print("Generating full core cross-section plot...")

# Core dimensions (from model.py)
ELEMENT_WIDTH_X = 7.7
ELEMENT_WIDTH_Y = 8.1
BE_THICKNESS = 7.7
GR_THICKNESS = 10.0

core_half_x = 5 * ELEMENT_WIDTH_X / 2.0
core_half_y = 6 * ELEMENT_WIDTH_Y / 2.0
total_half_x = core_half_x + BE_THICKNESS + GR_THICKNESS
total_half_y = core_half_y + BE_THICKNESS + GR_THICKNESS

fig1, ax1 = plt.subplots(figsize=(12, 14))

# Plot at Z=0 (midplane)
plot = openmc.Plot()
plot.origin = (0.0, 0.0, 0.0)
plot.width = (2 * total_half_x + 2, 2 * total_half_y + 2)  # Slight margin
plot.pixels = (800, 960)
plot.color_by = "material"
plot.colors = mat_colors
plot.basis = "xy"

# Use the openmc plotter
plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()

# Load and display the image
img = openmc.Plot.from_geometry(geometry)
img.origin = (0.0, 0.0, 0.0)
img.width = (2 * total_half_x + 2, 2 * total_half_y + 2)
img.pixels = (800, 960)
img.color_by = "material"
img.colors = mat_colors
img.basis = "xy"

# Use matplotlib to load the generated PNG
import glob
import os

# Generate plot via OpenMC's plotting
plot.filename = "full_core"
plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()

# Load the generated image
if os.path.exists("full_core.png"):
    img_data = plt.imread("full_core.png")
elif os.path.exists("plot_1.png"):
    img_data = plt.imread("plot_1.png")
else:
    # Find any generated PNG
    pngs = glob.glob("*.png")
    if pngs:
        img_data = plt.imread(pngs[0])
    else:
        print("Warning: No plot image found. Skipping full core plot.")
        img_data = None

if img_data is not None:
    extent = [-total_half_x - 1, total_half_x + 1,
              -total_half_y - 1, total_half_y + 1]
    ax1.imshow(img_data, extent=extent, origin="lower")
    ax1.set_xlabel("X (cm)", fontsize=12)
    ax1.set_ylabel("Y (cm)", fontsize=12)
    ax1.set_title("IAEA 10 MW MTR Benchmark - Full Core XY Cross Section\n"
                   "(5x6 lattice with Be + Graphite reflectors)", fontsize=13)

    # Add legend
    legend_patches = [
        mpatches.Patch(color="red", label="Fuel meat (UAl, HEU 93%)"),
        mpatches.Patch(color="gray", label="Aluminum cladding"),
        mpatches.Patch(color="cornflowerblue", label="Water (moderator/coolant)"),
        mpatches.Patch(color="gold", label="Beryllium reflector"),
        mpatches.Patch(color="dimgray", label="Graphite reflector"),
    ]
    ax1.legend(handles=legend_patches, loc="upper right", fontsize=10)

    # Draw grid lines showing element boundaries
    for i in range(6):  # 5 columns -> 6 boundaries
        x = -core_half_x + i * ELEMENT_WIDTH_X
        ax1.axvline(x=x, color="white", linewidth=0.5, alpha=0.3, linestyle="--")
    for j in range(7):  # 6 rows -> 7 boundaries
        y = -core_half_y + j * ELEMENT_WIDTH_Y
        ax1.axhline(y=y, color="white", linewidth=0.5, alpha=0.3, linestyle="--")

    fig1.tight_layout()
    fig1.savefig("full_core_xy.png", dpi=150, bbox_inches="tight")
    print("  Saved: full_core_xy.png")
    plt.close(fig1)


# =============================================================================
# Plot 2: Zoomed view of a single fuel element
# =============================================================================

print("Generating zoomed fuel element plot...")

fig2, ax2 = plt.subplots(figsize=(10, 10))

# Zoom into the center element (approximately at origin)
# The center of the core lattice is at (0, 0)
# A central SFE is at position (0, 0) in the lattice
zoom_width_x = ELEMENT_WIDTH_X + 1.0  # Slight margin around one element
zoom_width_y = ELEMENT_WIDTH_Y + 1.0

plot2 = openmc.Plot()
plot2.filename = "single_element"
plot2.origin = (0.0, 0.0, 0.0)  # Center of the core -> center SFE
plot2.width = (zoom_width_x, zoom_width_y)
plot2.pixels = (1200, 1200)  # High resolution to resolve plates
plot2.color_by = "material"
plot2.colors = mat_colors
plot2.basis = "xy"

plots2 = openmc.Plots([plot2])
plots2.export_to_xml()
openmc.plot_geometry()

# Load the zoomed image
if os.path.exists("single_element.png"):
    img_data2 = plt.imread("single_element.png")
elif os.path.exists("plot_1.png"):
    img_data2 = plt.imread("plot_1.png")
else:
    pngs = glob.glob("*.png")
    img_data2 = plt.imread(pngs[-1]) if pngs else None

if img_data2 is not None:
    extent2 = [-zoom_width_x/2, zoom_width_x/2,
               -zoom_width_y/2, zoom_width_y/2]
    ax2.imshow(img_data2, extent=extent2, origin="lower")
    ax2.set_xlabel("X (cm)", fontsize=12)
    ax2.set_ylabel("Y (cm)", fontsize=12)
    ax2.set_title("IAEA MTR - Single Standard Fuel Element (23 plates)\n"
                   "Plate structure: fuel meat (red) + Al clad (gray) + water channels (blue)",
                   fontsize=12)

    # Add dimension annotations
    # Plate pitch
    ax2.annotate("", xy=(0.175, -zoom_width_y/2 + 0.3),
                 xytext=(-0.175, -zoom_width_y/2 + 0.3),
                 arrowprops=dict(arrowstyle="<->", color="white", lw=1.5))
    ax2.text(0.0, -zoom_width_y/2 + 0.5, "0.35 cm pitch",
             ha="center", va="bottom", color="white", fontsize=8,
             bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.7))

    # Add legend
    legend_patches2 = [
        mpatches.Patch(color="red", label="Fuel meat (0.051 cm)"),
        mpatches.Patch(color="gray", label="Al clad (0.038 cm/side)"),
        mpatches.Patch(color="cornflowerblue", label="Water channel (0.223 cm)"),
    ]
    ax2.legend(handles=legend_patches2, loc="upper right", fontsize=10)

    fig2.tight_layout()
    fig2.savefig("single_element_xy.png", dpi=150, bbox_inches="tight")
    print("  Saved: single_element_xy.png")
    plt.close(fig2)


print("\nVisualization complete!")
