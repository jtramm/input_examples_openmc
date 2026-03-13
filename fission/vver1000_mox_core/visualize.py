#!/usr/bin/env python3
"""
Generate 2D slice plots of the VVER-1000 MOX full-core geometry.

This script creates an XY cross-section plot of the full core showing the
MOX/UOX loading pattern. The different fuel types are colored distinctly
so the checkerboard-like MOX loading pattern is clearly visible.

The VVER-1000 core has 163 hexagonal fuel assemblies in a hexagonal
arrangement. About 30% of the assemblies are loaded with MOX fuel
(weapons-grade plutonium in natural UO2), while the remainder are UOX
assemblies of varying enrichment.

The plot shows:
  - MOX fuel pins in distinct colors by zone (2.4/2.7/3.6% PuO2)
  - Fresh UOX (3.7%) pins at the core periphery
  - Burned UOX (2.0%) pins in the interior
  - Water moderator (blue), cladding (grey), and steel baffle

Prerequisites:
    Run model.py first to generate the XML input files.

Usage:
    python visualize.py
"""

import openmc
from model import build_model


def make_plots():
    """Create and export geometry visualization plots for the full core."""

    model = build_model()

    # Get material references by name for color assignment
    materials = {m.name: m for m in model.materials}

    # Color scheme: distinguish fuel types clearly
    # UOX fuels in shades of red/orange, MOX in shades of green/purple
    mat_colors = {
        materials['UO2 3.7%']:       (220, 60, 60),     # bright red - fresh UOX
        materials['UO2 2.0%']:       (255, 140, 60),    # orange - burned UOX
        materials['MOX 2.4% PuO2']:  (60, 180, 60),     # green - MOX inner
        materials['MOX 2.7% PuO2']:  (40, 140, 40),     # darker green - MOX mid
        materials['MOX 3.6% PuO2']:  (20, 100, 20),     # dark green - MOX outer
        materials['E110 Cladding']:  (180, 180, 180),   # grey - cladding
        materials['Borated Water']:  (100, 150, 255),   # blue - water
        materials['Steel Baffle']:   (100, 80, 60),     # brown - steel
    }

    # --- Plot 1: Full core XY cross-section ---
    # This shows all 163 assemblies with the MOX/UOX loading pattern.
    # The core diameter is ~316 cm, plus the steel baffle and water reflector.
    # Total plot width covers out to the reflector outer boundary (~210 cm radius).
    full_core = openmc.Plot(name='core_full_xy')
    full_core.filename = 'plot_core_full'
    full_core.origin = (0.0, 0.0, 0.0)
    full_core.width = (440.0, 440.0)
    full_core.pixels = (4000, 4000)
    full_core.color_by = 'material'
    full_core.colors = mat_colors

    # --- Plot 2: Zoomed view of central assemblies ---
    # Shows the inner few rings of assemblies with pin-level detail visible.
    # This view clearly shows the MOX three-zone pin layout (inner/mid/outer
    # Pu enrichment) versus the uniform UOX assemblies.
    zoom_center = openmc.Plot(name='core_zoom_center')
    zoom_center.filename = 'plot_core_zoom_center'
    zoom_center.origin = (0.0, 0.0, 0.0)
    zoom_center.width = (80.0, 80.0)
    zoom_center.pixels = (3000, 3000)
    zoom_center.color_by = 'material'
    zoom_center.colors = mat_colors

    # --- Plot 3: Single assembly detail ---
    # Zoomed to show individual pin structure within one MOX assembly.
    # The three concentric zones of different Pu content are visible:
    # inner (lighter green), middle, and outer (darker green).
    zoom_assy = openmc.Plot(name='core_zoom_assembly')
    zoom_assy.filename = 'plot_core_zoom_assembly'
    zoom_assy.origin = (0.0, 23.6, 0.0)  # offset to a MOX assembly position
    zoom_assy.width = (26.0, 26.0)
    zoom_assy.pixels = (2000, 2000)
    zoom_assy.color_by = 'material'
    zoom_assy.colors = mat_colors

    plots = openmc.Plots([full_core, zoom_center, zoom_assy])
    plots.export_to_xml()

    print("Plot XML exported. Run 'openmc --plot' to generate images.")
    print("Plots defined:")
    print(f"  1. {full_core.filename} - Full core ({full_core.pixels[0]}x"
          f"{full_core.pixels[1]} px, {full_core.width[0]} cm wide)")
    print(f"  2. {zoom_center.filename} - Central region ({zoom_center.pixels[0]}x"
          f"{zoom_center.pixels[1]} px)")
    print(f"  3. {zoom_assy.filename} - Single assembly detail ({zoom_assy.pixels[0]}x"
          f"{zoom_assy.pixels[1]} px)")


if __name__ == '__main__':
    make_plots()
