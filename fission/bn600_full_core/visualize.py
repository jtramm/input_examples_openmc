#!/usr/bin/env python3
"""
Generate visualization plots of the BN-600 full-core fast reactor geometry.

This script creates XY (radial) and XZ (axial) slice plots showing:
  1. Full-core XY cross-section at the midplane (z=0), colored by material
     to reveal the three enrichment zones and radial blanket
  2. Zoomed XY view of individual assemblies showing pin arrangement
  3. XZ axial cross-section through the core center, showing the active
     fuel region and upper/lower axial blankets

The distinct enrichment zones (LEZ/MEZ/HEZ) should be clearly visible as
concentric hexagonal rings with different fuel colors, surrounded by the
radial blanket ring.

Prerequisites:
    Run model.py first to generate the XML input files.

Usage:
    python visualize.py
"""

import openmc
from model import build_model


def make_plots():
    """Create and export geometry visualization plots."""

    model = build_model()

    # Get material references for custom color mapping
    materials = {m.name: m for m in model.materials}
    fuel_lez = materials['UO2 LEZ 17%']
    fuel_mez = materials['UO2 MEZ 21%']
    fuel_hez = materials['UO2 HEZ 26%']
    blanket = materials['UO2 Blanket 0.3%']
    steel = materials['SS Cladding/Wrapper']
    sodium = materials['Sodium Coolant']

    # Color scheme: enrichment zones in warm colors (lighter = more enriched),
    # blanket in green, sodium in blue, steel in grey
    mat_colors = {
        fuel_lez: (180, 60, 60),      # dark red - LEZ (17% enrichment)
        fuel_mez: (220, 120, 40),     # orange - MEZ (21% enrichment)
        fuel_hez: (240, 200, 50),     # yellow - HEZ (26% enrichment)
        blanket:  (60, 160, 60),      # green - depleted UO2 blanket
        steel:    (160, 160, 160),    # grey - stainless steel
        sodium:   (80, 130, 220),     # blue - sodium coolant
    }

    # Core dimensions for plot sizing
    assy_pitch = 9.82
    n_rings = 13
    core_extent = n_rings * assy_pitch  # ~128 cm

    # --- Plot 1: Full-core XY at midplane ---
    # This is the key plot showing the BN-600 core layout.
    # The three enrichment zones appear as concentric hexagonal rings
    # of assemblies, each with slightly different fuel color.
    # The radial blanket (green) surrounds the fuel zones.
    full_xy = openmc.Plot(name='core_xy_midplane')
    full_xy.filename = 'plot_core_xy'
    full_xy.basis = 'xy'
    full_xy.origin = (0.0, 0.0, 0.0)
    full_xy.width = (core_extent * 1.1, core_extent * 1.1)
    full_xy.pixels = (2000, 2000)
    full_xy.color_by = 'material'
    full_xy.colors = mat_colors

    # --- Plot 2: Zoomed XY showing assembly-level detail ---
    # This shows individual fuel pins within the hex assemblies.
    # The hex wrapper (grey steel) around each assembly and the
    # sodium inter-assembly gap (blue) should be visible.
    zoom_xy = openmc.Plot(name='assembly_detail_xy')
    zoom_xy.filename = 'plot_assembly_detail'
    zoom_xy.basis = 'xy'
    zoom_xy.origin = (0.0, 0.0, 0.0)
    zoom_xy.width = (25.0, 25.0)  # ~2.5 assemblies across
    zoom_xy.pixels = (1500, 1500)
    zoom_xy.color_by = 'material'
    zoom_xy.colors = mat_colors

    # --- Plot 3: XZ axial cross-section through core center ---
    # This slice shows the vertical structure of the core:
    #   - Lower axial blanket (green, 30 cm)
    #   - Active fuel (red/orange/yellow, 100 cm)
    #   - Upper axial blanket (green, 30 cm)
    # The compact nature of a fast reactor core (100 cm active height
    # vs ~400 cm for an LWR) should be apparent.
    xz_slice = openmc.Plot(name='core_xz_axial')
    xz_slice.filename = 'plot_core_xz'
    xz_slice.basis = 'xz'
    xz_slice.origin = (0.0, 0.0, 0.0)
    xz_slice.width = (core_extent * 1.1, 200.0)
    xz_slice.pixels = (2000, 1400)
    xz_slice.color_by = 'material'
    xz_slice.colors = mat_colors

    # --- Plot 4: Zoomed XZ showing single assembly axial structure ---
    # Close-up of axial blanket transition showing individual pins
    zoom_xz = openmc.Plot(name='assembly_axial_xz')
    zoom_xz.filename = 'plot_assembly_axial'
    zoom_xz.basis = 'xz'
    zoom_xz.origin = (0.0, 0.0, 0.0)
    zoom_xz.width = (12.0, 180.0)
    zoom_xz.pixels = (600, 1400)
    zoom_xz.color_by = 'material'
    zoom_xz.colors = mat_colors

    # Collect and export
    plots = openmc.Plots([full_xy, zoom_xy, xz_slice, zoom_xz])
    plots.export_to_xml()

    print("Plot XML exported. Run 'openmc --plot' to generate PNG images.")
    print("Plots defined:")
    print(f"  1. {full_xy.filename} - Full-core XY cross-section at midplane")
    print(f"     ({full_xy.pixels[0]}x{full_xy.pixels[1]} px, "
          f"{full_xy.width[0]:.0f}x{full_xy.width[1]:.0f} cm)")
    print(f"  2. {zoom_xy.filename} - Zoomed XY showing pin-level detail")
    print(f"     ({zoom_xy.pixels[0]}x{zoom_xy.pixels[1]} px, "
          f"{zoom_xy.width[0]:.0f}x{zoom_xy.width[1]:.0f} cm)")
    print(f"  3. {xz_slice.filename} - Full-core XZ axial cross-section")
    print(f"     ({xz_slice.pixels[0]}x{xz_slice.pixels[1]} px, "
          f"{xz_slice.width[0]:.0f}x{xz_slice.width[1]:.0f} cm)")
    print(f"  4. {zoom_xz.filename} - Single assembly axial structure")
    print(f"     ({zoom_xz.pixels[0]}x{zoom_xz.pixels[1]} px, "
          f"{zoom_xz.width[0]:.0f}x{zoom_xz.width[1]:.0f} cm)")


if __name__ == '__main__':
    make_plots()
