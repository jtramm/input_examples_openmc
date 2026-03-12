#!/usr/bin/env python3
"""
Generate 2D slice plots of the HTR-10 pebble bed reactor geometry.

This script creates visualization plots for both model modes:

  Pin mode (single pebble with explicit TRISO):
    - XY cross-section showing TRISO particles in the fuel zone
    - Zoomed-in view of individual TRISO particles
    - XZ cross-section showing pebble structure

  Core mode (cylindrical core with homogenized pebbles):
    - XY cross-section at core midplane
    - XZ cross-section showing core + reflector

Prerequisites:
    Run model.py first to generate the XML input files.

Usage:
    python visualize.py [--model {pin,core}]
"""

import argparse

import openmc
from model import build_model, PEBBLE_OUTER_RADIUS, FUEL_ZONE_RADIUS, \
    CORE_RADIUS, CORE_HEIGHT, REFLECTOR_THICKNESS, OPYC_OUTER_R


def make_plots(model_type='pin'):
    """Create and export geometry visualization plots."""

    model = build_model(model_type=model_type)

    # Get material references from the model
    materials = {m.name: m for m in model.materials}

    if model_type == 'pin':
        _plot_pin(model, materials)
    else:
        _plot_core(model, materials)

    print("Generated plot images.")


def _plot_pin(model, materials):
    """Create plots for the single-pebble (pin) model."""

    # Color map for TRISO and pebble materials
    colors = {}
    if 'UO2 Kernel' in materials:
        colors[materials['UO2 Kernel']] = (220, 30, 30)       # red
    if 'Buffer PyC' in materials:
        colors[materials['Buffer PyC']] = (255, 200, 100)      # tan
    if 'Inner PyC' in materials:
        colors[materials['Inner PyC']] = (100, 100, 100)       # dark grey
    if 'SiC' in materials:
        colors[materials['SiC']] = (50, 180, 50)               # green
    if 'Outer PyC' in materials:
        colors[materials['Outer PyC']] = (150, 150, 150)       # grey
    if 'Graphite Matrix' in materials:
        colors[materials['Graphite Matrix']] = (80, 80, 80)    # charcoal
    if 'Graphite Shell' in materials:
        colors[materials['Graphite Shell']] = (60, 60, 60)     # darker
    if 'Helium Coolant' in materials:
        colors[materials['Helium Coolant']] = (200, 220, 255)  # light blue

    # --- Plot 1: Full pebble XY cross-section ---
    # Shows the fuel zone with TRISO particles and the graphite shell
    pebble_plot = openmc.Plot(name='pebble_xy')
    pebble_plot.filename = 'plot_pebble_xy'
    pebble_plot.origin = (0.0, 0.0, 0.0)
    width = 2.2 * PEBBLE_OUTER_RADIUS
    pebble_plot.width = (width, width)
    pebble_plot.pixels = (800, 800)
    pebble_plot.color_by = 'material'
    pebble_plot.colors = colors

    # --- Plot 2: Zoomed in to see individual TRISO particles ---
    # At this zoom level, individual TRISO layers should be visible
    zoom_plot = openmc.Plot(name='triso_zoom')
    zoom_plot.filename = 'plot_triso_zoom'
    zoom_plot.origin = (0.8, 0.8, 0.0)  # off-center to catch some TRISOs
    zoom_width = 20 * OPYC_OUTER_R  # ~0.9 cm window
    zoom_plot.width = (zoom_width, zoom_width)
    zoom_plot.pixels = (800, 800)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = colors

    # --- Plot 3: XZ cross-section showing pebble structure ---
    xz_plot = openmc.Plot(name='pebble_xz')
    xz_plot.filename = 'plot_pebble_xz'
    xz_plot.basis = 'xz'
    xz_plot.origin = (0.0, 0.0, 0.0)
    xz_plot.width = (width, width)
    xz_plot.pixels = (800, 800)
    xz_plot.color_by = 'material'
    xz_plot.colors = colors

    plots = openmc.Plots([pebble_plot, zoom_plot, xz_plot])
    model.plots = plots
    model.export_to_model_xml()
    openmc.plot_geometry()


def _plot_core(model, materials):
    """Create plots for the cylindrical core model."""

    colors = {}
    if 'Homogenized Fuel Zone' in materials:
        colors[materials['Homogenized Fuel Zone']] = (220, 30, 30)
    if 'Graphite Shell' in materials:
        colors[materials['Graphite Shell']] = (100, 100, 100)
    if 'Graphite Reflector' in materials:
        colors[materials['Graphite Reflector']] = (60, 60, 60)
    if 'Helium Coolant' in materials:
        colors[materials['Helium Coolant']] = (200, 220, 255)

    # --- Plot 1: XY cross-section at core midplane ---
    xy_plot = openmc.Plot(name='core_xy')
    xy_plot.filename = 'plot_core_xy'
    xy_plot.origin = (0.0, 0.0, CORE_HEIGHT / 2.0)
    core_width = 2.2 * (CORE_RADIUS + REFLECTOR_THICKNESS)
    xy_plot.width = (core_width, core_width)
    xy_plot.pixels = (800, 800)
    xy_plot.color_by = 'material'
    xy_plot.colors = colors

    # --- Plot 2: XZ cross-section ---
    xz_plot = openmc.Plot(name='core_xz')
    xz_plot.filename = 'plot_core_xz'
    xz_plot.basis = 'xz'
    xz_plot.origin = (0.0, 0.0, CORE_HEIGHT / 2.0)
    z_width = 2.2 * (CORE_HEIGHT / 2.0 + REFLECTOR_THICKNESS)
    xz_plot.width = (core_width, z_width)
    xz_plot.pixels = (800, 600)
    xz_plot.color_by = 'material'
    xz_plot.colors = colors

    # --- Plot 3: Zoomed view showing individual pebbles ---
    zoom_plot = openmc.Plot(name='pebbles_zoom')
    zoom_plot.filename = 'plot_pebbles_zoom'
    zoom_plot.origin = (0.0, 0.0, CORE_HEIGHT / 2.0)
    zoom_plot.width = (40.0, 40.0)
    zoom_plot.pixels = (800, 800)
    zoom_plot.color_by = 'material'
    zoom_plot.colors = colors

    plots = openmc.Plots([xy_plot, xz_plot, zoom_plot])
    model.plots = plots
    model.export_to_model_xml()
    openmc.plot_geometry()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize HTR-10 model geometry'
    )
    parser.add_argument('--model', type=str, default='pin',
                        choices=['pin', 'core'],
                        help='Model type to visualize')
    args = parser.parse_args()
    make_plots(args.model)


if __name__ == '__main__':
    main()
