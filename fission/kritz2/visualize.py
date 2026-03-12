#!/usr/bin/env python3
"""
Visualize the KRITZ-2 pin cell geometry.

Produces an XY cross-section plot showing the fuel pellet, cladding,
and moderator regions.  The plot is saved as 'kritz2_pincell.png'.

Usage:
    python visualize.py [--config 2:1-cold]
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Import model builder to reuse the configuration data
from model import CONFIGS


def visualize(config_name):
    """
    Draw the XY cross section of the KRITZ-2 pin cell.

    The plot shows three concentric regions:
      - Blue circle: fuel pellet (UO2 or MOX)
      - Grey annulus: Zircaloy-2 cladding
      - Cyan square: moderator (borated light water)

    Dimensions are annotated on the plot.
    """
    cfg = CONFIGS[config_name]

    fuel_r = cfg['fuel_radius']
    clad_r = cfg['clad_outer_radius']
    half_pitch = cfg['pitch'] / 2.0

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    # Draw the moderator region (square background)
    moderator_patch = patches.Rectangle(
        (-half_pitch, -half_pitch), cfg['pitch'], cfg['pitch'],
        linewidth=2, edgecolor='black', facecolor='#B3E5FC',
        label=f'Moderator (H2O + {cfg["mod_boron_ppm"]} ppm B)'
    )
    ax.add_patch(moderator_patch)

    # Draw the cladding (grey annulus)
    clad_patch = patches.Circle(
        (0, 0), clad_r,
        linewidth=1.5, edgecolor='black', facecolor='#9E9E9E',
        label='Zircaloy-2 Cladding'
    )
    ax.add_patch(clad_patch)

    # Draw the fuel pellet
    fuel_color = '#E53935' if cfg['fuel_type'] == 'MOX' else '#FFA726'
    fuel_patch = patches.Circle(
        (0, 0), fuel_r,
        linewidth=1.5, edgecolor='black', facecolor=fuel_color,
        label=f'Fuel ({cfg["fuel_type"]})'
    )
    ax.add_patch(fuel_patch)

    # Annotate dimensions
    # Fuel radius
    ax.annotate('', xy=(fuel_r, 0), xytext=(0, 0),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(fuel_r / 2, 0.02, f'r_f = {fuel_r:.4f} cm',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Clad outer radius
    ax.annotate('', xy=(0, clad_r), xytext=(0, 0),
                arrowprops=dict(arrowstyle='<->', color='navy', lw=1.5))
    ax.text(0.02, clad_r / 2, f'r_c = {clad_r:.4f} cm',
            ha='left', va='center', fontsize=9, fontweight='bold', color='navy')

    # Pitch
    ax.annotate('', xy=(half_pitch, -half_pitch),
                xytext=(-half_pitch, -half_pitch),
                arrowprops=dict(arrowstyle='<->', color='green', lw=1.5))
    ax.text(0, -half_pitch + 0.03, f'pitch = {cfg["pitch"]:.4f} cm',
            ha='center', va='bottom', fontsize=10, fontweight='bold',
            color='green')

    # Set axis properties
    margin = 0.1
    ax.set_xlim(-half_pitch - margin, half_pitch + margin)
    ax.set_ylim(-half_pitch - margin, half_pitch + margin)
    ax.set_aspect('equal')
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'KRITZ-2 Pin Cell Cross Section\n{cfg["description"]}',
                 fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outfile = 'kritz2_pincell.png'
    plt.savefig(outfile, dpi=150)
    print(f"Saved plot to {outfile}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize KRITZ-2 pin cell geometry'
    )
    parser.add_argument(
        '--config', type=str, default='2:1-cold',
        choices=['2:1-cold', '2:1-hot', '2:13-cold', '2:13-hot',
                 '2:19-cold', '2:19-hot'],
        help='Configuration to visualize (default: 2:1-cold)'
    )
    args = parser.parse_args()

    visualize(args.config)


if __name__ == '__main__':
    main()
