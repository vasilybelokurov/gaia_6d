#!/usr/bin/env python3
"""
Plot histogram of Galactocentric azimuthal angle phi.

Convention: phi = 0 toward the Sun, increasing in prograde direction.

Author: VB
Date: 2025-10-18
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os


def main():
    # Load catalog
    cat_file = os.path.expanduser('~/data/gaia/gaia_dr3_6d_galactocentric.fits')
    print(f"Loading {cat_file}...")
    cat = Table.read(cat_file)
    print(f"Loaded {len(cat):,} stars")

    # Check if phi column exists
    if 'phi' not in cat.colnames:
        print("phi column not found, computing from X, Y...")
        phi = np.arctan2(cat['Y'], -cat['X'])
    else:
        phi = cat['phi']
        print("Using phi column from catalog")

    # Create histogram
    fig, ax = plt.subplots(figsize=(10, 6))

    # Histogram in radians (5-degree bins)
    bins = np.linspace(-np.pi, np.pi, 73)
    counts, edges, _ = ax.hist(phi, bins=bins, color='steelblue',
                                edgecolor='none', alpha=0.8)

    # Add reference lines
    ax.axvline(0, color='red', ls='--', lw=1.5, label=r'$\phi = 0$ (toward Sun)')
    ax.axvline(np.pi/2, color='orange', ls=':', lw=1.5, label=r'$\phi = \pm\pi/2$')
    ax.axvline(-np.pi/2, color='orange', ls=':', lw=1.5)

    # Labels
    ax.set_xlabel(r'Azimuthal angle $\phi$ [rad]', fontsize=12)
    ax.set_ylabel('Number of stars', fontsize=12)
    ax.set_yscale('log')
    ax.set_title(f'Distribution of Galactocentric azimuthal angle\n'
                 f'({len(cat):,} stars, $\phi=0$ toward Sun)', fontsize=13)

    # Secondary x-axis in degrees
    ax2 = ax.twiny()
    ax2.set_xlim(-180, 180)
    ax2.set_xlabel('Azimuthal angle [deg]', fontsize=11)

    ax.set_xlim(-np.pi, np.pi)
    ax.legend(loc='upper right')

    # Print statistics
    print(f"\nStatistics:")
    print(f"  phi range: {np.min(phi):.3f} to {np.max(phi):.3f} rad")
    print(f"  phi range: {np.degrees(np.min(phi)):.1f} to {np.degrees(np.max(phi)):.1f} deg")
    print(f"  median phi: {np.median(phi):.3f} rad ({np.degrees(np.median(phi)):.1f} deg)")

    # Save
    outfile = 'plots/phi_histogram.png'
    os.makedirs('plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, facecolor='white')
    print(f"\nSaved: {outfile}")
    plt.close()


if __name__ == '__main__':
    main()
