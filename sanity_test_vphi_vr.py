#!/usr/bin/env python3
"""
Sanity test: V_phi vs V_r density plot for quality-cut Gaia 6D sample.

Applies quality cuts:
- distance < 10 kpc
- ruwe < 1.4
- rv_error < 20 km/s

Converts to Galactocentric cylindrical velocities and plots V_phi vs V_r.

Author: VB
Date: 2025-10-18
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os


def cylindrical_velocities(X, Y, VX, VY):
    """
    Convert Galactocentric Cartesian velocities to cylindrical.

    Parameters
    ----------
    X, Y : array [kpc]
        Galactocentric Cartesian positions
    VX, VY : array [km/s]
        Galactocentric Cartesian velocities

    Returns
    -------
    V_R, V_phi : array [km/s]
        Cylindrical radial and azimuthal velocities
    """
    R_cyl = np.sqrt(X**2 + Y**2)

    # Radial velocity (in plane)
    V_R = (X * VX + Y * VY) / R_cyl

    # Azimuthal velocity (rotation)
    V_phi = (-X * VY + Y * VX) / R_cyl

    return V_R, V_phi


def main():
    """Main execution."""

    print("=" * 70)
    print("SANITY TEST: V_phi vs V_R DENSITY PLOT")
    print("=" * 70)

    # ========================================================================
    # Load catalog
    # ========================================================================

    input_file = os.path.expanduser('~/data/gaia/gaia_dr3_6d_galactocentric.fits')
    print(f"\nLoading: {input_file}")

    cat = Table.read(input_file)
    n_total = len(cat)
    print(f"  Total stars: {n_total:,}")

    # ========================================================================
    # Apply quality cuts
    # ========================================================================

    print("\nApplying quality cuts:")
    print("  - distance < 10 kpc")
    print("  - ruwe < 1.4")
    print("  - rv_error < 20 km/s")

    # Distance is in pc, convert to kpc
    mask_dist = (cat['distance'] / 1000.0) < 10.0
    mask_ruwe = cat['ruwe'] < 1.4
    mask_rv_err = cat['rv_error'] < 20.0

    mask_quality = mask_dist & mask_ruwe & mask_rv_err

    cat_clean = cat[mask_quality]
    n_clean = len(cat_clean)

    print(f"\n  Stars passing cuts: {n_clean:,} ({100*n_clean/n_total:.1f}%)")
    print(f"    distance < 10 kpc: {np.sum(mask_dist):,}")
    print(f"    ruwe < 1.4: {np.sum(mask_ruwe):,}")
    print(f"    rv_error < 20 km/s: {np.sum(mask_rv_err):,}")

    # ========================================================================
    # Convert to cylindrical velocities
    # ========================================================================

    print("\nConverting to Galactocentric cylindrical velocities...")

    V_R, V_phi = cylindrical_velocities(
        cat_clean['X'],
        cat_clean['Y'],
        cat_clean['VX'],
        cat_clean['VY']
    )

    print(f"  V_R range: {np.min(V_R):.1f} to {np.max(V_R):.1f} km/s")
    print(f"  V_phi range: {np.min(V_phi):.1f} to {np.max(V_phi):.1f} km/s")

    # ========================================================================
    # Create 2D histogram
    # ========================================================================

    print("\nCreating 2D density histogram...")

    v_min, v_max = -350, 350
    n_bins = 200

    # Mask to velocity range
    mask_vrange = (V_R >= v_min) & (V_R <= v_max) & (V_phi >= v_min) & (V_phi <= v_max)
    V_R_plot = V_R[mask_vrange]
    V_phi_plot = V_phi[mask_vrange]

    print(f"  Stars in plot range: {len(V_R_plot):,} ({100*len(V_R_plot)/n_clean:.1f}%)")

    # Create 2D histogram
    H, xedges, yedges = np.histogram2d(
        V_R_plot,
        V_phi_plot,
        bins=n_bins,
        range=[[v_min, v_max], [v_min, v_max]]
    )

    # ========================================================================
    # Plot
    # ========================================================================

    print("\nCreating plot...")

    fig, ax = plt.subplots(figsize=(10, 9))

    # Plot 2D histogram with log scale
    H_max = H.max() if H.max() > 0 else 1
    H_min = max(1, H[H > 0].min()) if np.any(H > 0) else 1

    im = ax.imshow(
        H.T,
        origin='lower',
        extent=[v_min, v_max, v_min, v_max],
        aspect='auto',
        cmap='viridis',
        norm=plt.matplotlib.colors.LogNorm(vmin=H_min, vmax=H_max),
        interpolation='nearest'
    )

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Number of stars per bin', fontsize=12)

    # Axes
    ax.set_xlabel('$V_R$ (km/s)', fontsize=14)
    ax.set_ylabel('$V_\\phi$ (km/s)', fontsize=14)
    ax.set_xlim(v_min, v_max)
    ax.set_ylim(v_min, v_max)

    # Grid
    ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Zero lines
    ax.axhline(0, color='white', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(0, color='white', linestyle='--', linewidth=1, alpha=0.5)

    # Title with sample info
    title = (
        f'Gaia DR3 6D: $V_\\phi$ vs $V_R$ density\n'
        f'$d < 10$ kpc, RUWE $< 1.4$, $\\sigma_{{RV}} < 20$ km/s '
        f'({n_clean:,} stars)'
    )
    ax.set_title(title, fontsize=13, pad=10)

    plt.tight_layout()

    # ========================================================================
    # Save plot
    # ========================================================================

    output_dir = 'plots'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'vphi_vr_density.png')

    plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_file}")

    plt.close()

    # ========================================================================
    # Print statistics
    # ========================================================================

    print("\n" + "=" * 70)
    print("SAMPLE STATISTICS")
    print("=" * 70)
    print(f"Total stars in catalog: {n_total:,}")
    print(f"Stars passing quality cuts: {n_clean:,} ({100*n_clean/n_total:.1f}%)")
    print(f"Stars in plot range [{v_min}, {v_max}] km/s: {len(V_R_plot):,}")
    print(f"\nMedian V_R: {np.median(V_R_plot):.1f} km/s")
    print(f"Median V_phi: {np.median(V_phi_plot):.1f} km/s")
    print(f"Std V_R: {np.std(V_R_plot):.1f} km/s")
    print(f"Std V_phi: {np.std(V_phi_plot):.1f} km/s")
    print("=" * 70)


if __name__ == '__main__':
    main()
