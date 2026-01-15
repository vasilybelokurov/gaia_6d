#!/usr/bin/env python3
"""
Figure 3 phase-space folds split by Galactocentric azimuthal angle phi.

Produces chevron plots for 9 bins of 40 degrees starting at phi = -20 deg.
Based on figure3_phase_space_folds.ipynb.

Convention: phi = 0 toward the Sun, increasing in prograde direction.

Author: VB
Date: 2025-01-15
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
from scipy.ndimage import gaussian_filter
import os


def column_normalize(H):
    """Normalize each column (r bin) to sum to 1."""
    H_norm = H.copy()
    for i in range(H.shape[1]):
        col_sum = H[:, i].sum()
        if col_sum > 0:
            H_norm[:, i] /= col_sum
    return H_norm


def make_chevron_plot(r, v_r, mask_pro, mask_ret, mask_gse, use,
                      phi_min, phi_max, n_stars_info, outfile):
    """
    Create the 2x3 chevron plot for a given phi bin.

    Parameters
    ----------
    r : array
        Galactocentric spherical radius (kpc)
    v_r : array
        Galactocentric radial velocity (km/s)
    mask_pro : array
        Prograde (L_z > 0) mask
    mask_ret : array
        Retrograde (L_z < 0) mask
    mask_gse : array
        GS/E selection mask (|L_z| < 700)
    use : array
        Quality cuts mask
    phi_min, phi_max : float
        Phi bin edges in degrees
    n_stars_info : dict
        Star counts for title
    outfile : str
        Output filename
    """
    # Axis ranges
    r_range = (2.0, 14.0)
    vr_range = (-500, 500)

    # Pixels
    n_r_top, n_vr_top = 75, 150
    n_r_bot, n_vr_bot = 50, 100

    # FWHM for background subtraction
    fwhm = {'prograde': 9, 'retrograde': 15, 'combined': 12}
    fwhm_smooth = 1.3

    # Samples
    samples = {
        'prograde': mask_pro,
        'retrograde': mask_ret,
        'combined': mask_gse
    }

    # Create histograms
    H_top = {}
    H_bot = {}

    for name, mask in samples.items():
        mask_final = mask & use

        H_top[name], _, _ = np.histogram2d(
            v_r[mask_final], r[mask_final],
            bins=[n_vr_top, n_r_top],
            range=[vr_range, r_range]
        )

        H_bot[name], _, _ = np.histogram2d(
            v_r[mask_final], r[mask_final],
            bins=[n_vr_bot, n_r_bot],
            range=[vr_range, r_range]
        )

    # Column normalize
    H_top_norm = {name: column_normalize(H_top[name]) for name in samples}
    H_bot_norm = {name: column_normalize(H_bot[name]) for name in samples}

    # Background subtraction
    H_bot_sub = {}
    for name in samples:
        sigma = fwhm[name] / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        background = gaussian_filter(H_bot_norm[name], sigma=sigma, mode='constant')

        sigma_smooth = fwhm_smooth / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        H_smoothed = gaussian_filter(H_bot_norm[name], sigma=sigma_smooth, mode='constant')

        H_bot_sub[name] = H_smoothed - background

    # Plotting
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    titles = [
        f'$L_z > 0$ (pro, N={n_stars_info["prograde"]:,})',
        f'$L_z < 0$ (ret, N={n_stars_info["retrograde"]:,})',
        f'Combined (N={n_stars_info["combined"]:,})'
    ]
    extent = [r_range[0], r_range[1], vr_range[0], vr_range[1]]

    # Top row: column-normalized
    for i, (name, title) in enumerate(zip(['prograde', 'retrograde', 'combined'], titles)):
        ax = axes[0, i]

        H_plot = np.log10(H_top_norm[name], where=(H_top_norm[name] > 0))
        H_plot[H_top_norm[name] == 0] = np.nan

        finite_vals = H_plot[np.isfinite(H_plot)]
        if len(finite_vals) > 0:
            vmin = np.percentile(finite_vals, 5)
            vmax = np.percentile(finite_vals, 95)
        else:
            vmin, vmax = -3, 0

        im = ax.imshow(H_plot, origin='lower', extent=extent, aspect='auto',
                       cmap='gray_r', vmin=vmin, vmax=vmax, interpolation='nearest')

        ax.set_xlabel('$r$ (kpc)', fontsize=12)
        if i == 0:
            ax.set_ylabel('$v_r$ (km s$^{-1}$)', fontsize=12)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('log(density)', fontsize=10)

    # Bottom row: background-subtracted
    for i, name in enumerate(['prograde', 'retrograde', 'combined']):
        ax = axes[1, i]

        H_plot = H_bot_sub[name]

        vmin = np.percentile(H_plot, 5)
        vmax = np.percentile(H_plot, 95)
        vlim = max(abs(vmin), abs(vmax))
        if vlim == 0:
            vlim = 0.01

        im = ax.imshow(H_plot, origin='lower', extent=extent, aspect='auto',
                       cmap='RdBu_r', vmin=-vlim, vmax=vlim, interpolation='nearest')

        ax.set_xlabel('$r$ (kpc)', fontsize=12)
        if i == 0:
            ax.set_ylabel('$v_r$ (km s$^{-1}$)', fontsize=12)
        ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Residual', fontsize=10)

    # Row labels
    fig.text(0.98, 0.75, 'Column-normalized', rotation=270, va='center', fontsize=14, weight='bold')
    fig.text(0.98, 0.30, 'Background-subtracted', rotation=270, va='center', fontsize=14, weight='bold')

    # Main title with phi range
    if phi_max > phi_min:
        phi_label = f'${phi_min:.0f}^\\circ < \\phi < {phi_max:.0f}^\\circ$'
    else:
        # Wrapping case
        phi_label = f'$\\phi \\in [{phi_min:.0f}^\\circ, {phi_max:.0f}^\\circ]$'

    fig.suptitle(f'GS/E Phase-Space Folds: {phi_label}', fontsize=14, weight='bold', y=1.02)

    plt.tight_layout(rect=[0, 0, 0.96, 1])
    plt.savefig(outfile, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {outfile}")


def main():
    print("=" * 70)
    print("FIGURE 3 PHASE-SPACE FOLDS BY PHI BIN")
    print("=" * 70)

    # Load catalog
    cat_file = os.path.expanduser('~/data/gaia/gaia_dr3_6d_galactocentric.fits')
    print(f"\nLoading {cat_file}...")
    cat = Table.read(cat_file)
    print(f"Total stars: {len(cat):,}")

    # Remove GC members
    vasiliev_file = os.path.expanduser('~/data/catalogues/gc_members_gaia_vasiliev.fits')
    if os.path.exists(vasiliev_file):
        dgc = fits.getdata(vasiliev_file, 1)
        gc_source_ids = dgc['SOURCE_ID'][dgc['prob'] > 0.0]
        mask_gc = np.isin(cat['source_id'], gc_source_ids)
        print(f"Removing {np.sum(mask_gc):,} GC members")
        cat = cat[~mask_gc]
    else:
        print(f"Warning: GC catalog not found: {vasiliev_file}")

    # Compute angular momentum L_z
    L_z = cat['X'] * cat['VY'] - cat['Y'] * cat['VX']
    L_z *= -1  # Sign convention

    # GS/E selection
    L_z_cut = 0.7e3
    mask_gse_base = np.abs(L_z) < L_z_cut
    mask_pro_base = (L_z > 0) & mask_gse_base
    mask_ret_base = (L_z < 0) & mask_gse_base

    # Galactocentric spherical coordinates
    r = np.sqrt(cat['X']**2 + cat['Y']**2 + cat['Z']**2)
    v_r = (cat['X']*cat['VX'] + cat['Y']*cat['VY'] + cat['Z']*cat['VZ']) / r

    # Compute phi (azimuthal angle, phi=0 toward Sun)
    phi_rad = np.arctan2(cat['Y'], -cat['X'])
    phi_deg = np.degrees(phi_rad)

    # Quality cuts
    d_helio = cat['distance'] / 1000.0
    mask_dist = d_helio < 15.0
    mask_dist_err = (cat['distance'] / cat['distance_error']) > 3.0
    parallax_sn = np.abs(cat['parallax'] / cat['parallax_error'])
    mask_parallax = parallax_sn > 10.0
    mask_ebv = cat['ebv'] < 100.0

    use_base = mask_dist & mask_dist_err & mask_parallax & mask_ebv

    print(f"\nGS/E sample (all phi, after cuts): {np.sum(mask_gse_base & use_base):,} stars")

    # Define phi bins: 40 degrees starting at -20
    # Bins: [-20,20], [20,60], [60,100], [100,140], [140,180], [-180,-140], [-140,-100], [-100,-60], [-60,-20]
    phi_edges = [-20, 20, 60, 100, 140, 180, -140, -100, -60, -20]

    # Create output directory
    os.makedirs('plots/phi_bins', exist_ok=True)

    print(f"\nProcessing 9 phi bins of 40 degrees each...")
    print("-" * 70)

    for i in range(9):
        phi_min = phi_edges[i]
        phi_max = phi_edges[i + 1]

        # Handle the wrap-around bin [140, 180] -> [-180, -140]
        if phi_min == 140 and phi_max == 180:
            # This is bin 5: [140, 180]
            mask_phi = (phi_deg >= phi_min) & (phi_deg <= phi_max)
            bin_label = f"phi_{phi_min:+04d}_to_{phi_max:+04d}"
        elif phi_min == 180 and phi_max == -140:
            # This shouldn't happen with current edge definition
            continue
        elif phi_min > phi_max:
            # Wrap case: e.g., [140, -140] would need special handling
            # But with our edge definition, this is [-180, -140]
            mask_phi = (phi_deg >= phi_min) | (phi_deg <= phi_max)
            bin_label = f"phi_{phi_min:+04d}_to_{phi_max:+04d}"
        else:
            # Normal case
            mask_phi = (phi_deg >= phi_min) & (phi_deg < phi_max)
            bin_label = f"phi_{phi_min:+04d}_to_{phi_max:+04d}"

        # Combine masks
        mask_gse = mask_gse_base & mask_phi
        mask_pro = mask_pro_base & mask_phi
        mask_ret = mask_ret_base & mask_phi
        use = use_base & mask_phi

        # Count stars
        n_pro = np.sum(mask_pro & use)
        n_ret = np.sum(mask_ret & use)
        n_tot = np.sum(mask_gse & use)

        print(f"Bin {i+1}: phi in [{phi_min:+4d}, {phi_max:+4d}] deg")
        print(f"  GS/E stars: {n_tot:,} (pro: {n_pro:,}, ret: {n_ret:,})")

        if n_tot < 100:
            print(f"  Skipping: too few stars")
            continue

        n_stars_info = {'prograde': n_pro, 'retrograde': n_ret, 'combined': n_tot}
        outfile = f"plots/phi_bins/figure3_{bin_label}.png"

        make_chevron_plot(r, v_r, mask_pro, mask_ret, mask_gse, use,
                          phi_min, phi_max, n_stars_info, outfile)

    print("-" * 70)
    print("\nDone.")


if __name__ == '__main__':
    main()
