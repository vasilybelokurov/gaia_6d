#!/usr/bin/env python3
"""
Transform Gaia 6D sample to Galactocentric coordinates.

Input: gaia_dr3_6d_full.fits (33M stars with Gaia observables)
Output: gaia_dr3_6d_galactocentric.fits (reduced catalog with Galactocentric coords)

Author: VB
Date: 2025-10-18
"""

import numpy as np
import time
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
import os


def transform_to_galactocentric(input_file, output_file):
    """
    Transform Gaia observables to Galactocentric coordinates.

    Parameters
    ----------
    input_file : str
        Input FITS catalog with Gaia observables
    output_file : str
        Output FITS catalog with Galactocentric coordinates
    """

    print("=" * 70)
    print("GAIA 6D → GALACTOCENTRIC TRANSFORMATION")
    print("=" * 70)

    # ========================================================================
    # Load catalog
    # ========================================================================

    print(f"\nLoading catalog: {input_file}")
    t0_load = time.time()

    cat = Table.read(input_file)
    n_total = len(cat)

    t_load = time.time() - t0_load
    print(f"  Loaded {n_total:,} stars in {t_load:.1f}s")

    # ========================================================================
    # Filter stars with full 6D data
    # ========================================================================

    print("\nFiltering stars with full 6D data...")

    has_pm = np.isfinite(cat['pmra']) & np.isfinite(cat['pmdec'])
    has_rv = np.isfinite(cat['radial_velocity'])
    has_dist = np.isfinite(cat['r_med_photogeo'])

    mask_6d = has_pm & has_rv & has_dist
    n_6d = np.sum(mask_6d)

    print(f"  Stars with proper motions: {np.sum(has_pm):,}")
    print(f"  Stars with radial velocity: {np.sum(has_rv):,}")
    print(f"  Stars with BJ distance: {np.sum(has_dist):,}")
    print(f"  Stars with full 6D: {n_6d:,} ({100*n_6d/n_total:.1f}%)")

    cat_6d = cat[mask_6d]

    # ========================================================================
    # Transform to Galactocentric
    # ========================================================================

    print(f"\nTransforming {n_6d:,} stars to Galactocentric coordinates...")
    print("  Using astropy default Galactocentric frame")

    t0_transform = time.time()

    # Create SkyCoord object (BJ distances are in parsecs!)
    c_icrs = coord.SkyCoord(
        ra=cat_6d['ra'] * u.deg,
        dec=cat_6d['dec'] * u.deg,
        distance=cat_6d['r_med_photogeo'] * u.pc,
        pm_ra_cosdec=cat_6d['pmra'] * u.mas/u.yr,
        pm_dec=cat_6d['pmdec'] * u.mas/u.yr,
        radial_velocity=cat_6d['radial_velocity'] * u.km/u.s,
        frame='icrs'
    )

    # Transform to Galactocentric (using astropy defaults)
    c_gal = c_icrs.transform_to(coord.Galactocentric())

    # Extract Cartesian coordinates and velocities
    X = c_gal.x.to_value(u.kpc)
    Y = c_gal.y.to_value(u.kpc)
    Z = c_gal.z.to_value(u.kpc)
    VX = c_gal.v_x.to_value(u.km/u.s)
    VY = c_gal.v_y.to_value(u.km/u.s)
    VZ = c_gal.v_z.to_value(u.km/u.s)

    t_transform = time.time() - t0_transform
    print(f"  Transformation complete in {t_transform:.1f}s ({n_6d/t_transform:.0f} stars/s)")

    # ========================================================================
    # Compute distance error (in pc)
    # ========================================================================

    print("\nComputing distance errors...")
    distance_error = 0.5 * (cat_6d['r_hi_photogeo'] - cat_6d['r_lo_photogeo'])
    distance_pc = cat_6d['r_med_photogeo']

    # ========================================================================
    # Build output catalog
    # ========================================================================

    print("\nBuilding output catalog...")

    out_cat = Table()

    # Identifier
    out_cat['source_id'] = cat_6d['source_id']

    # Observables
    out_cat['ra'] = cat_6d['ra']
    out_cat['dec'] = cat_6d['dec']
    out_cat['parallax'] = cat_6d['parallax']
    out_cat['parallax_error'] = cat_6d['parallax_error']
    out_cat['pmra'] = cat_6d['pmra']
    out_cat['pmdec'] = cat_6d['pmdec']
    out_cat['pmra_error'] = cat_6d['pmra_error']
    out_cat['pmdec_error'] = cat_6d['pmdec_error']
    out_cat['pmra_pmdec_corr'] = cat_6d['pmra_pmdec_corr']
    out_cat['distance'] = distance_pc
    out_cat['distance_error'] = distance_error
    out_cat['rv'] = cat_6d['radial_velocity']
    out_cat['rv_error'] = cat_6d['radial_velocity_error']

    # Galactocentric coordinates
    out_cat['X'] = X
    out_cat['Y'] = Y
    out_cat['Z'] = Z
    out_cat['VX'] = VX
    out_cat['VY'] = VY
    out_cat['VZ'] = VZ

    # Quality flags
    out_cat['ruwe'] = cat_6d['ruwe']
    out_cat['ebv'] = cat_6d['ebv']

    # Add units and descriptions
    out_cat['ra'].unit = u.deg
    out_cat['dec'].unit = u.deg
    out_cat['parallax'].unit = u.mas
    out_cat['parallax_error'].unit = u.mas
    out_cat['pmra'].unit = u.mas/u.yr
    out_cat['pmdec'].unit = u.mas/u.yr
    out_cat['pmra_error'].unit = u.mas/u.yr
    out_cat['pmdec_error'].unit = u.mas/u.yr
    out_cat['distance'].unit = u.pc
    out_cat['distance_error'].unit = u.pc
    out_cat['rv'].unit = u.km/u.s
    out_cat['rv_error'].unit = u.km/u.s
    out_cat['X'].unit = u.kpc
    out_cat['Y'].unit = u.kpc
    out_cat['Z'].unit = u.kpc
    out_cat['VX'].unit = u.km/u.s
    out_cat['VY'].unit = u.km/u.s
    out_cat['VZ'].unit = u.km/u.s

    print(f"  Output catalog: {len(out_cat.colnames)} columns, {len(out_cat):,} rows")

    # ========================================================================
    # Print statistics
    # ========================================================================

    print("\n=== Galactocentric Statistics ===")
    print(f"X range: {np.min(X):.2f} to {np.max(X):.2f} kpc")
    print(f"Y range: {np.min(Y):.2f} to {np.max(Y):.2f} kpc")
    print(f"Z range: {np.min(Z):.2f} to {np.max(Z):.2f} kpc")
    print(f"VX range: {np.min(VX):.1f} to {np.max(VX):.1f} km/s")
    print(f"VY range: {np.min(VY):.1f} to {np.max(VY):.1f} km/s")
    print(f"VZ range: {np.min(VZ):.1f} to {np.max(VZ):.1f} km/s")

    R_gal = np.sqrt(X**2 + Y**2)
    print(f"\nGalactocentric radius: {np.min(R_gal):.2f} to {np.max(R_gal):.2f} kpc")
    print(f"Height |Z|: {np.min(np.abs(Z)):.2f} to {np.max(np.abs(Z)):.2f} kpc")

    # ========================================================================
    # Save catalog
    # ========================================================================

    print(f"\nSaving catalog: {output_file}")
    t0_save = time.time()

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    out_cat.write(output_file, overwrite=True)

    t_save = time.time() - t0_save
    file_size_gb = os.path.getsize(output_file) / 1e9

    print(f"  Saved in {t_save:.1f}s ({file_size_gb:.2f} GB)")

    # ========================================================================
    # Summary
    # ========================================================================

    t_total = t_load + t_transform + t_save

    print("\n" + "=" * 70)
    print("TIMING SUMMARY")
    print("=" * 70)
    print(f"  Load catalog:     {t_load:8.1f}s")
    print(f"  Transform coords: {t_transform:8.1f}s ({n_6d/t_transform:,.0f} stars/s)")
    print(f"  Save catalog:     {t_save:8.1f}s")
    print(f"  TOTAL:            {t_total:8.1f}s ({t_total/60:.1f} min)")
    print("=" * 70)


def main():
    """Main execution."""

    input_file = os.path.expanduser('~/data/gaia/gaia_dr3_6d_full.fits')
    output_file = os.path.expanduser('~/data/gaia/gaia_dr3_6d_galactocentric.fits')

    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        return

    transform_to_galactocentric(input_file, output_file)

    print("\nDone.")


if __name__ == '__main__':
    main()
