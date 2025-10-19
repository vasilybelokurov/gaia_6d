#!/usr/bin/env python3
"""
Build full 6D stellar sample from Gaia DR3 only:
- Positions and proper motions from Gaia DR3 astrometry
- Radial velocities from Gaia RVS
- Distances from Bailer-Jones (photo-geometric)

Author: VB
Date: 2025-10-18
"""

import numpy as np
import sqlutilpy as sqlutil
from astropy.table import Table
import os


def build_gaia_6d_query(test_mode=False):
    """
    Build SQL query for full 6D Gaia sample.

    Parameters
    ----------
    test_mode : bool
        If True, limit to 1 degree cone for testing

    Returns
    -------
    query : str
        SQL query string
    """

    query = """
    SELECT
        -- Source identifier
        g.source_id,

        -- Position (ICRS, epoch 2016.0)
        g.ra,
        g.dec,

        -- Proper motions (mas/yr)
        g.pmra,
        g.pmdec,
        g.pmra_error,
        g.pmdec_error,
        g.pmra_pmdec_corr,

        -- Parallax (mas)
        g.parallax,
        g.parallax_error,
        g.parallax_over_error,

        -- Radial velocity (km/s)
        g.radial_velocity,
        g.radial_velocity_error,
        g.rv_nb_transits,
        g.rv_expected_sig_to_noise,

        -- Photometry
        g.phot_g_mean_mag,
        g.phot_bp_mean_mag,
        g.phot_rp_mean_mag,
        g.bp_rp,
        g.phot_g_mean_flux_over_error,
        g.phot_bp_mean_flux_over_error,
        g.phot_rp_mean_flux_over_error,
        g.phot_bp_rp_excess_factor,

        -- Astrometric quality indicators
        g.ruwe,
        g.astrometric_excess_noise,
        g.astrometric_excess_noise_sig,
        g.astrometric_chi2_al,
        g.astrometric_n_good_obs_al,
        g.astrometric_gof_al,
        g.astrometric_matched_transits,
        g.visibility_periods_used,

        -- Extinction
        g.ebv,

        -- Bailer-Jones distances (pc)
        bj.r_med_photogeo,
        bj.r_lo_photogeo,
        bj.r_hi_photogeo,
        bj.flag AS distance_flag

    FROM gaia_dr3.gaia_source AS g
    LEFT JOIN gaia_edr3_aux.distances_bj AS bj ON g.source_id = bj.source_id

    WHERE
        g.radial_velocity IS NOT NULL
        AND g.radial_velocity_error IS NOT NULL
    """

    if test_mode:
        # Test on 1 degree cone around arbitrary point (l=180, b=45)
        query += """
        AND q3c_radial_query(g.ra, g.dec, 180.0, 45.0, 1.0)
        """

    return query


def execute_query(query, test_mode=False):
    """
    Execute query via sqlutilpy and return astropy Table.

    Parameters
    ----------
    query : str
        SQL query
    test_mode : bool
        If True, print query info

    Returns
    -------
    table : astropy.table.Table
        Query results
    """
    if test_mode:
        print("Executing test query (1 degree cone)...")
        print(f"Query length: {len(query)} characters\n")

    data = sqlutil.get(query, asDict=True)

    n_stars = len(data['source_id']) if 'source_id' in data else 0
    print(f"Retrieved {n_stars:,} stars")

    return Table(data)


def save_catalog(table, filename, overwrite=True):
    """
    Save catalog to FITS file.

    Parameters
    ----------
    table : astropy.table.Table
        Data to save
    filename : str
        Output filename
    overwrite : bool
        Overwrite existing file
    """
    os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
    table.write(filename, overwrite=overwrite)
    print(f"Saved to: {filename}")


def print_sample_stats(table):
    """Print diagnostic statistics from the sample."""
    print("\n=== Sample Statistics ===")
    print(f"Total stars: {len(table):,}")

    # RV coverage
    has_rv = np.isfinite(table['radial_velocity'])
    print(f"Stars with RV: {np.sum(has_rv):,} ({100*np.sum(has_rv)/len(table):.1f}%)")

    # Distance coverage
    has_dist = np.isfinite(table['r_med_photogeo'])
    print(f"Stars with BJ distance: {np.sum(has_dist):,} ({100*np.sum(has_dist)/len(table):.1f}%)")

    # Full 6D
    has_6d = has_rv & has_dist & np.isfinite(table['pmra']) & np.isfinite(table['parallax'])
    print(f"Stars with full 6D: {np.sum(has_6d):,} ({100*np.sum(has_6d)/len(table):.1f}%)")

    # Quality indicators
    if np.sum(has_6d) > 0:
        print(f"\nQuality stats (6D sample):")
        print(f"  RUWE median: {np.median(table['ruwe'][has_6d]):.2f}")
        print(f"  RV error median: {np.median(table['radial_velocity_error'][has_6d]):.2f} km/s")
        print(f"  Parallax S/N median: {np.median(table['parallax_over_error'][has_6d]):.1f}")
        print(f"  Distance range: {np.min(table['r_med_photogeo'][has_6d]):.1f} - {np.max(table['r_med_photogeo'][has_6d]):.1f} pc")


def main():
    """Main execution."""

    # FULL SKY MODE
    print("=" * 60)
    print("Building Gaia 6D sample (FULL SKY)")
    print("=" * 60)

    query = build_gaia_6d_query(test_mode=False)
    table = execute_query(query, test_mode=False)

    # Save catalog
    outfile = os.path.expanduser('~/data/gaia/gaia_dr3_6d_full.fits')
    save_catalog(table, outfile)

    # Print diagnostics
    print_sample_stats(table)

    print("\n" + "=" * 60)
    print("Full-sky catalog complete")
    print("=" * 60)


if __name__ == '__main__':
    main()
