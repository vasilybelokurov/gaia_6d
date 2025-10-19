# Gaia DR3 Full 6D Stellar Sample

Complete 6D phase-space catalog of Milky Way stars from Gaia DR3, combining astrometry, radial velocities, and photo-geometric distances.

**Author:** Vasily Belokurov
**Date:** 2025-10-19

---

## Overview

This project builds the largest full 6D stellar sample from Gaia DR3 combining:
- **Positions & proper motions:** Gaia DR3 astrometry
- **Radial velocities:** Gaia RVS (Radial Velocity Spectrometer)
- **Distances:** Bailer-Jones et al. (2021) photo-geometric distances

The resulting catalogs enable precision Galactic dynamics studies including phase-space analysis of stellar streams, halo substructure, and merger debris (e.g., Gaia Sausage/Enceladus).

---

## Project Structure

```
gaia_6d/
├── README.md                          # This file
├── CLAUDE.md                          # Project working instructions
│
├── build_gaia_6d_sample.py           # Query Gaia DR3 + BJ distances via WSDB
├── transform_to_galactocentric.py    # Convert to Galactocentric coordinates
├── sanity_test_vphi_vr.py            # V_phi vs V_R diagnostic plot
│
├── figure3_phase_space_folds.ipynb   # Reproduce Energy Wrinkles Fig. 3
│
├── papers/                            # Reference papers
│   └── Energy_wrinkles_and_phase_space_folds_Resubmit/
│       └── gse_gdr3_rvs.tex          # Belokurov+ (2022) paper
│
└── plots/                             # Output figures
    ├── vphi_vr_density.png           # Sanity test: disk rotation
    └── figure3_phase_space_folds.png # GS/E phase-space chevrons
```

---

## Data Products

All catalogs saved to `~/data/gaia/`:

### 1. `gaia_dr3_6d_full.fits` (5.9 GB)
**33,812,183 stars** with Gaia RVS radial velocities

**Columns (36):**
- Identifier: `source_id`
- Astrometry: `ra`, `dec`, `parallax`, `pmra`, `pmdec` + errors, correlations
- Photometry: `phot_g_mean_mag`, `phot_bp_mean_mag`, `phot_rp_mean_mag`, `bp_rp`
- Radial velocity: `radial_velocity`, `radial_velocity_error`, `rv_nb_transits`
- Distances (Bailer-Jones): `r_med_photogeo`, `r_lo_photogeo`, `r_hi_photogeo` (**parsecs**)
- Quality: `ruwe`, `astrometric_excess_noise`, `phot_bp_rp_excess_factor`, etc.
- Extinction: `ebv`

**Selection:**
- `radial_velocity IS NOT NULL`
- `radial_velocity_error IS NOT NULL`

### 2. `gaia_dr3_6d_galactocentric.fits` (4.5 GB)
**33,581,727 stars** (99.3%) with full 6D kinematics

**Columns (22):**
- Identifier: `source_id`
- Observables: `ra`, `dec`, `parallax`, `pmra`, `pmdec`, `distance` (pc), `rv` + errors
- Galactocentric Cartesian: `X`, `Y`, `Z` (kpc), `VX`, `VY`, `VZ` (km/s)
- Quality: `ruwe`, `ebv`

**Galactocentric frame:**
- Astropy default (R☉ = 8.122 kpc, Z☉ = 20.8 pc, v☉ from Schönrich+ 2010)

**Build time:** ~30 seconds (transformation @ 10.6M stars/s)

---

## Scripts

### 1. `build_gaia_6d_sample.py`

Query Gaia DR3 + Bailer-Jones distances via WSDB using `sqlutilpy`.

**Output:** `gaia_dr3_6d_full.fits`

**Usage:**
```bash
source ~/Work/venvs/.venv/bin/activate
python build_gaia_6d_sample.py
```

**Query time:** ~3.5 hours for 33.8M stars

**Key features:**
- Joins `gaia_dr3.gaia_source` with `gaia_edr3_aux.distances_bj`
- Saves full astrometry + photometry + quality indicators
- No quality cuts applied (for maximum flexibility)

---

### 2. `transform_to_galactocentric.py`

Transform Gaia observables → Galactocentric Cartesian coordinates.

**Input:** `gaia_dr3_6d_full.fits`
**Output:** `gaia_dr3_6d_galactocentric.fits`

**Usage:**
```bash
python transform_to_galactocentric.py
```

**Runtime:** ~30 seconds

**Features:**
- Filters stars with complete 6D data
- Uses `astropy.coordinates` for transformation
- Saves reduced catalog (observables + Galactocentric coords + quality)

---

### 3. `sanity_test_vphi_vr.py`

Diagnostic plot: V_φ vs V_R density for nearby stars.

**Quality cuts:**
- `distance < 10 kpc`
- `ruwe < 1.4`
- `rv_error < 20 km/s`

**Output:** `plots/vphi_vr_density.png`

**Expected results:**
- Median V_φ ≈ 219 km/s (disk rotation)
- Median V_R ≈ 0 km/s
- Dispersions ~50 km/s

**Usage:**
```bash
python sanity_test_vphi_vr.py
```

---

### 4. `figure3_phase_space_folds.ipynb`

Jupyter notebook reproducing Figure 3 from Belokurov et al. (2022).

**Shows:** GS/E phase-space chevrons in (v_r, r) space

**3×2 panel layout:**
- **Columns:** L_z > 0 (prograde) | L_z < 0 (retrograde) | Combined
- **Rows:** Column-normalized density | Background-subtracted

**Sample selection:**
- Heliocentric distance < 15 kpc
- |L_z| < 700 kpc km/s (GS/E angular momentum range)
- Removes globular cluster members (Vasiliev+ 2021, prob > 0)

**Output:** `plots/figure3_phase_space_folds.png`

**Sample size:** ~2.9M stars in GS/E range

**Usage:**
Run cells sequentially in Jupyter notebook.

---

## Key Results

### Full 6D Sample
- **33.6M stars** with complete phase-space information
- **99.3%** have Bailer-Jones distances
- **Median RUWE:** 1.02 (excellent astrometry)
- **Median RV precision:** 3.27 km/s
- **Median parallax S/N:** 32

### GS/E Phase-Space Folds (Figure 3)
- **2.9M stars** in GS/E angular momentum range
- **5 chevron patterns** detected in (v_r, r) space
- Prograde/retrograde asymmetry confirms dynamical complexity
- Background subtraction reveals incomplete phase-mixing

---

## Dependencies

**Python packages:**
- `numpy`
- `astropy`
- `matplotlib`
- `scipy`
- `sqlutilpy` (for WSDB queries)

**Data requirements:**
- Access to WSDB (Institute of Astronomy, Cambridge)
- Vasiliev+2021 GC membership catalog: `~/data/catalogues/gc_members_gaia_vasiliev.fits`

**Environment:**
```bash
source ~/Work/venvs/.venv/bin/activate
```

---

## Important Notes

### Distance Units
⚠️ **Bailer-Jones distances are in PARSECS, not kiloparsecs!**

The `gaia_edr3_aux.distances_bj` table columns:
- `r_med_photogeo`: median distance [pc]
- `r_lo_photogeo`, `r_hi_photogeo`: 16th/84th percentiles [pc]

Convert to kpc: `distance_kpc = r_med_photogeo / 1000.0`

### Galactocentric Frame

Default astropy frame used (not Drimmel+2022 from paper):
- R☉ = 8.122 kpc
- Z☉ = 20.8 pc
- v☉ from Schönrich+2010

**To match paper exactly**, modify `transform_to_galactocentric.py` to use:
```python
gc_frame = coord.Galactocentric(
    galcen_distance=8.0*u.kpc,
    z_sun=0.0*u.kpc,
    galcen_v_sun=[-9.3, 251.5, 8.59]*u.km/u.s
)
```

---

## References

**Papers:**
- Belokurov et al. (2022), MNRAS, "Energy wrinkles and phase-space folds of the last major merger"
  See `papers/Energy_wrinkles_and_phase_space_folds_Resubmit/gse_gdr3_rvs.tex`

**Data:**
- Gaia Collaboration (2022), "Gaia Data Release 3"
- Bailer-Jones et al. (2021), AJ, 161, 147 "Estimating Distances from Parallaxes. V."
- Vasiliev & Baumgardt (2021), MNRAS, 505, 5978 "Gaia EDR3 view on galactic globular clusters"

---

## Contact

For questions about this code:
- Vasily Belokurov (vasily@ast.cam.ac.uk)
- Institute of Astronomy, University of Cambridge

---

## License

Research code for academic use.
