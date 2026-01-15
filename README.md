# Gaia DR3 Full 6D Stellar Sample

Complete 6D phase-space catalog of Milky Way stars from Gaia DR3, combining astrometry, radial velocities, and photo-geometric distances.

**Authors:** Vasily Belokurov, Alice Archer, Adam Dillamore, Hanyuan Zhang
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
├── figure3_phi_bins.py               # Phase-space folds by azimuthal angle
├── plot_phi_histogram.py             # Azimuthal angle distribution
│
├── papers/                            # Reference papers
│   └── Energy_wrinkles_and_phase_space_folds_Resubmit/
│       └── gse_gdr3_rvs.tex          # Belokurov+ (2022) paper
│
└── plots/                             # Output figures
    ├── vphi_vr_density.png           # Sanity test: disk rotation
    ├── figure3_phase_space_folds.png # GS/E phase-space chevrons
    ├── phi_histogram.png             # Azimuthal angle distribution
    └── phi_bins/                      # Chevrons by phi bin (40 deg)
```

---

## Data Products

### Download Links

The Galactocentric catalog is available for download:

**https://people.ast.cam.ac.uk/~vasily/data/gaia/dr3/**

- `gaia_dr3_6d_galactocentric.fits` (4.5 GB) - 33.6M stars with full 6D + Galactocentric coords

Also available at the same location:
- `gc_members_gaia_vasiliev.fits` (212 MB) - Vasiliev+2021 GC membership catalog
- `RRL_GDR3.fits` (442 MB) - Gaia DR3 RR Lyrae

Place downloaded files in `~/data/gaia/` to use with the scripts.

**Note:** `gaia_dr3_6d_full.fits` (5.9 GB) is not publicly hosted due to size. Build it locally using `build_gaia_6d_sample.py` if needed.

---

### Catalog Descriptions

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

**Columns (26):**
- Identifier: `source_id`
- Observables: `ra`, `dec`, `parallax`, `pmra`, `pmdec`, `distance` (pc), `rv` + errors
- Galactocentric Cartesian: `X`, `Y`, `Z` (kpc), `VX`, `VY`, `VZ` (km/s)
- Galactocentric Cylindrical: `R_cyl` (kpc), `phi` (rad), `V_R`, `V_phi` (km/s)
- Quality: `ruwe`, `ebv`

**Cylindrical coordinate convention:**
- `phi = 0` toward the Sun, increasing in prograde direction
- `V_phi > 0` for prograde rotation (median ≈ 220 km/s for disk stars)

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
python build_gaia_6d_sample.py
```

**Query time:** ~3.5 hours for 33.8M stars

**Key features:**
- Joins `gaia_dr3.gaia_source` with `gaia_edr3_aux.distances_bj`
- Saves full astrometry + photometry + quality indicators
- No quality cuts applied (for maximum flexibility)

---

### 2. `transform_to_galactocentric.py`

Transform Gaia observables → Galactocentric coordinates (Cartesian + cylindrical).

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
- Computes Cartesian (X, Y, Z, VX, VY, VZ) and cylindrical (R_cyl, phi, V_R, V_phi)
- Cylindrical convention: phi=0 toward Sun, V_phi>0 prograde

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

### 5. `plot_phi_histogram.py`

Diagnostic plot of Galactocentric azimuthal angle distribution.

**Output:** `plots/phi_histogram.png`

**Usage:**
```bash
python plot_phi_histogram.py
```

Shows strong concentration at phi ≈ 0 (toward Sun) due to RVS magnitude limit.

---

### 6. `figure3_phi_bins.py`

Phase-space folds split by azimuthal angle (40° bins starting at φ = −20°).

**Output:** `plots/phi_bins/figure3_phi_*.png`

**Usage:**
```bash
python figure3_phi_bins.py
```

**Results:** Only 3 of 9 bins populated due to RVS magnitude-limit geometry:
- φ ∈ [−20°, +20°]: 236k stars (strongest signal)
- φ ∈ [+20°, +60°]: 24k stars
- φ ∈ [−60°, −20°]: 31k stars

Stars at |φ| > 60° are too distant for RVS (median d > 8 kpc).

---

## Key Results

### Full 6D Sample
- **33.6M stars** with complete phase-space information
- **99.3%** have Bailer-Jones distances
- **Median RUWE:** 1.02 (excellent astrometry)
- **Median RV precision:** 3.27 km/s
- **Median parallax S/N:** 32

### GS/E Phase-Space Folds (Figure 3)
- **2.95M stars** in GS/E angular momentum range (|L_z| < 700 kpc km/s)
  - Prograde (L_z > 0): 428,725 stars
  - Retrograde (L_z < 0): 2,525,256 stars
- **292,504 stars** after quality cuts (9.9% retention)
  - Cuts: d_helio < 15 kpc, distance S/N > 3, parallax S/N > 10, E(B-V) < 100
  - Prograde: 59,195 stars
  - Retrograde: 233,309 stars
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

## Authors

- **Vasily Belokurov** (vasily@ast.cam.ac.uk) – Institute of Astronomy, University of Cambridge
- **Alice Archer** (aeaa2@cam.ac.uk) – Institute of Astronomy, University of Cambridge
- **Adam Dillamore** (a.dillamore@ucl.ac.uk) – University College London
- **Hanyuan Zhang** (hz420@cam.ac.uk) – Institute of Astronomy, University of Cambridge

---

## Contact

For questions about this code, contact Vasily Belokurov (vasily@ast.cam.ac.uk)

---

## License

Research code for academic use.
