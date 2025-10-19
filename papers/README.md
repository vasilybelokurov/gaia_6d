# Reference Papers

This directory contains LaTeX source and figures for reference papers related to this project.

---

## Energy_wrinkles_and_phase_space_folds_Resubmit/

**Paper:** Belokurov et al. (2022), MNRAS, "Energy wrinkles and phase-space folds of the last major merger"

**Contents:**
- `gse_gdr3_rvs.tex` - Main LaTeX source
- `references.bib` - Bibliography
- `mnras.cls`, `mnras.bst` - MNRAS journal style files
- `img/` - All paper figures (PDFs)
  - Figure 1 (`elz_reveal.pdf`) - Energy-angular momentum distribution
  - Figure 2 (`phase_space_fig1.pdf`, `phase_space_fig2.pdf`, `phase_space_fig3.pdf`) - Phase-space views
  - Figure 3 (**phase_space_fig2.pdf**) - **GS/E phase-space folds** (reproduced in this project)
  - Figure 4 (`maxvr_r.pdf`) - Phase-space structure sketch
  - Figure 5 (`phase_space_asym.pdf`) - Asymmetry in phase-space folds
  - Figure 6 (`phase_space_feh.pdf`) - Metallicity dependence
  - Auriga simulation figures (`Auriga_5_*.pdf`)
  - Merger simulation figures (`mergersim_*.pdf`)

**Key result reproduced:** Figure 3 (phase_space_fig2.pdf) showing chevron patterns in GS/E debris

**ADS Link:** https://ui.adsabs.harvard.edu/abs/2022MNRAS.514..689B

---

## Usage

To compile the LaTeX paper:
```bash
cd Energy_wrinkles_and_phase_space_folds_Resubmit/
pdflatex gse_gdr3_rvs.tex
bibtex gse_gdr3_rvs
pdflatex gse_gdr3_rvs.tex
pdflatex gse_gdr3_rvs.tex
```

This project reproduces **Figure 3** from this paper using the Jupyter notebook `figure3_phase_space_folds.ipynb`.
