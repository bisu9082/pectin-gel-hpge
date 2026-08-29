# pectin-gel-hpge — data and code package

Data and analysis code for: *Low-Methoxyl Pectin–Ca2+ Gelation as a Sample
Solidification Protocol for Marine Radioactivity Monitoring by HPGe Gamma
Spectrometry* (Journal of Environmental Radioactivity, JENVRAD-D-26-00429).

## Contents
- `spectra/before`, `spectra/after` — Genie-2000 CNF spectra of the original
  five paired preparations (S1–S5), acquired March 2026, 30,000 s live time.
- `spectra/replicates_before`, `spectra/replicates_after` — replicate paired
  preparations S6–S8 (2022-round KINS reference material), March 2026
  (S6–S8; S8 before-measurement 80,000 s).
- `spectra/efficiency_calibration` — certified aqueous multi-nuclide standard
  (Eckert & Ziegler via KINS, ref. 2025-04-01) measured in the same container
  at fill heights 2/5/10/20/30/40 mm (PM1EFF255D02–D40), May 2025.
- `calibration` — energy-calibration (CAL) and certificate (CTF) files.
- `raw_data` — `dataset_retention.csv` (definitive 26-pair dataset:
  activities ± u, MDA, masses, paired retention ± u), `peaks_all.csv`
  (independent net-peak areas, 6 lines × 26 spectra), `recovery_lines.csv`
  (line-by-line retentions), `efficiency_measured.json` (measured
  efficiencies per fill height).
- `scripts` — `cnf_reader.py` (CNF/CAM parser), `spectrum_analysis.py`
  (net-peak-area, Currie levels, line-by-line retention),
  `rebuild_stats.py` (dataset statistics), `make_figures_v2.py` (Figs. 4–5).

## Reproducing the manuscript values
1. `python scripts/spectrum_analysis.py` — independent peak analysis of all
   26 spectra (Table S2, MDA cross-checks).
2. `python scripts/rebuild_stats.py` — Table 1/2 statistics, k=2 consistency,
   activity-level comparison.
3. `python scripts/make_figures_v2.py` — Figures 4 and 5.

Python ≥3.10 with numpy, scipy, matplotlib, openpyxl.
