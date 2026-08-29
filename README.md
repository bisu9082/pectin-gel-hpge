# pectin-gel-hpge — data and code

Measurement data and analysis code supporting:

> **Low-Methoxyl Pectin–Ca²⁺ Gelation as a Sample Solidification Protocol for
> Marine Radioactivity Monitoring by HPGe Gamma Spectrometry**
> *Journal of Environmental Radioactivity* (manuscript JENVRAD-D-26-00429)

Every numerical result in the manuscript can be regenerated from the files in
this repository with the scripts in `scripts/`.

## Contents

| Path | Description |
|---|---|
| `spectra/before/` | Genie-2000 spectra before solidification, `S1`–`S8` (8 files) |
| `spectra/after/` | Genie-2000 spectra after solidification, `S1`–`S8` (8 files) |
| `spectra/efficiency_calibration/` | Certified aqueous multi-nuclide standard measured in the same container at fill heights 2, 5, 10, 20, 30 and 40 mm (`PM1EFF255D02`–`D40`) |
| `calibration/` | Energy-calibration file (`PM1STD25DH02.CAL`) and source certificate (`KINS2025_DH_2mm.CTF`) |
| `raw_data/dataset_retention.csv` | Definitive 26-pair dataset: activity concentrations ± *u* (*k* = 1), MDA, sample masses, paired retention ± *u* |
| `raw_data/peaks_all.csv` | Independent net-peak areas, Currie limits and count rates (6 gamma lines × 16 sample spectra) |
| `raw_data/recovery_lines.csv` | Line-by-line paired retention (primary and secondary gamma lines) |
| `raw_data/efficiency_measured.json` | Measured full-energy-peak efficiency at each fill height |
| `figures/` | Figures 4 and 5 of the manuscript and Figure S1 of the Supplementary Information |
| `scripts/` | Analysis code (see below) |

Samples `S1`–`S5` are the original preparations; `S6`–`S8` are the replicate
preparations of the 2022 KINS proficiency-test round. All measurements used
30,000 s live time except the `S8` before-solidification spectrum (80,000 s).

## Reproducing the manuscript

```bash
pip install -r requirements.txt
cd scripts
python spectrum_analysis.py   # peak analysis of all 26 pairs -> Table S2, MDA cross-check
python rebuild_stats.py       # Tables 1-3 and the statistics of Sections 3.3-3.5
python efficiency_curve.py    # measured efficiency calibration -> Figure S1
python make_figures.py        # Figures 4 and 5
```

`spectrum_analysis.py` re-integrates every spectrum independently of the
Genie-2000 results, so it serves as a check on the tabulated values rather
than a restatement of them. `make_figures.py` regenerates Figures 4 and 5
byte-for-byte as submitted.

Requires Python ≥ 3.10.

## Scripts

| File | Purpose |
|---|---|
| `_paths.py` | Repository-relative paths (no absolute paths; runs from any location) |
| `cnf_reader.py` | Reader for Canberra CNF/CAM spectrum files |
| `spectrum_analysis.py` | Net peak areas, counting uncertainties, Currie critical level and detection limit, line-by-line retention |
| `rebuild_stats.py` | Paired-retention dataset, per-nuclide statistics, *k* = 2 consistency, activity-level comparison, ⁶⁰Co sensitivity analysis |
| `efficiency_curve.py` | Empirical efficiency calibration from the certified standard |
| `make_figures.py` | Manuscript Figures 4 and 5 |

## License

Data: CC BY 4.0. Code: MIT. See `LICENSE`.
