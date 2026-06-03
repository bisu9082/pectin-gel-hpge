# pectin-gel-hpge

**Data and analysis code for:**

> Lee, Y., Kim, H., Shin, M., Choi, G. J., Yoo, J., & Kang, K.
> *Low-Methoxyl Pectin–Ca²⁺ Gelation as a Sample Solidification Protocol
> for Marine Radioactivity Monitoring by HPGe Gamma Spectrometry.*
> Journal of Environmental Radioactivity (submitted 2026).

Corresponding author: Ku Kang — bisu9082@gmail.com  
CBRN Defense Research Institute, Seoul 06796, Republic of Korea

---

## Repository contents

```
pectin-gel-hpge/
├── README.md
├── raw_data/
│   ├── manuscript_measurements.csv   ← activity concentrations (Table 1)
│   ├── processed_metrics.csv         ← computed R, Δ, B per sample
│   ├── recovery_summary.csv          ← per-nuclide summary statistics
│   └── co60_sensitivity_analysis.csv ← Co-60 power analysis output
├── scripts/
│   ├── rebuild_analysis.py           ← reproduces processed_metrics from raw_data
│   └── validate_submission.py        ← checks citation and figure completeness
└── validation/
    ├── calculation_audit.md          ← metric definitions and bias check
    ├── reference_audit.md            ← DOI resolver results and key corrections
    └── latex_build_status.md         ← static LaTeX validation log
```

---

## Data description

### `raw_data/manuscript_measurements.csv`

Activity concentrations (Bq kg⁻¹) for all radionuclide–sample combinations,
transcribed from the laboratory summary sheet used to prepare Table 1.

| Column | Description |
|---|---|
| `nuclide` | Radionuclide (137Cs, 134Cs, 241Am, 60Co) |
| `sample` | Sample identifier (S1–S5) |
| `energy_keV` | Principal gamma-ray energy |
| `KINS_certified_Bq_kg` | KINS proficiency-test certified reference value |
| `before_solidification_Bq_kg` | HPGe measurement before LM-pectin gelation |
| `after_solidification_Bq_kg` | HPGe measurement after LM-pectin gelation |
| `source` | Provenance note |

> **Important:** These values were transcribed from the manuscript. Before
> operational use, verify each value against the original spectroscopy export
> files from the Genie-2000 software. The spectroscopy exports should be
> deposited here alongside this table.

### `raw_data/processed_metrics.csv`

Derived metrics for each measured sample:
- `recovery_pct` = A_after / A_before × 100 (paired recovery R)
- `paired_change_pct` = R − 100 (Δ)
- `reference_bias_pct` = (A_after − A_KINS) / A_KINS × 100 (B)

### `raw_data/recovery_summary.csv`

Per-nuclide summary: mean, SD, min, max recovery; fraction within 80–120%.

### `raw_data/co60_sensitivity_analysis.csv`

Monte Carlo sensitivity analysis for the ⁶⁰Co n = 3 dataset
(500 000 simulations, seed 20260531). Reports 95% CI and minimum detectable
shift at 80% power.

---

## Reproducing the analysis

```bash
pip install pandas numpy scipy
python scripts/rebuild_analysis.py
```

Output files are written to `raw_data/`. Requires Python ≥ 3.9.

---

## Key metric distinction

This study reports **two distinct metrics** that must not be conflated:

| Metric | Formula | Interpretation |
|---|---|---|
| Paired recovery R | A_after / A_before × 100 | Effect of solidification step |
| KINS-reference bias B | (A_after − A_KINS) / A_KINS × 100 | Absolute calibration accuracy |

R evaluates whether the solidification protocol preserves activity.
B reflects all differences between the measurement system and KINS
certified values, including pre-existing calibration offsets.
Four B values exceed ±20% (see `validation/calculation_audit.md`);
this is discussed as a limitation in the manuscript.

---

## Spectroscopy exports (to be deposited)

The original Genie-2000 spectroscopy export files (Canberra format)
for all before- and after-solidification measurements should be
deposited in `raw_data/spectroscopy_exports/` before finalising the
repository for peer-review disclosure.

---

## License

Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)  
Code: [MIT License](https://opensource.org/licenses/MIT)

---

## Citation

If you use these data or scripts, please cite the manuscript above.
BibTeX entry will be added upon journal assignment of a DOI.
