# JER revision package

This folder contains a conservative revision prepared from the supplied
manuscript and the values printed in its original Table 1.

## Before submission

1. Upload this package and the original spectroscopy exports to GitHub or Zenodo.
2. Replace `[DATA_REPOSITORY_URL_TO_BE_ADDED_BEFORE_SUBMISSION]` in the manuscript
   and title page with the public URL.
3. Verify `raw_data/manuscript_measurements.csv` against the instrument exports.
4. Confirm that `figure/Fig1.png` was produced without generative AI.
   Elsevier does not permit generative-AI artwork in submitted manuscripts.
5. Confirm author affiliations, funding text, and the corresponding-author address.

## Rebuild derived files

```powershell
& "C:\Users\kkaan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" scripts\rebuild_analysis.py
```

## Important metric distinction

- Paired recovery measures retention during solidification.
- KINS-reference bias measures absolute agreement with certified values.
- The revised manuscript does not conflate these metrics.
