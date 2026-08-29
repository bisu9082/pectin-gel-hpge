"""Repository-relative paths (no absolute paths, runs anywhere)."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT / 'spectra'
RAW = ROOT / 'raw_data'
FIG = ROOT / 'figures'
RAW.mkdir(exist_ok=True)
