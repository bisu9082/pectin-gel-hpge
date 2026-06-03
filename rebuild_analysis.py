from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "raw_data"
FIG_DIR = ROOT / "figure"
VAL_DIR = ROOT / "validation"

ROWS = [
    ("137Cs", "S1", 45.1, 44.0, 44.4),
    ("137Cs", "S2", 45.1, 47.0, 43.4),
    ("137Cs", "S3", 24.6, 30.8, 30.7),
    ("137Cs", "S4", 24.6, 28.3, 30.6),
    ("137Cs", "S5", 24.8, 30.3, 29.9),
    ("134Cs", "S1", 45.1, 44.0, 40.2),
    ("134Cs", "S2", 45.1, 47.0, 47.3),
    ("134Cs", "S3", 24.6, 30.8, 28.4),
    ("134Cs", "S4", 24.6, 28.3, 20.1),
    ("241Am", "S1", 45.1, 44.0, 40.6),
    ("241Am", "S2", 45.1, 47.0, 45.8),
    ("241Am", "S3", 24.6, 30.8, 26.5),
    ("241Am", "S4", 24.6, 28.3, 26.1),
    ("241Am", "S5", 24.8, 30.3, 30.5),
    ("60Co", "S1", 45.1, 44.0, 40.5),
    ("60Co", "S2", 45.1, 47.0, 44.6),
    ("60Co", "S5", 24.8, 30.3, 29.5),
]

ENERGY_KEV = {"137Cs": 661.7, "134Cs": 604.7, "241Am": 59.5, "60Co": 1173.2}
ORDER = ["137Cs", "134Cs", "241Am", "60Co"]
COLORS = {"137Cs": "#C94F4A", "134Cs": "#E8943A", "241Am": "#4AACB0", "60Co": "#5B8DB8"}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def mean_sd(values: list[float]) -> tuple[float, float]:
    return float(np.mean(values)), float(np.std(values, ddof=1))


def t_critical_975(df: int) -> float:
    values = {1: 12.7062047364, 2: 4.3026527299, 3: 3.1824463053, 4: 2.7764451052}
    return values[df]


def build_data() -> tuple[list[dict], list[dict]]:
    raw_rows = []
    processed_rows = []
    for nuclide, sample, certified, before, after in ROWS:
        recovery = after / before * 100
        paired_change = recovery - 100
        reference_bias = (after - certified) / certified * 100
        raw_rows.append(
            {
                "nuclide": nuclide,
                "sample": sample,
                "energy_keV": ENERGY_KEV[nuclide],
                "KINS_certified_Bq_kg": certified,
                "before_solidification_Bq_kg": before,
                "after_solidification_Bq_kg": after,
                "source": "Values transcribed from manuscript Table 1; verify against instrument export before submission",
            }
        )
        processed_rows.append(
            {
                "nuclide": nuclide,
                "sample": sample,
                "recovery_percent_after_over_before": round(recovery, 3),
                "paired_change_percent": round(paired_change, 3),
                "reference_bias_percent_after_vs_KINS": round(reference_bias, 3),
                "paired_recovery_within_80_120": "yes" if abs(paired_change) <= 20 else "no",
                "reference_bias_within_plusminus20": "yes" if abs(reference_bias) <= 20 else "no",
            }
        )
    return raw_rows, processed_rows


def build_summary(processed_rows: list[dict]) -> list[dict]:
    summary_rows = []
    for nuclide in ORDER:
        values = [
            float(row["recovery_percent_after_over_before"])
            for row in processed_rows
            if row["nuclide"] == nuclide
        ]
        avg, sd = mean_sd(values)
        summary_rows.append(
            {
                "nuclide": nuclide,
                "energy_keV": ENERGY_KEV[nuclide],
                "n": len(values),
                "mean_recovery_percent": round(avg, 3),
                "sample_sd_percent": round(sd, 3),
                "min_recovery_percent": round(min(values), 3),
                "max_recovery_percent": round(max(values), 3),
            }
        )
    return summary_rows


def build_co60_sensitivity(summary_rows: list[dict]) -> list[dict]:
    row = next(row for row in summary_rows if row["nuclide"] == "60Co")
    mean = float(row["mean_recovery_percent"])
    sd = float(row["sample_sd_percent"])
    n = int(row["n"])
    df = n - 1
    critical = t_critical_975(df)
    ci_half_width = critical * sd / math.sqrt(n)
    rng = np.random.default_rng(20260531)
    draws = rng.normal(0, sd, size=(500_000, n))

    def simulated_power(delta: float) -> float:
        samples = draws + delta
        test_statistics = samples.mean(axis=1) / (samples.std(axis=1, ddof=1) / math.sqrt(n))
        return float(np.mean(np.abs(test_statistics) > critical))

    low, high = 0.0, 100.0
    for _ in range(40):
        mid = (low + high) / 2
        if simulated_power(mid) < 0.80:
            low = mid
        else:
            high = mid
    detectable_shift = high
    return [
        {
            "nuclide": "60Co",
            "n": n,
            "mean_recovery_percent": round(mean, 3),
            "sample_sd_percent": round(sd, 3),
            "ci95_lower_percent": round(mean - ci_half_width, 3),
            "ci95_upper_percent": round(mean + ci_half_width, 3),
            "two_sided_alpha": 0.05,
            "target_power": 0.80,
            "minimum_detectable_absolute_shift_percentage_points": round(detectable_shift, 3),
            "method": "Monte Carlo one-sample two-sided t-test; 500000 draws; seed 20260531",
            "interpretation": "Descriptive sensitivity analysis based on observed SD; not a substitute for prospective replication",
        }
    ]


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
            "figure.dpi": 150,
        }
    )


def make_fig4(raw_rows: list[dict], summary_rows: list[dict]) -> None:
    style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax = axes[0]
    cs_rows = [row for row in raw_rows if row["nuclide"] == "137Cs"]
    x = np.arange(len(cs_rows))
    before = np.array([row["before_solidification_Bq_kg"] for row in cs_rows])
    after = np.array([row["after_solidification_Bq_kg"] for row in cs_rows])
    ax.scatter(x - 0.07, before, s=55, facecolors="white", edgecolors=COLORS["137Cs"], linewidth=1.6, label="Before solidification")
    ax.scatter(x + 0.07, after, s=55, color=COLORS["137Cs"], label="After solidification")
    for xpos, b, a in zip(x, before, after):
        ax.plot([xpos - 0.07, xpos + 0.07], [b, a], color="#888888", linewidth=0.9)
    ax.set_xticks(x, [row["sample"] for row in cs_rows])
    ax.set_ylabel(r"$^{137}$Cs activity concentration (Bq kg$^{-1}$)")
    ax.set_xlabel("Sample")
    ax.legend(frameon=False, loc="upper right")
    ax.text(-0.10, 1.04, "(a)", transform=ax.transAxes, fontsize=15, fontweight="bold")

    ax = axes[1]
    means = [row["mean_recovery_percent"] for row in summary_rows]
    sds = [row["sample_sd_percent"] for row in summary_rows]
    labels = [row["nuclide"] for row in summary_rows]
    colors = [COLORS[label] for label in labels]
    xx = np.arange(len(labels))
    ax.bar(xx, means, yerr=sds, color=colors, width=0.58, capsize=4)
    ax.axhline(100, color="#555555", linestyle="--", linewidth=1.1, label="100% recovery")
    ax.axhline(80, color="#888888", linestyle=":", linewidth=1.0)
    ax.axhline(120, color="#888888", linestyle=":", linewidth=1.0, label="Paired recovery range (80-120%)")
    ax.set_xticks(xx, labels)
    ax.set_ylim(60, 130)
    ax.set_ylabel("Paired recovery (%)")
    for xpos, row in zip(xx, summary_rows):
        ax.text(xpos, row["mean_recovery_percent"] + row["sample_sd_percent"] + 2, f'{row["mean_recovery_percent"]:.1f}%\n(n={row["n"]})', ha="center", fontsize=9, fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.text(-0.10, 1.04, "(b)", transform=ax.transAxes, fontsize=15, fontweight="bold")
    fig.savefig(FIG_DIR / "Fig4_revised.png", bbox_inches="tight")
    plt.close(fig)


def make_fig5(processed_rows: list[dict], co60_rows: list[dict]) -> None:
    style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), constrained_layout=True)
    matrix = np.full((len(ORDER), 5), np.nan)
    sample_idx = {f"S{i}": i - 1 for i in range(1, 6)}
    nuclide_idx = {nuclide: i for i, nuclide in enumerate(ORDER)}
    for row in processed_rows:
        matrix[nuclide_idx[row["nuclide"]], sample_idx[row["sample"]]] = float(row["paired_change_percent"])
    ax = axes[0]
    image = ax.imshow(matrix, cmap="RdYlBu", vmin=-30, vmax=30, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            ax.text(j, i, "N/M" if np.isnan(value) else f"{value:+.1f}%", ha="center", va="center", color="#777777" if np.isnan(value) else "black", fontsize=9, fontweight="bold" if not np.isnan(value) else "normal")
    ax.set_xticks(range(5), [f"S{i}" for i in range(1, 6)])
    ax.set_yticks(range(4), [f"{n}\n({ENERGY_KEV[n]:.1f} keV)" for n in ORDER])
    ax.set_xlabel("Sample")
    cb = fig.colorbar(image, ax=ax, fraction=0.047, pad=0.03)
    cb.set_label("Paired change after solidification (%)")
    ax.text(-0.12, 1.04, "(a)", transform=ax.transAxes, fontsize=15, fontweight="bold")

    ax = axes[1]
    co = co60_rows[0]
    mean = co["mean_recovery_percent"]
    lower = co["ci95_lower_percent"]
    upper = co["ci95_upper_percent"]
    ax.errorbar([0], [mean], yerr=[[mean - lower], [upper - mean]], fmt="o", color=COLORS["60Co"], capsize=5, markersize=8)
    ax.axhspan(80, 120, color="#4AACB0", alpha=0.10, label="Paired recovery range (80-120%)")
    ax.axhline(100, color="#555555", linestyle="--", linewidth=1.1)
    ax.set_xlim(-0.7, 0.7)
    ax.set_ylim(75, 125)
    ax.set_xticks([0], [r"$^{60}$Co"])
    ax.set_ylabel("Paired recovery (%)")
    ax.text(0.06, mean, f"Mean {mean:.1f}%\n95% CI {lower:.1f}-{upper:.1f}%\n(n={co['n']})", va="center", fontsize=9)
    ax.legend(frameon=False, loc="lower center", fontsize=8)
    ax.text(-0.10, 1.04, "(b)", transform=ax.transAxes, fontsize=15, fontweight="bold")
    fig.savefig(FIG_DIR / "Fig5_revised.png", bbox_inches="tight")
    plt.close(fig)


def make_fig2() -> None:
    source = FIG_DIR / "Fig2.png"
    output = FIG_DIR / "Fig2_revised.png"
    image = Image.open(source).convert("RGB")
    draw = ImageDraw.Draw(image)
    width, height = image.size
    font = ImageFont.truetype("arial.ttf", 46) if Path("C:/Windows/Fonts/arial.ttf").exists() else ImageFont.load_default()
    x0, y0, x1, y1 = int(width * 0.255), int(height * 0.855), int(width * 0.495), int(height * 0.915)
    draw.rectangle((x0, y0, x1, y1), fill="white")
    draw.text((x0 + 8, y0 + 4), "CaCl2\u00b72H2O mass measurement", fill="black", font=font)
    image.save(output, dpi=(300, 300))


def write_validation_note(processed_rows: list[dict], co60_rows: list[dict]) -> None:
    failures = [row for row in processed_rows if row["reference_bias_within_plusminus20"] == "no"]
    note = [
        "# Calculation audit",
        "",
        "Generated from manuscript Table 1 values. Instrument-export raw files were not available.",
        "Before submission, replace or verify `manuscript_measurements.csv` against the original spectroscopy exports.",
        "",
        "## Key distinction",
        "",
        "- `recovery_percent_after_over_before`: paired preservation metric for the solidification step.",
        "- `reference_bias_percent_after_vs_KINS`: agreement of the post-solidification result with the certified KINS value.",
        "- These metrics are not interchangeable.",
        "",
        f"Reference-bias values outside +/-20%: {len(failures)}",
    ]
    note.extend(f"- {row['nuclide']} {row['sample']}: {row['reference_bias_percent_after_vs_KINS']}%" for row in failures)
    note.extend(
        [
            "",
            "## 60Co sensitivity analysis",
            "",
            f"- Observed mean recovery: {co60_rows[0]['mean_recovery_percent']}%",
            f"- 95% CI: {co60_rows[0]['ci95_lower_percent']}% to {co60_rows[0]['ci95_upper_percent']}%",
            f"- Minimum detectable absolute shift at 80% power: {co60_rows[0]['minimum_detectable_absolute_shift_percentage_points']} percentage points",
            "- This is a descriptive sensitivity analysis based on the observed SD, not prospective evidence that n=3 is sufficient for all uses.",
            "",
            "## Submission TODO",
            "",
            "- Replace `[DATA_REPOSITORY_URL_TO_BE_ADDED_BEFORE_SUBMISSION]` with a public GitHub or Zenodo URL.",
            "- Confirm whether Fig1 was created without generative AI. Elsevier does not permit generative-AI artwork in submitted manuscripts.",
        ]
    )
    (VAL_DIR / "calculation_audit.md").write_text("\n".join(note) + "\n", encoding="utf-8")


def main() -> None:
    for directory in (RAW_DIR, FIG_DIR, VAL_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    raw_rows, processed_rows = build_data()
    summary_rows = build_summary(processed_rows)
    co60_rows = build_co60_sensitivity(summary_rows)
    write_csv(RAW_DIR / "manuscript_measurements.csv", list(raw_rows[0]), raw_rows)
    write_csv(RAW_DIR / "processed_metrics.csv", list(processed_rows[0]), processed_rows)
    write_csv(RAW_DIR / "recovery_summary.csv", list(summary_rows[0]), summary_rows)
    write_csv(RAW_DIR / "co60_sensitivity_analysis.csv", list(co60_rows[0]), co60_rows)
    make_fig2()
    make_fig4(raw_rows, summary_rows)
    make_fig5(processed_rows, co60_rows)
    write_validation_note(processed_rows, co60_rows)
    print("Analysis package rebuilt successfully.")


if __name__ == "__main__":
    main()
