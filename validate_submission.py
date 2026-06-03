from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "main.tex"


def report(level: str, message: str) -> None:
    print(f"[{level}] {message}")


def main() -> int:
    text = MAIN.read_text(encoding="utf-8")
    errors = 0
    warnings = 0

    citations = set()
    for group in re.findall(r"\\cite\{([^}]+)\}", text):
        citations.update(item.strip() for item in group.split(","))
    bibitems = set(re.findall(r"\\bibitem\{([^}]+)\}", text))
    missing = sorted(citations - bibitems)
    unused = sorted(bibitems - citations)
    if missing:
        errors += 1
        report("ERROR", f"Missing bibliography items: {missing}")
    else:
        report("OK", f"All {len(citations)} cited keys have bibliography items.")
    if unused:
        warnings += 1
        report("WARN", f"Unused bibliography items: {unused}")

    figures = re.findall(r"\\includegraphics(?:\[[^\]]+\])?\{([^}]+)\}", text)
    missing_figures = [name for name in figures if not (ROOT / "figure" / name).exists()]
    if missing_figures:
        errors += 1
        report("ERROR", f"Missing figure files: {missing_figures}")
    else:
        report("OK", f"All {len(figures)} manuscript figures exist.")

    highlights = (ROOT / "highlight.tex").read_text(encoding="utf-8")
    items = re.findall(r"\\item\s+(.+)", highlights)
    too_long = [item for item in items if len(item) > 85]
    if too_long:
        errors += 1
        report("ERROR", f"Highlights longer than 85 characters: {too_long}")
    else:
        report("OK", f"All {len(items)} highlights are 85 characters or fewer.")

    repository_placeholder = "[DATA_REPOSITORY_URL_TO_BE_ADDED_BEFORE_SUBMISSION]"
    if repository_placeholder in text:
        warnings += 1
        report("WARN", "Replace the data-repository URL placeholder before submission.")

    if "Fig1.png" in figures:
        warnings += 1
        report("WARN", "Confirm Fig1 provenance: submitted Elsevier artwork must not use generative AI.")

    if text.count("{") != text.count("}"):
        errors += 1
        report("ERROR", "Brace counts differ in main.tex.")
    else:
        report("OK", "Brace counts match in main.tex.")

    report("SUMMARY", f"errors={errors}, warnings={warnings}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
