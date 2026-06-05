#!/usr/bin/env python3
"""
TMT Plex Log Analyser
=====================
Walks a directory tree, finds all .log files, and for each file extracts:
  - Plex determination from every calibration-method section (noCorrection,
    mean, median, Cubic Spline, …) via the "SIGNIFICANT MEDIANS with BH" array
  - Plex determination from the param-medic / crux output section
  - Ground-truth plex from the folder name (e.g. "TMT16-Plex(PXD…)" → 16plex)

Writes a human-readable report with per-file details and aggregate statistics.

Usage
-----
  python analyse_tmt_logs.py /path/to/root -o report.txt
  python analyse_tmt_logs.py /path/to/root --ignore-dirs crux-output alllogs
"""

from __future__ import annotations

import os
import re
import sys
import argparse
from datetime import datetime
from pathlib import Path
from collections import defaultdict
from typing import Optional

# ── Plex / channel configuration ──────────────────────────────────────────────
#
# The final SIGNIFICANT MEDIANS with BH array maps to this channel order
# (35 entries total):
#
#   Idx  0-17  : standard 18plex TMTpro channels (126 … 135N)
#   Idx  18-34 : additional 35plex channels
#
# 6plex channels are the ones marked "## 6" in the original channel list:
#   idx 0  126       idx 1  127N
#   idx 4  128C      idx 5  129N
#   idx 8  130C      idx 9  131N
#
# PLEX_CONFIGS must be ordered SMALLEST → LARGEST so that the subset-match
# loop returns the tightest fitting label.

PLEX_CONFIGS: list[tuple[str, frozenset[int]]] = [
    ("2plex",  frozenset([0, 2])),
    ("6plex",  frozenset([0, 1, 4, 5, 8, 9])),
    ("10plex", frozenset(range(10))),
    ("11plex", frozenset(range(11))),
    ("16plex", frozenset(range(16))),
    ("18plex", frozenset(range(18))),
    ("35plex", frozenset(range(35))),
]

N_CHANNELS = 35


def bool_array_to_plex(bool_array: list[bool]) -> str:
    """Return the tightest plex label that is a superset of all True indices."""
    true_idx = frozenset(i for i, v in enumerate(bool_array) if v)
    if not true_idx:
        return "none"
    for label, valid in PLEX_CONFIGS:      # smallest first
        if true_idx <= valid:              # subset check
            return label
    return f"unknown(>{N_CHANNELS}ch)"


# ── Folder-name plex extraction ───────────────────────────────────────────────

_FOLDER_TMT_RE  = re.compile(r"TMT[-_]?(\d+)", re.IGNORECASE)
_FOLDER_UNK_RE  = re.compile(r"TMT[\s_-]+UNKNOWN", re.IGNORECASE)


def plex_from_path(path: Path) -> str:
    """
    Walk all components of *path* (from deepest to shallowest) and return the
    plex label encoded in the first TMT-containing component found.

    Examples
    --------
    .../TMT16-Plex(PXD075635)/file.log  →  "16plex"
    .../TMT-16PLEX(PXD014750)/file.log  →  "16plex"
    .../TMT10-Plex(PXD067886)/file.log  →  "10plex"
    .../(TMT UNKNOWN)Thuebingen/...     →  "unknown"
    .../LFQ(PXD073864)/file.log         →  "not_found"
    """
    for part in reversed(path.parts):
        if _FOLDER_UNK_RE.search(part):
            return "unknown"
        m = _FOLDER_TMT_RE.search(part)
        if m:
            return f"{m.group(1)}plex"
    return "not_found"


# ── Log parsing ───────────────────────────────────────────────────────────────

_ANSI_RE = re.compile(r"(?:\x1b|\x033)?\[[\d;]+m")

# Method separator (three lines):
#   ############…
#   ##  ↑PREV_METHOD | NEXT_METHOD ↓ ####…
#   ############…
# Group 1 = name of section that just ended  (above the separator)
# Group 2 = name of section that is starting (below the separator)
_METHOD_SEP_RE = re.compile(
    r"#{40,}\n"
    r"##\s*↑(\S+)[^|]*\|\s*(.+?)\s*↓[^\n]*\n"
    r"#{40,}",
)

# Separator between our output and crux / param-medic output
_CRUX_SEP_RE = re.compile(
    r"#{40,}\n"
    r"#\s*↓CRUX OUTPUT[^\n]*\n"
    r"#{40,}",
    re.IGNORECASE,
)

# "SIGNIFICANT MEDIANS with BH" boolean array  (the BH-corrected final call)
_BH_ARRAY_RE = re.compile(
    r"SIGNIFICANT MEDIANS with BH[^\n]*\n\s*(\[.*?\])",
    re.DOTALL,
)

# Crux / param-medic TMT detection patterns
_CRUX_NO_TMT_RE    = re.compile(r"TMT:\s*no reporter ions detected",          re.I)
_CRUX_TMT_PRES_RE  = re.compile(r"TMT(\d+)\s+is\s+present",                   re.I)
_CRUX_NPLEX_DET_RE = re.compile(r"TMT[:\s]*(\d+)[-\s]*plex\s+reporter\s+ions\s+detected", re.I)


def _strip_ansi(text: str) -> str:
    return _ANSI_RE.sub("", text)


def _parse_bool_array(text: str) -> Optional[list[bool]]:
    """Extract the last 'SIGNIFICANT MEDIANS with BH' boolean array in *text*."""
    matches = list(_BH_ARRAY_RE.finditer(text))
    if not matches:
        return None
    tokens = re.findall(r"True|False", matches[-1].group(1))
    return [t == "True" for t in tokens] if tokens else None


def _parse_crux_plex(text: str) -> str:
    """
    Parse the crux/param-medic output section and return the detected plex label.

    Priority order:
      1. "TMT: no reporter ions detected"  → "none"
      2. "TMT10 is present"                → "10plex"   (specific reagent check)
      3. "TMT: N-plex reporter ions …"     → "Nplex"
      4. fallback                          → "unknown"
    """
    if _CRUX_NO_TMT_RE.search(text):
        return "none"
    m = _CRUX_TMT_PRES_RE.search(text)
    if m:
        return f"{m.group(1)}plex"
    m = _CRUX_NPLEX_DET_RE.search(text)
    if m:
        return f"{m.group(1)}plex"
    return "unknown"


def parse_log(path: Path) -> dict:
    """
    Parse one .log file.

    Returns a dict with:
        methods      : OrderedDict[method_name → plex_label]
        param_medic  : plex_label string, or None if no crux section found
        folder_plex  : plex_label derived from folder name
        error        : error message string if parsing failed badly, else None
    """
    try:
        raw  = path.read_text(errors="replace")
        text = _strip_ansi(raw)
    except OSError as exc:
        return {"methods": {}, "param_medic": None,
                "folder_plex": plex_from_path(path), "error": str(exc)}

    # ── separate crux section ────────────────────────────────────────────────
    crux_match = _CRUX_SEP_RE.search(text)
    if crux_match:
        our_text        = text[: crux_match.start()]
        crux_text       = text[crux_match.end() :]
        param_medic_plex = _parse_crux_plex(crux_text)
    else:
        our_text         = text
        param_medic_plex = None

    # ── split our_text into named method sections ────────────────────────────
    sep_matches = list(_METHOD_SEP_RE.finditer(our_text))
    sections: list[tuple[str, str]] = []     # (name, section_text)

    if not sep_matches:
        sections.append(("method_0", our_text))
    else:
        # Section before the first separator is named by group(1) of sep[0]
        sections.append((sep_matches[0].group(1),
                         our_text[: sep_matches[0].start()]))
        for i, sm in enumerate(sep_matches):
            start = sm.end()
            end   = sep_matches[i + 1].start() if i + 1 < len(sep_matches) else len(our_text)
            sections.append((sm.group(2), our_text[start:end]))

    # ── extract BH array from each section ──────────────────────────────────
    methods: dict[str, str] = {}
    for name, sec_text in sections:
        arr = _parse_bool_array(sec_text)
        if arr is not None:
            methods[name] = bool_array_to_plex(arr)
        else:
            methods[name] = "parse_error"

    return {
        "methods":     methods,
        "param_medic": param_medic_plex,
        "folder_plex": plex_from_path(path),
        "error":       None,
    }


# ── Directory walking ─────────────────────────────────────────────────────────

def walk_logs(root: Path, ignore_dirs: set[str]) -> list[Path]:
    """Collect all .log files under *root*, skipping directories in *ignore_dirs*."""
    found: list[Path] = []
    for dirpath_str, dirnames, filenames in os.walk(root):
        dirpath = Path(dirpath_str)
        # Prune in-place so os.walk does not descend
        dirnames[:] = [d for d in dirnames if d not in ignore_dirs]
        for fn in filenames:
            if fn.lower().endswith(".log"):
                found.append(dirpath / fn)
    return sorted(found)


# ── Report generation ─────────────────────────────────────────────────────────

_W = 80   # report width


def _hr(char: str = "=", w: int = _W) -> str:
    return char * w


def _hdr(text: str, char: str = "=", w: int = _W) -> str:
    return f" {text} ".center(w, char)


def _match_str(a: str, b: str) -> str:
    return "✓ MATCH   " if a == b else f"✗ MISMATCH  ({a} vs {b})"


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def generate_report(
    results: list[tuple[Path, dict]],
    root: Path,
    ignore_dirs: set[str],
) -> str:
    out: list[str] = []

    out.append(_hr())
    out.append("TMT Plex Analysis Report")
    out.append(f"Generated   : {datetime.now():%Y-%m-%d %H:%M:%S}")
    out.append(f"Root        : {root}")
    out.append(f"Log files   : {len(results)}")
    if ignore_dirs:
        out.append(f"Ignored dirs: {', '.join(sorted(ignore_dirs))}")
    out.append(_hr())
    out.append("")

    # ── per-file blocks ───────────────────────────────────────────────────────
    for path, res in results:
        out.append(_hr())
        out.append(f"FILE: {_rel(path, root)}")
        out.append(_hr("-"))

        folder_plex  = res["folder_plex"]
        param_medic  = res["param_medic"]
        methods      = res["methods"]
        error        = res.get("error")

        if error:
            out.append(f"  ERROR reading file: {error}")
            out.append("")
            continue

        out.append(f"  {'Ground truth (folder name):':<34}{folder_plex}")
        out.append("")

        # ── our methods ───────────────────────────────────────────────────────
        out.append("  Our Methods:")
        NW = max((len(n) for n in methods), default=12) + 2
        for mname, mplex in methods.items():
            vs = _match_str(mplex, folder_plex)
            out.append(f"    {mname:<{NW}}{mplex:<14}  {vs}")

        out.append("")

        # ── param-medic ───────────────────────────────────────────────────────
        if param_medic is not None:
            vs = _match_str(param_medic, folder_plex)
            out.append(f"  {'Param-Medic (crux):':<34}{param_medic:<14}  {vs}")
            out.append("")
            out.append("  Methods vs Param-Medic:")
            for mname, mplex in methods.items():
                vs = _match_str(mplex, param_medic)
                out.append(f"    {mname:<{NW}}{vs}")
        else:
            out.append("  Param-Medic (crux): (section not found in log)")

        out.append("")

    # ── statistics ────────────────────────────────────────────────────────────
    out.append(_hr())
    out.append(_hdr("STATISTICS"))
    out.append(_hr())
    out.append("")

    # Collect ordered list of unique method names (preserving first-seen order)
    method_names: list[str] = []
    _seen: set[str] = set()
    for _, res in results:
        for mname in res["methods"]:
            if mname not in _seen:
                method_names.append(mname)
                _seen.add(mname)

    n_total   = len(results)
    n_with_pm = sum(1 for _, r in results if r["param_medic"] is not None)
    n_errors  = sum(1 for _, r in results if r.get("error"))

    out.append(f"  Total log files analysed   : {n_total}")
    out.append(f"  Files with param-medic data: {n_with_pm}")
    out.append(f"  Files with read errors     : {n_errors}")
    out.append("")

    # ── accuracy vs folder name ───────────────────────────────────────────────
    out.append("  Accuracy vs folder name (ground truth):")
    cols = f"    {'Method':<24}{'Correct':>9}{'Total':>8}{'Accuracy':>11}"
    out.append(cols)
    out.append("    " + "-" * (len(cols) - 4))

    for mname in method_names:
        pairs   = [(r["methods"].get(mname), r["folder_plex"])
                   for _, r in results
                   if mname in r["methods"] and not r.get("error")]
        correct = sum(1 for mp, fp in pairs if mp == fp)
        total   = len(pairs)
        pct     = 100 * correct / total if total else 0.0
        out.append(f"    {mname:<24}{correct:>9}{total:>8}{pct:>10.1f}%")

    if n_with_pm:
        pm_correct = sum(1 for _, r in results
                         if r["param_medic"] is not None
                         and r["param_medic"] == r["folder_plex"])
        pct = 100 * pm_correct / n_with_pm
        out.append(f"    {'param-medic':<24}{pm_correct:>9}{n_with_pm:>8}{pct:>10.1f}%")

    out.append("")

    # ── plex distribution ─────────────────────────────────────────────────────
    out.append("  Plex distribution across all files:")

    col_names = ["folder"] + method_names
    if n_with_pm:
        col_names.append("param-medic")

    dists: dict[str, dict[str, int]] = {c: defaultdict(int) for c in col_names}
    for _, r in results:
        if r.get("error"):
            continue
        dists["folder"][r["folder_plex"]] += 1
        for mname in method_names:
            if mname in r["methods"]:
                dists[mname][r["methods"][mname]] += 1
        if r["param_medic"] is not None:
            dists["param-medic"][r["param_medic"]] += 1

    all_labels_seen: set[str] = set()
    for d in dists.values():
        all_labels_seen |= d.keys()

    _ORDERED_LABELS = [
        "none", "2plex", "6plex", "10plex", "11plex",
        "16plex", "18plex", "35plex", "unknown", "not_found", "parse_error",
    ]
    label_order = [l for l in _ORDERED_LABELS if l in all_labels_seen] + \
                  sorted(all_labels_seen - set(_ORDERED_LABELS))

    CW = 14  # column width
    hrow = f"    {'Label':<{CW}}" + "".join(f"{c:>{CW}}" for c in col_names)
    out.append(hrow)
    out.append("    " + "-" * (len(hrow) - 4))
    for lbl in label_order:
        row = f"    {lbl:<{CW}}" + "".join(
            f"{dists[c].get(lbl, 0):>{CW}}" for c in col_names
        )
        out.append(row)

    # ── method agreement with param-medic ─────────────────────────────────────
    if n_with_pm:
        out.append("")
        out.append("  Method agreement with param-medic:")
        for mname in method_names:
            pairs   = [(r["methods"].get(mname), r["param_medic"])
                       for _, r in results
                       if mname in r["methods"]
                       and r["param_medic"] is not None
                       and not r.get("error")]
            correct = sum(1 for mp, pm in pairs if mp == pm)
            total   = len(pairs)
            pct     = 100 * correct / total if total else 0.0
            out.append(f"    {mname:<24}{correct:>4} / {total:<5} ({pct:.1f}%)")

    out.append("")
    out.append(_hr())
    return "\n".join(out)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Walk a directory tree and analyse TMT plex .log files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument(
        "root",
        type=Path,
        help="Root directory to walk",
    )
    ap.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("tmt_plex_report.txt"),
        help="Output report file (default: tmt_plex_report.txt)",
    )
    ap.add_argument(
        "--ignore-dirs",
        nargs="*",
        default=[],
        metavar="DIR",
        help=(
            "Directory *names* (not paths) to skip during walking. "
            "Example: --ignore-dirs crux-output alllogs"
        ),
    )
    args = ap.parse_args(argv)

    root = args.root.resolve()
    if not root.is_dir():
        print(f"Error: '{root}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    ignore = set(args.ignore_dirs)

    print(f"Walking '{root}' …")
    if ignore:
        print(f"  (ignoring directories named: {', '.join(sorted(ignore))})")

    log_files = walk_logs(root, ignore)
    print(f"Found {len(log_files)} .log file(s).")

    if not log_files:
        print("Nothing to analyse.")
        sys.exit(0)

    results: list[tuple[Path, dict]] = []
    max_rel_len = max(len(_rel(f, root)) for f in log_files)

    for lf in log_files:
        rel = _rel(lf, root)
        print(f"  Parsing {rel:<{max_rel_len}} … ", end="", flush=True)
        res = parse_log(lf)
        results.append((lf, res))
        if res.get("error"):
            print(f"ERROR: {res['error']}")
        else:
            method_calls = ", ".join(
                f"{k}={v}" for k, v in res["methods"].items()
            )
            pm = res["param_medic"] or "—"
            folder = res["folder_plex"]
            print(f"OK  folder={folder}  pm={pm}  [{method_calls}]")

    report = generate_report(results, root, ignore)

    args.output.write_text(report, encoding="utf-8")
    print(f"\nReport written to '{args.output}'")


if __name__ == "__main__":
    main()
