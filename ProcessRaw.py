#!/usr/bin/env python3
"""
process_raw.py – traverse a folder tree, optionally convert .raw → .mzML,
                 then run analyses and write log files.

Default (no flags): analyze all .mzML files and write logs.

Flags:
    --convert-raw       Convert .raw → .mzML (skips if .mzML already exists)
    --reconvert         Re-convert even if .mzML already exists (implies --convert-raw)
    --delete-converted  Delete .raw after a successful conversion
    --no-analyze        Skip the analysis step (findTMTChannels + param-medic)
    --dry-run           Print what would be done without doing anything
"""

import os
import sys
import glob
import argparse
import subprocess

# ── colour helper ─────────────────────────────────────────────────────────────
def col(inp, r=0, g=0, b=0):
    return f'\033[38;2;{r};{g};{b}m{inp}\033[0m'

# ── config ────────────────────────────────────────────────────────────────────
FILE_CONVERTER  = os.path.expanduser('~/openms-development/openms_build/bin/FileConverter')
THERMO_EXE      = os.path.expanduser('~/WindowsApps/ThermoFileParserWin/ThermoRawFileParser.exe')
FIND_TMT_SCRIPT = os.path.join(os.path.dirname(__file__), 'findTMTChannels.py')
CRUX_BIN        = os.path.expanduser('~/openms-development/crux-py/crux-toolkit/src/crux')
VENV_PYTHON     = os.path.expanduser('~/openms-development/py/.venv/bin/python3.14')   # ← adjust to your venv
DELIMITER       = '#'*60 + '\n# ↓CRUX OUTPUT   |   OUR OUTPUT↑  \n' + '#'*60 + '\n'

# ── conversion ────────────────────────────────────────────────────────────────
def convert_raw(raw_path: str, reconvert: bool, dry_run: bool) -> bool:
    """Convert a single .raw file to .mzML. Returns True on successful conversion."""
    mzml_path = os.path.splitext(raw_path)[0] + '.mzML'

    if os.path.exists(mzml_path) and not reconvert:
        print(f"{col('skip (exists):', r=180, g=180, b=0)} {mzml_path}")
        return False

    if dry_run:
        print(f"{col('[dry-run] would convert:', g=200, b=255)} {raw_path}")
        if delete_converted:
            print(f"{col('[dry-run] would delete:', r=255, g=100, b=100)} {raw_path}")
        return False

    command = [
        FILE_CONVERTER,
        '-in',  raw_path,
        '-out', mzml_path,
        '-RawToMzML:ThermoRaw_executable', THERMO_EXE,
    ]
    try:
        subprocess.run(command, check=True)
        print(f"{col('converted:', g=255, b=255)} {mzml_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"{col('Error converting', r=255)} {raw_path}: {e}", file=sys.stderr)
        return False

# ── analysis ──────────────────────────────────────────────────────────────────
def analyze_mzml(mzml_path: str, dry_run: bool):
    """Run findTMTChannels and param-medic; write output to <basename>.log."""
    log_path = os.path.splitext(mzml_path)[0] + '.log'

    if dry_run:
        print(f"{col('[dry-run] would analyze:', g=200, b=255)} {mzml_path}  →  {log_path}")
        return

    # ── Program 1: findTMTChannels (venv python, stdout → log) ───────────────
    cmd1 = [VENV_PYTHON, FIND_TMT_SCRIPT, mzml_path]
    print(f"{col('analyzing [1/2]:', g=180, b=255)} {mzml_path}")
    try:
        result1 = subprocess.run(
            cmd1,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        with open(log_path, 'w') as log:
            log.write(result1.stdout)
        print(f"{col('  → log written:', g=180)} {log_path}")
    except subprocess.CalledProcessError as e:
        print(f"{col('Error [findTMTChannels]', r=255)} {mzml_path}: {e}", file=sys.stderr)
        return

    # ── delimiter ─────────────────────────────────────────────────────────────
    with open(log_path, 'a') as log:
        log.write(DELIMITER)

    # ── Program 2: crux param-medic (stdout appended to log) ─────────────────
    cmd2 = [CRUX_BIN, 'param-medic', '--overwrite', 'T', mzml_path]
    print(f"{col('analyzing [2/2]:', g=180, b=255)} {mzml_path}")
    try:
        result2 = subprocess.run(
            cmd2,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        with open(log_path, 'a') as log:
            log.write(result2.stderr)
        print(f"{col('  → log appended:', g=180)} {log_path}")
    except subprocess.CalledProcessError as e:
        print(f"{col('Error [param-medic]', r=255)} {mzml_path}: {e}", file=sys.stderr)

# ── folder processing ─────────────────────────────────────────────────────────
def process_folder(folder: str, args):
    # conversion pass (only if requested)
    if args.convert_raw or args.reconvert:
        for raw_path in sorted(glob.glob(os.path.join(folder, '*.raw'))):
            convert_raw(
                raw_path,
                reconvert=args.reconvert,
                dry_run=args.dry_run,
            )

    # analysis pass (default behaviour, skip with --no-analyze)
    if not args.no_analyze:
        for mzml_path in sorted(glob.glob(os.path.join(folder, '*.mzML'))):
            analyze_mzml(mzml_path, dry_run=args.dry_run)

    if args.delete_raws:
        for raw_path in sorted(glob.glob(os.path.join(folder, '*.raw'))):
            if not args.dry_run:
                os.remove(raw_path)
            print(f"{col('deleted:', r=255, g=100, b=100)} {raw_path}")

# ── entry point ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description=(
            'Traverse folders and analyze .mzML files (default). '
            'Optionally convert .raw → .mzML first.'
        )
    )
    parser.add_argument('root', help='Root directory to traverse')
    parser.add_argument('--convert-raw', action='store_true',
                        help='Convert .raw → .mzML (skip if .mzML already exists)')
    parser.add_argument('--reconvert', action='store_true',
                        help='Re-convert even if .mzML already exists (implies --convert-raw)')
    parser.add_argument('--delete-raws', action='store_true',
                        help='Delete .raw after a successful conversion')
    parser.add_argument('--no-analyze', action='store_true',
                        help='Skip the analysis step (findTMTChannels + param-medic)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print actions without executing them')
    args = parser.parse_args()

    root = os.path.expanduser(args.root)
    if not os.path.isdir(root):
        print(f"{col('Error:', r=255)} '{root}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    print(f"{col('Root:', g=200, b=200)} {root}")
    if args.dry_run:
        print(col('[DRY RUN – nothing will be changed]', r=255, g=200, b=0))

    for dirpath, dirnames, _ in os.walk(root):
        dirnames.sort()
        print(f"\n{col('entering:', b=200, g=200)} {dirpath}")
        process_folder(dirpath, args)

    print(f"\n{col('done.', g=255)}")


if __name__ == '__main__':
    main()
