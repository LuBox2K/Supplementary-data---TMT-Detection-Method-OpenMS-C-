import re
from pathlib import Path
from collections import defaultdict



## ── Plex definitions ────────────────────────────────────────────────────────
PLEX_CONFIGS = [
    ('35plex', list(range(35))),
    ('18plex', list(range(18))),
    ('16plex', list(range(16))),
    ('11plex', list(range(11))),
    ('10plex', list(range(10))),
    ('6plex', [0, 1, 4, 5, 8, 9]),
    ('2plex', [0, 2]),
]


def bool_array_to_plex(bool_array):
    true_indices = {i for i, v in enumerate(bool_array) if v}
    if not true_indices:
        return 'none'
    for label, indices in reversed(PLEX_CONFIGS):
        if true_indices.issubset(set(indices)):
            return label
    return 'unknown'


## ── Parsing ─────────────────────────────────────────────────────────────────

## Production separator:
##   ############################################################
##   # ↓CRUX OUTPUT   |   OUR OUTPUT↑
##   ############################################################
SEPARATOR_RE = re.compile(r'#{10,}\n#[^\n]*\n#{10,}')


def split_sections(text):
    match = SEPARATOR_RE.search(text)
    if match:
        return text[:match.start()].strip(), text[match.end():].strip()
    return '', text.strip()


def parse_python_section(text):
    """Extract 'SIGNIFICANT MEDIANS with BH' boolean array."""
    match = re.search(r'SIGNIFICANT MEDIANS with BH.*?\n\s*(\[.*?\])', text, re.DOTALL)
    if not match:
        return None
    tokens = re.findall(r'True|False', match.group(1))
    return [t == 'True' for t in tokens] if tokens else None


def parse_crux_section(text):
    if 'TMT: no reporter ions detected' in text:
        return 'none'
    if 'TMT10 is present' in text:
        return '10plex'
    if 'TMT: 6-plex reporter ions detected' in text:
        return '6plex'
    return 'unknown'


def parse_file(filepath):
    text = filepath.read_text(errors='replace')
    python_text, crux_text = split_sections(text)
    crux_call  = parse_crux_section(crux_text)
    bool_array = parse_python_section(python_text)
    python_call = bool_array_to_plex(bool_array) if bool_array is not None else 'no_array'
    return python_call, crux_call


## ── Main ────────────────────────────────────────────────────────────────────

def main():
    directory = Path(__file__).parent
    all_files = sorted(list(directory.glob('*.txt')) + list(directory.glob('*.log')))
    files = [f for f in all_files
             if 'param-medic' in f.read_text(errors='replace').lower()
             or 'crux' in f.read_text(errors='replace').lower()]

    if not files:
        print(f'No param-medic log files found in {directory}')
        return

    python_counts = defaultdict(int)
    crux_counts   = defaultdict(int)
    agreement     = 0
    results       = []

    for f in files:
        py_call, crux_call = parse_file(f)
        python_counts[py_call] += 1
        crux_counts[crux_call] += 1
        if py_call == crux_call:
            agreement += 1
        results.append((f.name, py_call, crux_call))

    col_w = 55
    print(f'\n{"File":<{col_w}} {"Python call":<14} {"Crux call"}')
    print('-' * (col_w + 30))
    for name, py_call, crux_call in results:
        marker = '' if py_call == crux_call else '  <- MISMATCH'
        print(f'{name:<{col_w}} {py_call:<14} {crux_call}{marker}')

    all_labels = sorted(set(list(python_counts) + list(crux_counts)))
    print('\n' + '=' * 60)
    print('SUMMARY')
    print('=' * 60)
    print(f'\n  {"Plex call":<16} {"Python":>8} {"Crux":>8}')
    print(f'  {"-"*16} {"-"*8} {"-"*8}')
    for label in all_labels:
        print(f'  {label:<16} {python_counts[label]:>8} {crux_counts[label]:>8}')
    print(f'\n  Total files:   {len(files):>4}')
    print(f'  Agreement:     {agreement:>4} / {len(files)}  ({100*agreement/len(files):.1f}%)')


if __name__ == '__main__':
    main()
