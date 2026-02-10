#!/usr/bin/env python3
r"""
check_bib_usage.py

Scan a .bib file for entry keys and search the workspace for usages (\cite, \nocite, etc.).
Usage: python scripts/check_bib_usage.py --bib bibliography/literature.bib
"""
from pathlib import Path
import re
import argparse
import sys


def extract_bib_keys(bib_path: Path):
    text = bib_path.read_text(encoding="utf8")
    # Capture entry blocks: @type{key, ...}
    entry_pat = re.compile(r"@(\w+)\s*{\s*([^,\n\r]+)\s*,(.*?)(?=\n@|\Z)", re.IGNORECASE | re.DOTALL)
    entries = []
    for m in entry_pat.finditer(text):
        etype = m.group(1).strip()
        key = m.group(2).strip()
        body = m.group(3).strip()
        # try to extract title and doi
        title = None
        doi = None
        t = re.search(r"title\s*=\s*[{\"](.*?)[}\"]\s*,?\n", body, re.IGNORECASE | re.DOTALL)
        if t:
            title = t.group(1).strip()
        d = re.search(r"doi\s*=\s*[{\"](.*?)[}\"]\s*,?\n", body, re.IGNORECASE | re.DOTALL)
        if d:
            doi = d.group(1).strip()
        entries.append({'type': etype, 'key': key, 'title': title, 'doi': doi, 'body': body})
    return entries


def find_files(root: Path, exts=None):
    if exts is None:
        exts = {'.tex', '.md', '.Rmd', '.sty'}
    files = [p for p in root.rglob('*') if p.suffix.lower() in exts]
    return files


def scan_for_usage(keys, files, bib_path: Path = None):
    used = {k: 0 for k in keys}

    # Precompile patterns for faster search
    key_patterns = {k: re.compile(re.escape(k)) for k in keys}
    missing_cites = {}

    for f in files:
        # skip the bibliography file itself
        try:
            if bib_path and f.resolve() == bib_path.resolve():
                continue
        except Exception:
            pass

        try:
            txt = f.read_text(encoding='utf8', errors='ignore')
        except Exception:
            continue

        # Special handling for \nocite{*} or \nocite{key1,key2}
        for nocite in re.findall(r'\\nocite\s*{([^}]*)}', txt):
            content = nocite.strip()
            if content == '*':
                for k in keys:
                    used[k] = max(used[k], 1)
            else:
                for part in content.split(','):
                    p = part.strip()
                    if p in used:
                        used[p] = max(used[p], 1)

        # Find citation commands and record missing keys (per-line)
        for lineno, line in enumerate(txt.splitlines(), start=1):
            for cite_match in re.finditer(r"\\(cite|citep|citet|parencite|autocite|footcite)\s*{([^}]*)}", line):
                group = cite_match.group(2)
                for part in group.split(','):
                    p = part.strip()
                    if p and p not in used:
                        missing_cites.setdefault(p, []).append((str(f), lineno, line.strip()))

        # Generic search for key occurrences
        for k, pat in key_patterns.items():
            if pat.search(txt):
                used[k] += 1

    return used, missing_cites


def main():
    ap = argparse.ArgumentParser(description='Check .bib keys for usage in workspace files')
    ap.add_argument('--bib', type=Path, default=Path('bibliography/literature.bib'))
    ap.add_argument('--root', type=Path, default=Path('.'), help='Workspace root to search')
    ap.add_argument('--ext', action='append', help='Additional file extensions to search (e.g. .sty)')
    args = ap.parse_args()

    bib_path = args.bib
    if not bib_path.exists():
        print(f'Error: bib file not found: {bib_path}', file=sys.stderr)
        sys.exit(2)

    entries = extract_bib_keys(bib_path)
    if not entries:
        print('No entries found in', bib_path)
        sys.exit(0)

    keys = [e['key'] for e in entries]

    # detect duplicate keys
    dup_keys = [k for k, cnt in __import__('collections').Counter(keys).items() if cnt > 1]

    # detect duplicate titles and DOIs
    title_map = {}
    doi_map = {}
    for e in entries:
        if e.get('title'):
            t = re.sub(r"\s+", " ", e['title'].strip().lower())
            title_map.setdefault(t, []).append(e['key'])
        if e.get('doi'):
            d = e['doi'].strip().lower()
            doi_map.setdefault(d, []).append(e['key'])

    dup_titles = {t: ks for t, ks in title_map.items() if len(ks) > 1}
    dup_dois = {d: ks for d, ks in doi_map.items() if len(ks) > 1}

    exts = None
    if args.ext:
        exts = {e if e.startswith('.') else f'.{e}' for e in args.ext}
    files = find_files(args.root, exts=exts)

    used, missing_cites = scan_for_usage(keys, files, bib_path=bib_path)

    total = len(keys)
    used_keys = [k for k, c in used.items() if c > 0]
    unused_keys = [k for k, c in used.items() if c == 0]

    print(f'Total bib entries: {total}')
    print(f'Used entries: {len(used_keys)}')
    print(f'Unused entries: {len(unused_keys)}')

    if unused_keys:
        print('\nList of unused keys:')
        for k in unused_keys:
            print('-', k)
    else:
        print('\nAll bib entries appear to be referenced (or covered by \nocite{*}).')

    if dup_keys:
        print('\nDuplicate bib keys in file:')
        for k in dup_keys:
            print('-', k)

    if dup_titles:
        print('\nDuplicate titles (different keys with same title):')
        for t, ks in dup_titles.items():
            print('-', ks)

    if dup_dois:
        print('\nDuplicate DOIs (different keys with same DOI):')
        for d, ks in dup_dois.items():
            print('-', ks)

    if missing_cites:
        print('\nCitation keys referenced in source but missing from bibliography:')
        for key, locs in missing_cites.items():
            print(f'- {key}')
            for fn, ln, line in locs:
                print(f'    {fn}:{ln}: {line}')


if __name__ == '__main__':
    main()
