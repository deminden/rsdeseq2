#!/usr/bin/env python3
"""Structural metrics tracked by docs/reports/hardening-progress.md.

All counts are static and dependency-free, so they are reproducible from any
checkout. They are trend indicators, not proofs: the parser strips comments,
string literals and `#[cfg(test)]` items with a brace matcher rather than a
full Rust parser.

Metrics:
  production_unwrap / _expect / _unreachable  `.unwrap()`, `.expect(` and
      `unreachable!` in crate sources outside comments, strings and
      `#[cfg(test)]` items.
  include_sites        `include!(` invocations in crate sources.
  pub_fn               `pub fn` declarations (any `pub(...)` restriction excluded).
  builder_fit_methods  `pub fn fit*` methods in `core/implementation`.
  mutual_module_pairs  unordered top-level module pairs that reference each
      other through `crate::` paths or `use super::`-free imports.
  integration_test_binaries  Rust files directly under `tests/`.
  platform_cfg_sites   `all(target_os = "linux", target_arch = "x86_64")` in
      `cfg`, `cfg!` and `not(...)` forms.
  clone_share          share of production lines inside a repeated window of
      CLONE_WINDOW consecutive normalized non-trivial lines.
  r_entry_points / r_unwrapped_entry_points  `extern "C" fn rsdeseq2_r_*`
      definitions in the R bridge, and those whose body does not call
      `run_backend` or `guard_value`.
  r_unsafe_uses        `unsafe` tokens in the R bridge sources.
"""
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CRATE = ROOT / 'crates/rsdeseq2'
BRIDGE = ROOT / 'r-pkg/rsdeseq2/src/rust/src'
CLONE_WINDOW = 10
PLATFORM_CFG = re.compile(r'all\(\s*target_os\s*=\s*"linux"\s*,\s*target_arch\s*=\s*"x86_64"\s*\)')
CRATE_REF = re.compile(r'\bcrate::([a-z_][a-z0-9_]*)')


def strip_comments_and_strings(text: str) -> str:
    """Blank out comments, string and char literals, keeping line structure."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith('//', i):
            j = text.find('\n', i)
            j = n if j < 0 else j
            i = j
        elif text.startswith('/*', i):
            depth, j = 1, i + 2
            while j < n and depth:
                if text.startswith('/*', j):
                    depth, j = depth + 1, j + 2
                elif text.startswith('*/', j):
                    depth, j = depth - 1, j + 2
                else:
                    j += 1
            out.append('\n' * text.count('\n', i, j))
            i = j
        elif c == 'r' and re.match(r'r#*"', text[i:i + 10]) and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == '_')):
            hashes = len(re.match(r'r(#*)"', text[i:]).group(1))
            end = text.find('"' + '#' * hashes, i + 2 + hashes)
            end = n if end < 0 else end + 1 + hashes
            out.append('""' + '\n' * text.count('\n', i, end))
            i = end
        elif c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == '\\' else 1
            out.append('""' + '\n' * text.count('\n', i, j))
            i = j + 1
        elif c == "'" and re.match(r"'(\\.|[^\\'])'", text[i:i + 4]):
            m = re.match(r"'(\\.|[^\\'])'", text[i:i + 4])
            out.append("' '")
            i += len(m.group(0))
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def strip_line_comments(text: str) -> str:
    """Drop `//` comment lines (sufficient for attribute/`cfg!` counting)."""
    return '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('//'))


def strip_cfg_test(text: str) -> str:
    """Remove items annotated `#[cfg(test)]` (modules, functions, impls)."""
    result = []
    i = 0
    marker = re.compile(r'#\[cfg\(test\)\]')
    while True:
        m = marker.search(text, i)
        if not m:
            result.append(text[i:])
            break
        result.append(text[i:m.start()])
        j = m.end()
        brace = text.find('{', j)
        semi = text.find(';', j)
        if brace < 0 or (0 <= semi < brace):
            i = semi + 1 if semi >= 0 else len(text)
            continue
        depth, k = 1, brace + 1
        while k < len(text) and depth:
            depth += {'{': 1, '}': -1}.get(text[k], 0)
            k += 1
        result.append('\n' * text.count('\n', m.start(), k))
        i = k
    return ''.join(result)


def production_sources() -> dict[Path, str]:
    sources = {}
    for path in sorted((CRATE / 'src').rglob('*.rs')):
        text = strip_cfg_test(strip_comments_and_strings(path.read_text()))
        if path.name == 'tests.rs' or '/tests/' in path.as_posix():
            text = ''
        sources[path] = text
    return sources


def top_module(path: Path) -> str:
    rel = path.relative_to(CRATE / 'src').parts
    name = rel[0]
    return name[:-3] if name.endswith('.rs') else name


def mutual_pairs(sources: dict[Path, str]) -> list[tuple[str, str]]:
    modules = {top_module(p) for p in sources} - {'lib', 'main'}
    edges = defaultdict(set)
    for path, text in sources.items():
        source = top_module(path)
        for target in CRATE_REF.findall(text):
            if target in modules and target != source:
                edges[source].add(target)
    pairs = {tuple(sorted((a, b))) for a in edges for b in edges[a] if a in edges.get(b, set())}
    return sorted(pairs)


def clone_share(sources: dict[Path, str]) -> float:
    lines = []
    for path, text in sources.items():
        for number, line in enumerate(text.splitlines()):
            norm = re.sub(r'\s+', ' ', line).strip()
            if len(norm) > 3 and norm not in ('}', '};', '{', ')', '),', '});', ']', '],'):
                lines.append((path, number, norm))
    windows = defaultdict(list)
    for start in range(len(lines) - CLONE_WINDOW + 1):
        chunk = lines[start:start + CLONE_WINDOW]
        if len({p for p, _, _ in chunk}) != 1:
            continue
        windows[tuple(n for _, _, n in chunk)].append(start)
    cloned = set()
    for starts in windows.values():
        if len(starts) > 1:
            for s in starts:
                cloned.update(range(s, s + CLONE_WINDOW))
    return len(cloned) / len(lines) if lines else 0.0


def bridge_metrics() -> dict[str, int]:
    # String literals are blanked, so the ABI string appears as `""`.
    entry = re.compile(r'extern\s+""\s+fn\s+(rsdeseq2_r_[a-z0-9_]+)\s*\(')
    entries, unwrapped, unsafe = 0, 0, 0
    for path in sorted(BRIDGE.rglob('*.rs')):
        text = strip_comments_and_strings(path.read_text())
        unsafe += len(re.findall(r'\bunsafe\b', text))
        for m in entry.finditer(text):
            entries += 1
            brace = text.find('{', m.end())
            depth, k = 1, brace + 1
            while k < len(text) and depth:
                depth += {'{': 1, '}': -1}.get(text[k], 0)
                k += 1
            if 'run_backend' not in text[brace:k] and 'guard_value' not in text[brace:k]:
                unwrapped += 1
    return dict(r_entry_points=entries, r_unwrapped_entry_points=unwrapped, r_unsafe_uses=unsafe)


def collect() -> dict:
    sources = production_sources()
    raw = {p: p.read_text() for p in sorted((CRATE / 'src').rglob('*.rs'))}
    all_code = '\n'.join(sources.values())
    stripped_all = {p: strip_comments_and_strings(t) for p, t in raw.items()}
    builder_dir = CRATE / 'src/core/implementation'
    metrics = dict(
        production_unwrap=len(re.findall(r'\.unwrap\(\)', all_code)),
        production_expect=len(re.findall(r'\.expect\(', all_code)),
        production_unreachable=len(re.findall(r'\bunreachable!', all_code)),
        include_sites=sum(len(re.findall(r'\binclude!\s*\(', t)) for t in stripped_all.values()),
        pub_fn=sum(len(re.findall(r'\bpub\s+(?:const\s+)?(?:unsafe\s+)?fn\b', t)) for t in sources.values()),
        builder_fit_methods=sum(len(re.findall(r'\bpub\s+fn\s+fit(?!_type\b)[a-z0-9_]*\s*[<(]', t))
                                for p, t in sources.items() if builder_dir in p.parents),
        integration_test_binaries=len(list((CRATE / 'tests').glob('*.rs'))),
        platform_cfg_sites=sum(len(PLATFORM_CFG.findall(strip_line_comments(t))) for t in raw.values()),
        production_lines=sum(len([l for l in t.splitlines() if l.strip()]) for t in sources.values()),
    )
    metrics['production_panics'] = (metrics['production_unwrap'] + metrics['production_expect']
                                    + metrics['production_unreachable'])
    pairs = mutual_pairs(sources)
    metrics['mutual_module_pairs'] = len(pairs)
    metrics['mutual_module_pair_list'] = [f'{a}<->{b}' for a, b in pairs]
    metrics['clone_share_percent'] = round(100 * clone_share(sources), 2)
    metrics.update(bridge_metrics())
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--json', action='store_true', help='print JSON instead of a table')
    args = parser.parse_args()
    metrics = collect()
    if args.json:
        print(json.dumps(metrics, indent=2))
    else:
        for key, value in metrics.items():
            if isinstance(value, list):
                value = ', '.join(value)
            print(f'{key}\t{value}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
