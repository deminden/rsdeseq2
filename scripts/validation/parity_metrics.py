#!/usr/bin/env python3
"""Significance-aware parity metrics shared by the validation scripts.

Absolute differences alone hide what matters downstream. These metrics add:

- log10 differences of p-values and adjusted p-values where the reference is
  below 1e-3 (both values positive);
- Spearman correlation of the Wald statistic;
- overlap of the top 100, 500 and 2,000 genes ranked by p-value;
- significance flips at padj 0.05 and 0.1 (a missing padj counts as not
  significant), which `require_no_flips` turns into a hard gate.
"""
from __future__ import annotations

import math
from typing import Iterable, Mapping

SMALL_P = 1e-3
TOP_K = (100, 500, 2000)
ALPHAS = (0.05, 0.1)


def parse(value: str | float | None) -> float:
    if value is None:
        return math.nan
    if isinstance(value, float):
        return value
    if value in ('', 'NA', 'NaN'):
        return math.nan
    return float(value)


def _ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float:
    pairs = [(a, b) for a, b in zip(x, y) if math.isfinite(a) and math.isfinite(b)]
    if len(pairs) < 2:
        return math.nan
    rx = _ranks([a for a, _ in pairs])
    ry = _ranks([b for _, b in pairs])
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else math.nan


def _top(rows: Mapping[str, Mapping[str, str]], k: int) -> set[str]:
    ranked = sorted(
        (parse(row.get('pvalue')), gene) for gene, row in rows.items()
        if math.isfinite(parse(row.get('pvalue'))))
    return {gene for _, gene in ranked[:k]}


def significance_metrics(got: Mapping[str, Mapping[str, str]],
                         expected: Mapping[str, Mapping[str, str]]) -> dict[str, float | int]:
    """Metrics over the genes present in `expected`; `got` must cover them."""
    missing = [gene for gene in expected if gene not in got]
    if missing:
        raise ValueError(f'{len(missing)} reference genes missing, e.g. {missing[0]}')
    out: dict[str, float | int] = {}
    for column in ('pvalue', 'padj'):
        worst = 0.0
        n = 0
        for gene, row in expected.items():
            e, g = parse(row.get(column)), parse(got[gene].get(column))
            if math.isfinite(e) and 0 < e < SMALL_P and math.isfinite(g) and g > 0:
                worst = max(worst, abs(math.log10(g) - math.log10(e)))
                n += 1
        out[f'{column}_small_max_log10'] = worst
        out[f'{column}_small_n'] = n
    genes = list(expected)
    out['stat_spearman'] = spearman([parse(got[g].get('stat')) for g in genes],
                                    [parse(expected[g].get('stat')) for g in genes])
    for k in TOP_K:
        reference = _top(expected, k)
        candidate = _top({g: got[g] for g in genes}, k)
        out[f'top{k}_overlap'] = len(reference & candidate) / len(reference) if reference else math.nan
    for alpha in ALPHAS:
        flips = 0
        for gene, row in expected.items():
            e, g = parse(row.get('padj')), parse(got[gene].get('padj'))
            if (math.isfinite(e) and e < alpha) != (math.isfinite(g) and g < alpha):
                flips += 1
        out[f'flips_padj_{alpha:g}'] = flips
    return out


def require_no_flips(metrics: Mapping[str, float | int], context: str) -> None:
    """Hard gate: any significance flip against the reference fails."""
    flips = {key: value for key, value in metrics.items() if key.startswith('flips_') and value}
    if flips:
        raise SystemExit(f'{context}: significance flips against the reference: {flips}')


def rows_by_gene(rows: Iterable[Mapping[str, str]]) -> dict[str, Mapping[str, str]]:
    out = {}
    for row in rows:
        if row['gene'] in out:
            raise ValueError(f"duplicate gene {row['gene']}")
        out[row['gene']] = row
    return out
