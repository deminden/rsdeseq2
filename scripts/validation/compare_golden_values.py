#!/usr/bin/env python3
"""Compare two golden value dumps field by field.

Each dump comes from `GOLDEN_DUMP=1 cargo golden --write --manifest PATH`.
For every field whose values differ, print the number of values, the largest
absolute and relative difference, and whether NaN or length patterns differ.

Usage: compare_golden_values.py BEFORE AFTER [--summary]
"""
import math
import struct
import sys


def load(path):
    out = {}
    for line in open(path):
        if line.startswith('#') or '\t' not in line:
            continue
        key, _, rest = line.rstrip('\n').rpartition('\t')
        if '|' not in rest:
            out[key] = None
            continue
        digest, values = rest.split('|', 1)
        out[key] = [struct.unpack('>d', bytes.fromhex(v))[0] for v in values.split(',') if v]
    return out


def main():
    before, after = load(sys.argv[1]), load(sys.argv[2])
    summary = '--summary' in sys.argv
    by_field = {}
    for key in sorted(set(before) | set(after)):
        a, b = before.get(key), after.get(key)
        if a == b:
            continue
        field = key.split('\t')[-1]
        workflow = key.split('\t')[1]
        name = f'{workflow}\t{field}'
        if a is None or b is None or len(a) != len(b):
            by_field.setdefault(name, ['shape', 0, 0.0, 0.0])
            continue
        entry = by_field.setdefault(name, ['values', 0, 0.0, 0.0])
        for x, y in zip(a, b):
            if math.isnan(x) and math.isnan(y):
                continue
            if math.isnan(x) != math.isnan(y) or math.isinf(x) or math.isinf(y):
                if x != y:
                    entry[0] = 'nan/inf'
                continue
            d = abs(x - y)
            if d > 0:
                entry[1] += 1
                entry[2] = max(entry[2], d)
                entry[3] = max(entry[3], d / max(abs(x), abs(y)))
    print('workflow\tfield\tkind\tvalues_changed\tmax_abs\tmax_rel')
    for name, (kind, n, da, dr) in sorted(by_field.items(), key=lambda kv: -kv[1][3]):
        print(f'{name}\t{kind}\t{n}\t{da:.3e}\t{dr:.3e}')
    if summary:
        print(f'# {len(by_field)} fields differ', file=sys.stderr)


if __name__ == '__main__':
    main()
