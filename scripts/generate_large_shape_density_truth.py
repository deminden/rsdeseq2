#!/usr/bin/env python3
"""Independent 100-decimal density truths from elementary log-gamma identities."""
import argparse
import hashlib
import json
import struct
from pathlib import Path

import mpmath as mp


def number(value):
    return struct.unpack('>d', bytes.fromhex(value.removeprefix('bits:')))[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('fresh_output', type=Path)
    args = parser.parse_args()
    if args.fresh_output.exists():
        parser.error('use a fresh output path')
    mp.mp.dps = 100
    source = json.loads(args.reference.read_text())
    rows = []
    for case in source['kernel_cases']:
        model = case['builtin_likelihood']
        kind = model['type']
        values = []
        for encoded in case['arguments']:
            x = mp.mpf(number(encoded))
            if kind == 'pois':
                y = number(model['y'][0])
                rate = mp.mpf(number(model['scale'][0]))
                log_link = model['link'] == 'log'
                shape = mp.mpf(y + 1e-5 if log_link else y + 1.0)
                argument = mp.exp(-x) if log_link else abs(x)
                density = shape * mp.log(rate) + (shape - 1) * mp.log(argument) - rate * argument - mp.loggamma(shape)
                value = density - mp.log(shape if log_link else rate)
            elif kind == 'binom':
                y = number(model['y'][0]); trials = number(model['trials'][0])
                logit = model['link'] == 'logit'
                a = mp.mpf(y + 1e-5 if logit else y + 1.0)
                b = mp.mpf((trials + 2e-5) - (y + 1e-5) if logit else trials - y + 1.0)
                argument = 1 / (1 + mp.exp(-x)) if logit else abs(x)
                normalization = mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b)
                density = (a - 1) * mp.log(argument) + (b - 1) * mp.log1p(-argument) - normalization
                value = density + (mp.log(mp.mpf(trials + 2e-5) / (a * b)) if logit else -mp.log(mp.mpf(trials + 1.0)))
            else:
                a = mp.mpf(number(model['df1'][0])) / 2
                b = mp.mpf(number(model['df2'][0])) / 2
                if mp.isinf(a):
                    value = b * mp.log(b) - mp.loggamma(b) - b * x - b * mp.exp(-x)
                elif mp.isinf(b):
                    value = a * mp.log(a) - mp.loggamma(a) + a * x - a * mp.exp(x)
                else:
                    z = mp.log(a / b) + x
                    value = a * z - (a + b) * mp.log1p(mp.exp(z)) - (mp.loggamma(a) + mp.loggamma(b) - mp.loggamma(a + b))
            values.append(mp.nstr(value, 90))
        rows.append(dict(name=case['name'], log_density=values))
    result = dict(source=dict(mpmath=mp.__version__, decimal_precision=mp.mp.dps,
        input_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
        method='independent density equations with high-precision log-gamma; f64-derived shape parameters'), cases=rows)
    args.fresh_output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
