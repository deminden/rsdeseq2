#!/usr/bin/env python3
"""Independent 90-digit normalizers and beta CDFs for the Student-t leaf panel.

Only the input bits are read from the R panel. Gamma uses the integrated
exponential series and positive recurrence. Beta uses an exact rational
normalizer and the integrated binomial series.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import struct


def number(value):
    return struct.unpack('>d', bytes.fromhex(value.removeprefix('bits:')))[0]


def rational_beta(a, b):
    a, b = Fraction(a), Fraction(b)
    if a == 1:
        return 1 / b
    if b == 1:
        return 1 / a
    if a in (Fraction(1, 2), Fraction(3, 2)):
        return rational_beta(b, a)
    if a.denominator == 1 and b in (Fraction(1, 2), Fraction(3, 2)):
        value = Fraction(2)
        for k in range(1, int(a)):
            value *= Fraction(2 * k, 2 * k + 1)
        return value if b == Fraction(1, 2) else value / (2 * a + 1)
    raise ValueError('This oracle supports only its rational beta shape family')


def gamma(x):
    # Reduce to 0 < x <= 1. The integral from 256 to infinity is at most
    # exp(-256) there. Leave 111 digits for cancellation in the series.
    with localcontext() as context:
        context.prec = 240
        scale = Decimal(1)
        while x > 1:
            x -= 1
            scale *= x
        cutoff = Decimal(256)
        coefficient = Decimal(1)
        total = 1 / x
        for n in range(1, 2049):
            coefficient *= -cutoff / n
            term = coefficient / (x + n)
            total += term
            if n > cutoff and abs(term) < Decimal('1e-115'):
                break
        else:
            raise ArithmeticError('Independent Gamma integral did not converge')
        result = scale * cutoff ** x * total
    return +result


def cdf(a, b, x):
    normalizer = rational_beta(a, b)
    normalizer = Decimal(normalizer.numerator) / Decimal(normalizer.denominator)
    a, b = Decimal(a), Decimal(b)
    reflect = x >= (a + 1) / (a + b + 2)
    if reflect:
        a, b, x = b, a, 1 - x
    coefficient = Decimal(1)
    total = 1 / a
    for n in range(1, 10001):
        coefficient *= (Decimal(n) - b) * x / n
        term = coefficient / (a + n)
        total += term
        if not term or abs(term) < abs(total) * Decimal('1e-85'):
            break
    else:
        raise ArithmeticError('Independent binomial series did not converge')
    exponent = int(a)
    power = x ** exponent
    if a != exponent:
        if a - exponent != Decimal('.5'):
            raise ValueError('Only integer/half-integer powers are supported')
        power *= x.sqrt()
    value = power * total / normalizer
    return 1 - value if reflect else value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a fresh output path; references are never overwritten')
    source = json.loads(args.input.read_text())
    cases = []
    with localcontext() as context:
        context.prec = 90
        arguments = source['gamma_arguments']
        gamma_values = [format(float(gamma(Decimal.from_float(number(x)))), '.17g')
                        for x in arguments]
        for case in source['cases']:
            a, b, x = (Decimal.from_float(number(case[k])) for k in ('a', 'b', 'p'))
            value = cdf(a, b, x)
            cases.append({k: case[k] for k in ('a', 'b', 'p')} | {
                'beta': format(float(rational_beta(a, b)), '.17g'),
                'regularized': format(float(value), '.17g'),
                'decimal': str(value),
            })
    args.output.write_text(json.dumps({'source': {
        'precision_decimal_digits': 90,
        'algorithm': 'integrated exponential series, exact rational beta and integrated binomial series',
        'inputs': args.input.name,
    }, 'gamma_arguments': arguments, 'gamma_values': gamma_values, 'cases': cases}, indent=2) + '\n')


if __name__ == '__main__':
    main()
