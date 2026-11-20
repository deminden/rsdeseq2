#!/usr/bin/env python3
"""Independent Gaussian posterior references for finite uniform priors."""

import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from generate_student_uniform_mathematical_reference import decimal_value, encode, encoded


def inverse_atan(n):
    inverse_square = 1 / Decimal(n * n)
    power = total = 1 / Decimal(n)
    for k in range(1, 1025):
        power *= -inverse_square
        term = power / (2 * k + 1)
        total += term
        if abs(term) < Decimal("1e-358"):
            return total
    raise ArithmeticError("Arctangent series did not converge")


@lru_cache
def normalization():
    with localcontext() as ctx:
        ctx.prec = 360
        return (32 * inverse_atan(5) - 8 * inverse_atan(239)).sqrt()


def density(q):
    return (-q * q / 2).exp()


@lru_cache
def tail(q):
    if q < 24:
        # Integrate the exponential series about zero. Use extra precision
        # for cancellation both inside the series and when forming the tail.
        with localcontext() as ctx:
            ctx.prec = 360
            power = total = q
            square_half = q * q / 2
            for k in range(1, 4097):
                power *= -square_half / k
                term = power / (2 * k + 1)
                total += term
                if k > square_half and abs(term) < Decimal("1e-345"):
                    result = normalization() / 2 - total
                    break
            else:
                raise ArithmeticError("Gaussian integral series did not converge")
        return +result
    # Repeated integration by parts bounds the scaled tail by successive
    # partial sums. At 24 and above the terms reach the tolerance before
    # their magnitudes start increasing.
    square = q * q
    term = total = Decimal(1)
    for k in range(1, 4097):
        term *= -Decimal(2 * k - 1) / square
        total += term
        if abs(term) < Decimal("1e-100"):
            return density(q) * total / q
    raise ArithmeticError("Gaussian tail expansion did not converge")


def interval(a, b):
    if a >= 0:
        return tail(a) - tail(b)
    if b <= 0:
        return tail(-b) - tail(-a)
    return normalization() - tail(-a) - tail(b)


def observation(components, estimate, se, thresholds):
    log_likelihoods, weights = [], []
    first = second = Decimal(0)
    cumulative = [Decimal(0) for _ in thresholds]
    for prior, lower, upper in components:
        a = (lower - estimate) / se
        if lower == upper:
            likelihood = density(a) / (normalization() * se)
            mean, raw_second = lower, lower * lower
            probabilities = [Decimal(t >= upper) for t in thresholds]
        else:
            b = (upper - estimate) / se
            mass = interval(a, b)
            likelihood = mass / (normalization() * (upper - lower))
            left, right = density(a), density(b)
            mean_delta = se * (left - right) / mass
            mean = estimate + mean_delta
            raw_second = estimate * (estimate + 2 * mean_delta) + se * se * (1 + (a * left - b * right) / mass)
            probabilities = []
            for threshold in thresholds:
                if threshold <= lower:
                    probabilities.append(Decimal(0))
                elif threshold >= upper:
                    probabilities.append(Decimal(1))
                else:
                    t = (threshold - estimate) / se
                    probabilities.append(interval(a, t) / mass)
        log_likelihoods.append(likelihood.ln())
        weight = prior * likelihood
        weights.append(weight)
        first += weight * mean
        second += weight * raw_second
        for i, probability in enumerate(probabilities):
            cumulative[i] += weight * probability
    total = sum(weights)
    mean = first / total
    variance = max(second / total - mean * mean, Decimal(0))
    return {
        "log_likelihoods": encoded(log_likelihoods),
        "posterior_weights": encoded(w / total for w in weights),
        "posterior_mean": encode(mean),
        "posterior_sd": encode(variance.sqrt()),
        "cdf": encoded(p / total for p in cumulative),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    cases = {}
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.Emin = -999999999
        for case in source["cases"]:
            components = [
                tuple(decimal_value(value) for value in values)
                for values in zip(case["weights"], case["lower"], case["upper"], strict=True)
            ]
            thresholds = [decimal_value(t) for t in case["thresholds"]]
            rows = []
            for estimate, se in zip(case["estimates"], case["standard_errors"], strict=True):
                row = observation(components, decimal_value(estimate), decimal_value(se), thresholds)
                rows.append({"estimate": estimate, "standard_error": se, **row})
            cases[case["name"]] = rows
    fixture = {
        "source": {
            "inputs": args.input.name,
            "input_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "working_decimal_digits": 100,
        },
        "cases": cases,
    }
    with args.output.open("x") as output:
        output.write(json.dumps(fixture, indent=2) + "\n")


if __name__ == "__main__":
    main()
