#!/usr/bin/env python3
"""Independent 90-digit Student-t posterior references for finite uniform priors."""

import argparse
from decimal import Decimal, localcontext
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from generate_beta_mathematical_reference import gamma, number


def decimal_value(value):
    if value.startswith("bits:"):
        return Decimal.from_float(number(value))
    return Decimal(value)


@lru_cache
def beta(a, b):
    return gamma(a) * gamma(b) / gamma(a + b)


def incomplete_beta(a, b, x):
    if x == 1:
        return beta(a, b)
    if x > Decimal("0.5"):
        return beta(a, b) - incomplete_beta(b, a, 1 - x)
    # Integrate the binomial series for (1-t)^(b-1) term by term.
    coefficient = Decimal(1)
    total = 1 / a
    for n in range(1, 2049):
        coefficient *= (Decimal(n) - b) * x / n
        term = coefficient / (a + n)
        total += term
        if abs(term) <= abs(total) * Decimal("1e-92"):
            return x ** a * total
    raise ArithmeticError("Student-t beta series did not converge")


@lru_cache
def tail(q, df, order):
    a = (df - order) / 2
    b = Decimal(order + 1) / 2
    x = df / (df + q * q)
    return df.sqrt() ** order * incomplete_beta(a, b, x)


def interval(a, b, df, order):
    if a >= 0:
        return tail(a, df, order) - tail(b, df, order)
    if b <= 0:
        return (-1) ** order * (tail(-b, df, order) - tail(-a, df, order))
    center = tail(Decimal(0), df, order)
    return center - tail(b, df, order) + (-1) ** order * (center - tail(-a, df, order))


def observation(components, estimate, se, df, thresholds):
    normalization = tail(Decimal(0), df, 0)
    root_df = df.sqrt()
    log_likelihoods, weights = [], []
    first = second = Decimal(0)
    cumulative = [Decimal(0) for _ in thresholds]
    for prior, lower, upper in components:
        a = (lower - estimate) / se
        if lower == upper:
            shape = 1 + a * a / df
            likelihood = 1 / (root_df * normalization * se * shape ** ((df + 1) / 2))
            mean, raw_second = lower, lower * lower
            probabilities = [Decimal(t >= upper) for t in thresholds]
        else:
            b = (upper - estimate) / se
            mass = interval(a, b, df, 0)
            likelihood = mass / (2 * normalization * (upper - lower))
            # Preserve the undefined moments in the public R-compatible
            # contract when the full Student-t moment does not exist.
            mean = raw_second = Decimal("NaN")
            if df > 1:
                mean_delta = se * interval(a, b, df, 1) / mass
                mean = estimate + mean_delta
                if df > 2:
                    raw_second = estimate * (estimate + 2 * mean_delta) + se * se * interval(a, b, df, 2) / mass
            probabilities = []
            for threshold in thresholds:
                if threshold <= lower:
                    probabilities.append(Decimal(0))
                elif threshold >= upper:
                    probabilities.append(Decimal(1))
                else:
                    t = (threshold - estimate) / se
                    probabilities.append(interval(a, t, df, 0) / mass)
        log_likelihoods.append(likelihood.ln())
        weight = prior * likelihood
        weights.append(weight)
        first += weight * mean
        second += weight * raw_second
        for i, probability in enumerate(probabilities):
            cumulative[i] += weight * probability
    total = sum(weights)
    mean = first / total
    variance = second / total - mean * mean
    sd = variance if variance.is_nan() else max(variance, Decimal(0)).sqrt()
    return {
        "log_likelihoods": encoded(log_likelihoods),
        "posterior_weights": encoded(w / total for w in weights),
        "posterior_mean": encode(mean),
        "posterior_sd": encode(sd),
        "cdf": encoded(p / total for p in cumulative),
    }


def encode(value):
    return "NaN" if value.is_nan() else format(float(value), ".17g")


def encoded(values):
    return [encode(value) for value in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    source_bytes = args.input.read_bytes()
    source = json.loads(source_bytes)
    cases = {}
    with localcontext() as context:
        context.prec = 100
        for case in source["cases"]:
            df = decimal_value(case["df"])
            components = [
                tuple(decimal_value(value) for value in values)
                for values in zip(case["weights"], case["lower"], case["upper"], strict=True)
            ]
            thresholds = [decimal_value(t) for t in case["thresholds"]]
            rows = []
            for estimate, se in zip(case["estimates"], case["standard_errors"], strict=True):
                row = observation(components, decimal_value(estimate), decimal_value(se), df, thresholds)
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
