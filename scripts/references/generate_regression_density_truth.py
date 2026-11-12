#!/usr/bin/env python3
"""Independent high precision density checks for R tail/cancellation differences.

Requires mpmath (development reference generation only). Never imports candidate
code or uses its outputs. Preserves the separate original-R observations.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import struct
import mpmath as mp


def decode(value):
    if value.startswith("bits:"):
        return mp.mpf(struct.unpack(">d", bytes.fromhex(value[5:]))[0])
    return mp.mpf(value)


def student(x, df, ncp):
    # Integrate the defining independent Normal/chi-square mixture. Centering
    # the transformed chi variable keeps the df=1e10 case numerically resolved.
    a = 1 + x*x/df
    k = ncp*x/mp.sqrt(df+x*x)
    mode = (k + mp.sqrt(k*k+4*df))/2
    peak = df*mp.log(mode)-mode*mode/2+k*mode
    def f(u):
        z = mode+u
        return mp.exp(df*mp.log(z)-z*z/2+k*z-peak) if z > 0 else mp.mpf(0)
    area = mp.quad(f, [-mode, -min(mode, 20), 0, 20, mp.inf])
    return (mp.log(2)-mp.log(2*mp.pi)/2-df*mp.log(2)/2
            -mp.loggamma(df/2)-mp.log(df)/2-(df+1)*mp.log(a)/2
            -ncp*ncp/2+peak+mp.log(area))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Use a fresh output path")
    mp.mp.dps = 90
    reference = json.loads(args.reference.read_text())
    leaves = []
    for c in reference["leaves"]:
        name = c["name"]
        if name.startswith("nb_"):
            size, mu = {"nb_large_size": (10**8, 4), "nb_large_count": (30, 10**6),
                        "nb_probability": (3, mp.mpf(3)*(1-mp.mpf(float(.4)))/mp.mpf(float(.4)))}[name]
            size, mu = mp.mpf(size), mp.mpf(mu)
            def density(x):
                return (mp.loggamma(x+size)-mp.loggamma(size)-mp.loggamma(x+1)
                        -size*mp.log1p(mu/size)+x*(mp.log(mu)-mp.log(size+mu)))
        elif name.startswith("student_") and name.endswith("_0.7") and "Inf" not in name:
            df = mp.mpf(name.split("_")[1])
            ncp = mp.mpf(float(.7))
            def density(x):
                return student(x, df, ncp)
        else:
            continue
        leaves.append({"name": name, "expected": [mp.nstr(density(decode(v)), 80) for v in c["observations"]]})
    args.output.write_text(json.dumps({"source": {"method": "defining Normal/chi-square integral and log-gamma PMF",
        "mpmath": mp.__version__, "decimal_precision": mp.mp.dps}, "leaves": leaves}, indent=2)+"\n")


if __name__ == "__main__":
    main()
