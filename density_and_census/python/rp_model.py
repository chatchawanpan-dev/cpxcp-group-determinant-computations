#!/usr/bin/env python3
"""Limiting distribution of r_p in the independence model for irregular pairs.

Model: each half-index i in {1,...,(p-3)/2} is irregular independently with
probability 1/p.  As p -> infinity the points i/p of irregular half-indices form
a Poisson process of intensity 1 on (0, 1/2), with residues i mod L uniformly
distributed for every L.  An integer r is inadmissible exactly when some
irregular i satisfies r | i and i/(p-1) <= theta_r = r/(2r+1).  Hence

  P(r_p > k) = sum_{S subset {1..k}} (-1)^{|S|} exp(-mu(S)),
  mu(S) = int_0^{1/2} dens{n : r | n for some r in S with theta_r >= x} dx,

computed exactly in rationals; exp is evaluated with mpmath at 60 digits.

Usage: rp_model.py KMAX OUT_JSON
"""
import json
import sys
from fractions import Fraction
from math import lcm

import mpmath
import numpy as np

mpmath.mp.dps = 60


def main():
    kmax = int(sys.argv[1])
    out = sys.argv[2]
    L = 1
    for r in range(1, kmax + 1):
        L = lcm(L, r)
    n = np.arange(L, dtype=np.int64)
    masks = np.zeros(L, dtype=np.int64)
    for r in range(1, kmax + 1):
        masks |= np.where(n % r == 0, 1 << (r - 1), 0)
    dens = {}
    for A in range(1 << kmax):
        dens[A] = Fraction(int(np.count_nonzero(masks & A)), L)
    theta = {r: Fraction(r, 2 * r + 1) for r in range(1, kmax + 1)}

    def mu(S):
        members = [r for r in range(1, kmax + 1) if S >> (r - 1) & 1]
        total, prev = Fraction(0), Fraction(0)
        for j, r in enumerate(members):
            suffix = 0
            for q in members[j:]:
                suffix |= 1 << (q - 1)
            total += (theta[r] - prev) * dens[suffix]
            prev = theta[r]
        return total

    tail = []      # tail[k] = P(r_p > k)
    for k in range(0, kmax + 1):
        acc = mpmath.mpf(0)
        for S in range(1 << k):
            sign = -1 if bin(S).count("1") % 2 else 1
            acc += sign * mpmath.e ** (-mpmath.mpf(mu(S).numerator) / mu(S).denominator)
        tail.append(acc)
    probs = {str(k): float(tail[k - 1] - tail[k]) for k in range(1, kmax + 1)}
    res = {
        "model": "Poisson limit of independent irregular half-indices",
        "kmax": kmax,
        "P(r_p = k)": probs,
        "P(r_p > k)": {str(k): float(tail[k]) for k in range(1, kmax + 1)},
        "check_P(r_p=1)=exp(-1/3)": float(mpmath.e ** (-mpmath.mpf(1) / 3)),
    }
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
