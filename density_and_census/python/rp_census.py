#!/usr/bin/env python3
"""Census of irregular pairs and of r_p for primes 7 <= p < PMAX.

Method (independent of python/tail_exponent_census.py, which inverts
(e^x-1)/x): with y = x^2/4,
    (x/2) coth(x/2) = cosh(x/2) / (sinh(x/2)/(x/2)) = sum_k B_{2k} 4^k/(2k)! y^k,
    sinh(x/2)/(x/2) = sum_k y^k/(2k+1)!,   cosh(x/2) = sum_k y^k/(2k)!.
For 1 <= k <= (p-3)/2 the coefficient of y^k is a unit multiple of B_{2k}
modulo p, so (p, 2k) is an irregular pair exactly when it vanishes mod p.
Inverse factorials come from Wilson's theorem: 1/n! = (-1)^(n+1) (p-1-n)!.

r_p is the least r >= 1 with (2r+1) i > r (p-1) for every irregular
half-index i divisible by r (manuscript, Section 1).

Usage: rp_census.py PMAX OUT_JSON OUT_CSV [JOBS]
"""
import csv
import json
import sys
from multiprocessing import Pool

import numpy as np
from flint import nmod_poly


def primes_below(n):
    s = np.ones(n, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(q) for q in np.nonzero(s)[0]]


def factorials_mod(p):
    """Array F with F[n] = n! mod p for 0 <= n <= p-1 (vectorized prefix product)."""
    a = np.arange(p, dtype=np.int64)
    a[0] = 1
    d = 1
    while d < p:
        a[d:] = (a[d:] * a[:-d]) % p
        d *= 2
    return a


def irregular_half_indices(p):
    F = factorials_mod(p)
    n = np.arange(p - 1)
    inv = F[p - 1 - n] % p                       # (p-1-n)!
    inv = np.where(n % 2 == 1, inv, (p - inv) % p)   # sign (-1)^(n+1)
    L = (p - 1) // 2                             # need y^k for k <= (p-3)/2
    sinh_c = inv[1::2][:L]                       # 1/(2k+1)!
    cosh_c = inv[0::2][:L]                       # 1/(2k)!
    s = nmod_poly(sinh_c.tolist(), p)
    c = nmod_poly(cosh_c.tolist(), p)
    g = (c * s.inverse_series_trunc(L)).coeffs()
    g = g + [0] * (L - len(g))
    return [k for k in range(1, (p - 3) // 2 + 1) if int(g[k]) == 0]


def r_min(p, I):
    r = 1
    while any(i % r == 0 and (2 * r + 1) * i <= r * (p - 1) for i in I):
        r += 1
    return r


def work(p):
    I = irregular_half_indices(p)
    return p, I, r_min(p, I)


def main():
    pmax = int(sys.argv[1])
    out_json, out_csv = sys.argv[2], sys.argv[3]
    jobs = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    ps = [q for q in primes_below(pmax) if q >= 7]
    with Pool(jobs) as pool:
        rows = pool.map(work, ps, chunksize=16)
    rows.sort()
    hist = {}
    for p, I, r in rows:
        hist[r] = hist.get(r, 0) + 1
    irregular = [row for row in rows if row[1]]
    rmax = max(r for _, _, r in rows)
    summary = {
        "range": f"7 <= p < {pmax}",
        "primes": len(rows),
        "irregular_primes": len(irregular),
        "irregular_pairs": sum(len(I) for _, I, _ in rows),
        "r_p_counts": {str(k): hist[k] for k in sorted(hist)},
        "max_r_p": rmax,
        "primes_with_max_r_p": [[p, [2 * i for i in I]] for p, I, r in rows if r == rmax],
        "max_index_of_irregularity": max(len(I) for _, I, _ in rows),
        "r_p_le_(p-5)/2_for_all": all(r <= (p - 5) // 2 for p, _, r in rows),
    }
    # counts below 20000, for comparison with Section 8 of Revision 3
    h2 = {}
    for p, I, r in rows:
        if p < 20000:
            h2[r] = h2.get(r, 0) + 1
    summary["below_20000"] = {
        "primes": sum(h2.values()),
        "irregular_primes": sum(1 for p, I, r in rows if p < 20000 and I),
        "r_p_counts": {str(k): h2[k] for k in sorted(h2)},
    }
    json.dump(summary, open(out_json, "w"), indent=1)
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow(["p", "irregular_indices_2i", "r_p", "e_p_minus_p_if_vandiver"])
        for p, I, r in rows:
            w.writerow([p, " ".join(str(2 * i) for i in I), r, 3 + 2 * r])
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
