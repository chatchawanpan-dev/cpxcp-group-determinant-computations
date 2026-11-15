#!/usr/bin/env python3
"""Census of the least full-tail exponent e(p) = p + 3 + 2 r_p for primes 7 <= p <= PMAX.

Irregular half-indices are computed from B_k/k! = [x^k] x/(e^x-1) modulo p (series inversion
in F_p[x]); r_p is the least r >= 1 such that (2r+1) i > r(p-1) for every irregular half-index i
divisible by r.  Optionally cross-checks against exact Bernoulli numbers for p <= CHECK.
"""
import sys, json, math, csv
from flint import nmod_poly, fmpq
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
CHECK = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = sys.argv[3] if len(sys.argv) > 3 else None
ALL_ROWS_OUT = sys.argv[4] if len(sys.argv) > 4 else None
if PMAX < 103 or CHECK < 0:
    raise ValueError('PMAX must be at least 103 for the first-exception summary; CHECK must be nonnegative')

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def irregular_half_indices(p):
    n = p - 2                                   # need coefficients up to x^(p-3)
    fact = 1; coeffs = []
    for k in range(n):                          # (e^x - 1)/x = sum x^k/(k+1)!
        fact = fact * (k + 1) % p
        coeffs.append(pow(fact, p - 2, p))
    f = nmod_poly(coeffs, p)
    g = f.inverse_series_trunc(n)               # x/(e^x-1) = sum B_k x^k / k!
    gc = g.coeffs() + [0] * n
    return [k // 2 for k in range(2, p - 2, 2) if int(gc[k]) == 0]

def r_min(p, I):
    r = 1
    while not all((2 * r + 1) * i > r * (p - 1) for i in I if i % r == 0):
        r += 1
    return r

if CHECK:
    for p in primes_upto(CHECK):
        if p < 7: continue
        exact = [k // 2 for k in range(2, p - 2, 2) if int(fmpq.bernoulli(k).p % p) == 0]
        assert exact == irregular_half_indices(p), p
rows = []
all_rows = []
hist = {}
for p in primes_upto(PMAX):
    if p < 7: continue
    I = irregular_half_indices(p)
    r = r_min(p, I)
    hist[r] = hist.get(r, 0) + 1
    if ALL_ROWS_OUT:
        all_rows.append([p, json.dumps([2*i for i in I]), r, p+3+2*r, int(p < 2**31)])
    if I:
        rows.append({"p": p, "irregular_indices": [2 * i for i in I], "r_p": r, "tail_exponent": p + 3 + 2 * r})
tot = sum(hist.values())
summary = {"PMAX": PMAX, "primes": tot, "irregular": len(rows),
           "r_p_histogram": {str(k): hist[k] for k in sorted(hist)},
           "fraction_r_p_gt_1": (tot - hist.get(1, 0)) / tot,
           "heuristic_1_minus_exp(-1/3)": 1 - math.exp(-1 / 3),
           "max_r_p": max(hist), "max_r_p_le_(p-5)/2_for_all": all(rw["r_p"] <= (rw["p"] - 5) // 2 for rw in rows),
           "smallest_with_r_p_gt_1": next(rw["p"] for rw in rows if rw["r_p"] > 1),
           "max_index_of_irregularity": max(len(rw["irregular_indices"]) for rw in rows)}
print(json.dumps(summary, indent=1))
if OUT:
    json.dump({"summary": summary, "irregular_primes": rows}, open(OUT, "w"), indent=0)
if ALL_ROWS_OUT:
    with open(ALL_ROWS_OUT, 'w', newline='') as stream:
        writer = csv.writer(stream, delimiter=';')
        writer.writerow(['p','irregular_indices','r_p','p_plus_3_plus_2r_p','within_cited_vandiver_range'])
        writer.writerows(all_rows)
