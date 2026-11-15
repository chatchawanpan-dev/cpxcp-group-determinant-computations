#!/usr/bin/env python3
"""Explicit instances of the construction in the proof of Theorem D.

Start from v0 with v_p(Theta(v0)) = p + 3 in the coordinates of Lemma 2.1
(exceptional direction infinity), put y_b = f_b/varpi and
y_inf = g_t/varpi^2 with g_t = f_inf + (tau0 - t) varpi^(p-1).  For each
direction search y'_d = y_d + p w (w in R small) until N(y'_d) is a rational
prime l_d, distinct from the others; y'_d R is then a prime ideal in the ray
class of y_d R modulo pR.  With m = m0 prod l_d and a = N - p - 2, the
augmentation s = p^a m0 and the character values varpi y'_b, varpi^2 y'_inf
give, by (2.12) and (2.5), an integral vector v with Theta(v) = p^N m
(exact determinant).

Usage: theorem_d_construction.py SEED OUT_JSON
"""
import json
import random
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
from flint import fmpz  # noqa: E402
from cpxcp_lib import (add, augmentation, canon, char_value, directions,  # noqa: E402
                       div_by_varpi, group_det, lift, mul, norm, projection,
                       reconstruct, scal, varpi_pow, vp)
from mechanism_test import random_v0  # noqa: E402


def construct(p, N, m0, rng):
    a = N - p - 2
    assert a >= 1 and m0 % p
    v0, th0 = random_v0(p, rng)
    s0 = augmentation(v0)
    tau0 = (s0 // p) % p
    t = m0 % p if a == 1 else 0
    f = {d: char_value(v0, d, p) for d in directions(p)}
    y = {b: div_by_varpi(f[b], p) for b in range(p)}
    g = add(f["inf"], scal((tau0 - t) % p, varpi_pow(p - 1, p), p), p)
    y["inf"] = div_by_varpi(div_by_varpi(g, p), p)
    primes, ynew = {}, {}
    used = set()
    for d in directions(p):
        while True:
            w = [rng.randint(-1, 1) for _ in range(p)]
            cand = add(y[d], scal(p, w, p), p)
            n = norm(cand, p)
            if n > 1 and n not in used and n % p and fmpz(n).is_prime():
                break
        used.add(n)
        primes[str(d)] = n
        ynew[d] = cand
    m = m0
    for n in primes.values():
        m *= n
    s = p ** a * m0
    fv = {b: mul(varpi_pow(1, p), ynew[b], p) for b in range(p)}
    fv["inf"] = mul(varpi_pow(2, p), ynew["inf"], p)
    pis = {d: lift(fv[d], s, p) for d in directions(p)}
    v, ok = reconstruct(pis, s, p)
    rec = {"p": p, "N": N, "m0": m0, "t": t, "tau0": tau0,
           "primes_l_delta": primes, "all_l_delta_1_mod_p": all(n % p == 1 for n in primes.values()),
           "integral": ok}
    if ok:
        th = group_det(v, p)
        rec.update({
            "augmentation_matches": augmentation(v) == s,
            "character_values_match": all(canon(projection(v, d, p), p) == canon(fv[d], p)
                                          for d in directions(p)),
            "determinant_equals_p^N_m": th == p ** N * m,
            "theta_valuation": vp(th, p),
            "max_abs_entry": max(abs(x) for x in v.values()),
            "m": str(m), "augmentation": s, "theta": str(th),
            "vector_row_major": [v[(i, j)] for i in range(p) for j in range(p)],
        })
    return rec


def main():
    seed, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(seed)
    cases = [(5, 8, 1), (5, 9, 2), (7, 10, 1), (7, 10, 3), (7, 10, -2), (7, 11, 1),
             (7, 11, 5), (11, 14, 1), (11, 14, 7), (11, 15, -1), (11, 16, 2)]
    rows = [construct(p, N, m0, rng) for p, N, m0 in cases]
    passed = all(r["integral"] and r["all_l_delta_1_mod_p"] and r["augmentation_matches"]
                 and r["character_values_match"] and r["determinant_equals_p^N_m"] for r in rows)
    summary = {"cases": len(rows), "seed": seed, "all_passed": passed}
    json.dump({"summary": summary, "cases": rows}, open(out, "w"), indent=1)
    print(json.dumps(summary))
    for r in rows:
        print(r["p"], r["N"], r["m0"], r.get("determinant_equals_p^N_m"), r.get("max_abs_entry"))


if __name__ == "__main__":
    main()
