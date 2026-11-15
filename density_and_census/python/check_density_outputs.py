#!/usr/bin/env python3
"""Independent checks of the Revision 4 computations (Table 2, Section 8).

Checks, by algorithms different from the producers:
  1. census: for every prime in the CSV, r_p is recomputed from the listed
     irregular indices; the r_p histogram and totals are recomputed; the
     irregular indices of a deterministic sample of primes (every 25th prime,
     plus 103, 157, 12613, 17681, 78233, 172243) are recomputed by inverting
     (e^x - 1)/x modulo p, and for p < 1500 from exact Bernoulli numbers.
  2. model: P(r_p = 1) = exp(-1/3) and P(r_p = 2) = exp(-1/5) - exp(-11/30)
     (closed forms of the inclusion-exclusion sums for k <= 2).
  3. Table 2 rows recomputed from census and model outputs.
  4. mechanism and construction cases: each stored vector v is re-read, its
     projections and character values are recomputed, every reconstruction
     congruence is checked, and Theta(v) is recomputed as
     s * prod_d Res(Phi_p, F_d) with SymPy resultants (Fourier factorization),
     independently of the FLINT determinant used by the producers.

Usage: check_density_outputs.py OUTPUT_DIR REPORT_JSON
"""
import csv
import json
import math
import sys

from flint import fmpq, nmod_poly
from sympy import Poly, cyclotomic_poly, resultant, symbols

X = symbols("X")


def primes_below(n):
    sieve = bytearray([1]) * n
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n) if sieve[i]]


def irregular_by_exp_inverse(p):
    n = p - 2
    fact, coeffs = 1, []
    for k in range(n):
        fact = fact * (k + 1) % p
        coeffs.append(pow(fact, p - 2, p))
    g = nmod_poly(coeffs, p).inverse_series_trunc(n).coeffs() + [0] * n
    return [k // 2 for k in range(2, p - 2, 2) if int(g[k]) == 0]


def r_min(p, half):
    r = 1
    while any(i % r == 0 and (2 * r + 1) * i <= r * (p - 1) for i in half):
        r += 1
    return r


def theta_by_factorization(vec, p, s):
    v = {(i, j): vec[i * p + j] for i in range(p) for j in range(p)}
    phi = Poly(cyclotomic_poly(p, X), X)
    total = s
    pis = {}
    for d in list(range(p)) + ["inf"]:
        pi = [0] * p
        for (i, j), x in v.items():
            pi[j % p if d == "inf" else (i + d * j) % p] += x
        pis[d] = pi
        total *= int(resultant(phi, Poly(sum(c * X ** t for t, c in enumerate(pi)), X)))
    congruences = all((pis["inf"][j] + sum(pis[b][(i + b * j) % p] for b in range(p)) - s) % p == 0
                      for i in range(p) for j in range(p))
    return total, congruences, sum(vec) == s


def main():
    out, report = sys.argv[1], sys.argv[2]
    res = {}
    # 1. census
    summ = json.load(open(f"{out}/census_200000.json"))
    rows = list(csv.reader(open(f"{out}/census_200000.csv"), delimiter=";"))[1:]
    hist, bad_r, irr = {}, [], 0
    half_of = {}
    for p_s, idx, r_s, _ in rows:
        p = int(p_s)
        half = [int(x) // 2 for x in idx.split()]
        half_of[p] = half
        irr += bool(half)
        r = r_min(p, half)
        if r != int(r_s):
            bad_r.append(p)
        hist[r] = hist.get(r, 0) + 1
    primes = [q for q in primes_below(200000) if q >= 7]
    sample = sorted(set(primes[::25]) | {103, 157, 12613, 17681, 78233, 172243})
    mism = [p for p in sample if irregular_by_exp_inverse(p) != half_of[p]]
    exact = [p for p in primes if p < 1500 and
             [k // 2 for k in range(2, p - 2, 2) if int(fmpq.bernoulli(k).p % p) == 0] != half_of[p]]
    res["census"] = {
        "rows_cover_all_primes_7_to_200000": [int(r[0]) for r in rows] == primes,
        "r_p_recomputed_mismatches": bad_r,
        "histogram_matches_summary": {str(k): hist[k] for k in sorted(hist)} == summ["r_p_counts"],
        "irregular_count_matches": irr == summ["irregular_primes"],
        "sample_size_exp_inverse": len(sample), "sample_mismatches": mism,
        "exact_bernoulli_p_lt_1500_mismatches": exact,
        "p172243": rows[[int(r[0]) for r in rows].index(172243)],
    }
    # 2. model
    model = json.load(open(f"{out}/model_k12.json"))
    p1, p2 = model["P(r_p = k)"]["1"], model["P(r_p = k)"]["2"]
    res["model"] = {
        "P1_minus_exp(-1/3)": abs(p1 - math.exp(-1 / 3)),
        "P2_minus_closed_form": abs(p2 - (math.exp(-1 / 5) - math.exp(-11 / 30))),
    }
    # 3. Table 2
    n = summ["primes"]
    table = [[str(k), summ["r_p_counts"].get(str(k), 0), round(n * model["P(r_p = k)"][str(k)], 1)]
             for k in range(1, 10)]
    table.append([">=10", sum(v for k, v in summ["r_p_counts"].items() if int(k) >= 10),
                  round(n * model["P(r_p > k)"]["9"], 1)])
    res["table_rows"] = table
    # 4. vectors
    vec_checks = []
    for name in ("mechanism_test.json", "theorem_d_construction.json"):
        data = json.load(open(f"{out}/{name}"))
        for c in data["cases"]:
            if "vector_row_major" not in c:
                vec_checks.append({"file": name, "ok": False, "reason": "no vector"})
                continue
            p = c["p"]
            s = c.get("augmentation", c.get("s_prime"))
            th, cong, aug = theta_by_factorization(c["vector_row_major"], p, s)
            ok = cong and aug and th == int(c["theta"])
            if "m" in c:
                ok = ok and th == p ** c["N"] * int(c["m"])
            vec_checks.append({"file": name, "p": p, "ok": ok})
    res["vectors"] = {"checked": len(vec_checks), "all_ok": all(v["ok"] for v in vec_checks)}
    res["all_passed"] = (res["census"]["rows_cover_all_primes_7_to_200000"]
                         and not bad_r and not mism and not exact
                         and res["census"]["histogram_matches_summary"]
                         and res["census"]["irregular_count_matches"]
                         and res["model"]["P1_minus_exp(-1/3)"] < 1e-12
                         and res["model"]["P2_minus_closed_form"] < 1e-12
                         and res["vectors"]["all_ok"])
    json.dump(res, open(report, "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("all_passed", "vectors", "model")}, indent=1))
    print("table rows:", table)


if __name__ == "__main__":
    main()
