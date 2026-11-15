#!/usr/bin/env python3
"""Finite check of the lifting mechanism behind Theorem D (Steps 1-3).

For p in {5, 7, 11}: draw v0 in Z^G with augmentation p, exceptional
direction infinity (sum v_ij j = 0 mod p, sum v_ij i != 0 mod p) and
v_p(Theta(v0)) = p + 3.  Replace every character value f_d by
f_d (1 + p beta_d) with random beta_d in R, add c varpi^(p-1) to f_inf, and
take an augmentation s' with s'/p = tau0 - c (mod p).  The vectors given by
(2.12) and (2.5) must be integral, have the new character values, and
Theta(v') = s' prod N(f'_d) (exact p^2 x p^2 determinant).  As a control,
the residue s'/p = tau0 - c + 1 must break integrality.

Usage: mechanism_test.py TRIALS_PER_PRIME SEED OUT_JSON
"""
import json
import random
import sys

sys.path.insert(0, __import__("os").path.dirname(__file__))
from cpxcp_lib import (add, augmentation, canon, char_value, directions,  # noqa: E402
                       group_det, lift, mul, norm, one, projection,
                       reconstruct, scal, varpi_pow, vp, vp_norm)


def random_v0(p, rng):
    while True:
        v = {(i, j): rng.randint(-2, 2) for i in range(p) for j in range(p)}
        d = (-sum(x * j for (i, j), x in v.items())) % p
        v[(0, 1)] += d
        v[(0, 0)] -= d
        if sum(x * i for (i, j), x in v.items()) % p == 0:
            v[(1, 0)] += 1
            v[(0, 0)] -= 1
        v[(0, 0)] += p - augmentation(v)
        if vp_norm(char_value(v, "inf", p), p) != 2:
            continue
        th = group_det(v, p)
        if th != 0 and vp(th, p) == p + 3:
            return v, th


def trial(p, rng):
    v0, th0 = random_v0(p, rng)
    s0 = augmentation(v0)
    tau0 = (s0 // p) % p
    f = {d: char_value(v0, d, p) for d in directions(p)}
    c = rng.randrange(p)
    fp = {}
    for d in directions(p):
        beta = [rng.randint(-3, 3) for _ in range(p)]
        fp[d] = mul(f[d], add(one(p), scal(p, beta, p), p), p)
    fp["inf"] = add(fp["inf"], scal(c, varpi_pow(p - 1, p), p), p)
    k0 = rng.randint(-3, 3)
    s1 = p * ((tau0 - c) % p + p * k0)
    if s1 == 0:
        s1 = p * ((tau0 - c) % p + p)
    pis = {d: lift(fp[d], s1, p) for d in directions(p)}
    v1, ok = reconstruct(pis, s1, p)
    rec = {"p": p, "theta_v0_valuation": vp(th0, p), "c": c, "s_prime": s1,
           "integral": ok}
    if ok:
        same = all(canon(projection(v1, d, p), p) == canon(fp[d], p) for d in directions(p))
        th1 = group_det(v1, p)
        prod = s1
        for d in directions(p):
            prod *= norm(fp[d], p)
        rec.update({"character_values_match": same, "determinant_matches": th1 == prod,
                    "theta_valuation": vp(th1, p) if th1 else None,
                    "theta": str(th1),
                    "vector_row_major": [v1[(i, j)] for i in range(p) for j in range(p)]})
    s2 = s1 + p                                # wrong residue of s/p
    pis2 = {d: lift(fp[d], s2, p) for d in directions(p)}
    _, ok2 = reconstruct(pis2, s2, p)
    rec["wrong_residue_integral"] = ok2
    return rec


def main():
    trials, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    rng = random.Random(seed)
    rows = []
    for p in (5, 7, 11):
        for _ in range(trials):
            rows.append(trial(p, rng))
    passed = all(r["integral"] and r["character_values_match"] and r["determinant_matches"]
                 and not r["wrong_residue_integral"] for r in rows)
    summary = {"cases": len(rows), "seed": seed, "all_passed": passed,
               "per_prime": {str(p): sum(1 for r in rows if r["p"] == p) for p in (5, 7, 11)}}
    json.dump({"summary": summary, "cases": rows}, open(out, "w"), indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
