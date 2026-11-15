#!/usr/bin/env python3
"""Exact checks of Proposition 4.5 and Lemma 4.4: the least full tail of S(C_p x C_p) at irregular primes.

The script works directly in R/p^2 R with R = Z[zeta_p], represented as (Z/p^2)[x]/(x^p-1)
followed by reduction modulo Phi_p.  It
  1. computes the irregular half-indices I_p = {i : p | num(B_{2i}), 1<=i<=(p-3)/2};
  2. computes the logarithmic residues of the circular units u_a=(1-zeta^a)/(1-zeta)
     directly from their expansions (independently of any Bernoulli formula) and compares
     them with the Bernoulli formula;
  3. for r = 1, 2, ... builds the transforms of Proposition 4.5,
         H_b = T P(b^(2r+1) T^(2r)),   H_inf = T^(2r+1),
     where P(X) = 1 + X + sum_{nu in U} c_nu X^nu and the c_nu are chosen so that
     [X^nu] log P = 0 for nu in U = {nu >= 2 : r nu in I_p}; it then solves for circular
     units u_b, u_inf with the required residues (impossible only when r itself lies in I_p);
  4. for a solvable r, forms the actual rays of f_b = (1-zeta) u_b and f_inf=(1-zeta)^(2r+1) u_inf
     modulo p from their residues modulo p^2, and checks all p^2 integral reconstruction
     congruences  c_inf(j) + sum_b c_b(i+bj) == s (mod p)  with s = p^2;
  5. records which congruences fail and which nonconstant moment conditions
     (deg h_n <= n and h_{inf,n} = [b^n] h_n, 1 <= n <= p-2) fail.
If the congruences hold, integral reconstruction gives v in Z^G with augmentation s and
these character values, so Theta(v) = s * p^p * p^(2r+1) = p^(p+3+2r).
"""
import sys, json, hashlib
from pathlib import Path
import numpy as np
from flint import fmpq, fmpz

p = int(sys.argv[1]) if len(sys.argv) > 1 else 103
if p < 7 or not fmpz(p).is_prime():
    raise ValueError('p must be a prime at least 7')
if p**4 - p*p > np.iinfo(np.int64).max:
    raise ValueError('p exceeds the proven int64 multiplication bound')
witness_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else None
if witness_dir is not None:
    witness_dir.mkdir(parents=True, exist_ok=True)
q = p * p
M = (p - 3) // 2
m = (p - 1) // 2

# ---------- F_p[T]/(T^(p-1)) helpers (exact, python ints) ----------
N = p - 1                      # precision: coefficients of T^0..T^(p-2)
fact = [1] * p
for i in range(1, p):
    fact[i] = fact[i - 1] * i % p
invfact = [pow(f, p - 2, p) for f in fact]

def ser_mul(A, B, n=N):
    C = [0] * n
    for i, a in enumerate(A[:n]):
        if a:
            for j in range(0, n - i):
                if B[j]:
                    C[i + j] = (C[i + j] + a * B[j]) % p
    return C

def ser_exp_scaled(t, n=N):          # exp(tT) truncated
    return [pow(t, k, p) * invfact[k] % p for k in range(n)]

def ser_log1p(X, n=N):              # log(1+X), X has zero constant term, n <= p-1
    res = [0] * n
    P = [1] + [0] * (n - 1)
    for k in range(1, n):
        P = ser_mul(P, X, n)
        if not any(P):
            break
        c = pow(k, p - 2, p) * (1 if k % 2 == 1 else p - 1) % p
        res = [(r + c * x) % p for r, x in zip(res, P)]
    return res

def normalized_log(A, n=N):
    c0 = A[0] % p
    assert c0 != 0
    inv = pow(c0, p - 2, p)
    X = [a * inv % p for a in A[:n]]
    X[0] = 0
    return ser_log1p(X, n)

# ---------- Bernoulli data ----------
kappa = {}
I_p = []
for i in range(1, M + 1):
    B = fmpq.bernoulli(2 * i)
    num, den = int(B.p), int(B.q)
    if num % p == 0:
        I_p.append(i)
    # kappa_i = B_{2i} / (2i (2i)!) mod p ; denominators prime to p for 2i <= p-3
    d = den * (2 * i) % p * fact[2 * i] % p
    kappa[i] = num % p * pow(d, p - 2, p) % p

# ---------- circular unit residues computed directly ----------
def unit_residue(a):                 # u_a = 1 + zeta + ... + zeta^(a-1) -> sum_k exp(kT)
    S = [0] * N
    for k in range(a):
        E = ser_exp_scaled(k)
        S = [(s + e) % p for s, e in zip(S, E)]
    return S

ulog = {}
formula_ok = True
for a in range(2, m + 1):
    L = normalized_log(unit_residue(a))
    ulog[a] = L
    # compare with (a-1)/2 T + sum kappa_i (a^{2i}-1) T^{2i}
    inv2 = pow(2, p - 2, p)
    expect = [0] * N
    expect[1] = (a - 1) * inv2 % p
    for i in range(1, M + 1):
        expect[2 * i] = kappa[i] * (pow(a, 2 * i, p) - 1) % p
    if L != expect:
        formula_ok = False

# ---------- linear algebra over F_p ----------
def solve_mod_p(A, bvec):
    """Solve A x = b over F_p (A list of rows). Return a solution or None if inconsistent."""
    rows, cols = len(A), len(A[0])
    Mx = [row[:] + [bv] for row, bv in zip(A, bvec)]
    piv_cols = []
    r = 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if Mx[i][c] % p), None)
        if pr is None:
            continue
        Mx[r], Mx[pr] = Mx[pr], Mx[r]
        inv = pow(Mx[r][c], p - 2, p)
        Mx[r] = [x * inv % p for x in Mx[r]]
        for i in range(rows):
            if i != r and Mx[i][c] % p:
                f = Mx[i][c]
                Mx[i] = [(x - f * y) % p for x, y in zip(Mx[i], Mx[r])]
        piv_cols.append(c)
        r += 1
    for i in range(r, rows):
        if Mx[i][cols] % p:
            return None
    x = [0] * cols
    for i, c in enumerate(piv_cols):
        x[c] = Mx[i][cols]
    return x

avec = list(range(2, m + 1))         # 49 or 50 circular units u_2..u_m
def even_rows(maxdeg):
    return [i for i in range(1, M + 1) if 2 * i <= maxdeg]

def lift(target_log, maxdeg):
    """Find exponents n_a (mod p) and k with sum n_a log(u_a) + k T = target through degree maxdeg
    (even coordinates and the linear one). Odd coordinates >=3 of target must vanish."""
    for j in range(3, maxdeg + 1, 2):
        assert target_log[j] == 0, "odd coordinate in target"
    rows = even_rows(maxdeg)
    A = [[ulog[a][2 * i] for a in avec] for i in rows]
    bvec = [target_log[2 * i] for i in rows]
    x = solve_mod_p(A, bvec)
    if x is None:
        return None
    lin = sum(n * ulog[a][1] for n, a in zip(x, avec)) % p
    k = (target_log[1] - lin) % p
    return dict(zip(avec, x)), k

# ---------- arithmetic in (Z/p^2)[x]/(x^p - 1) ----------
def cmul(A, B):
    C = np.zeros(p, dtype=np.int64)
    for i in range(p):
        if A[i]:
            C += A[i] * np.roll(B, i)
            C %= q
    return C

def cpow(A, e):
    R = np.zeros(p, dtype=np.int64); R[0] = 1
    Bp = A.copy()
    while e:
        if e & 1:
            R = cmul(R, Bp)
        Bp = cmul(Bp, Bp)
        e >>= 1
    return R

U = {}
for a in avec:
    v = np.zeros(p, dtype=np.int64); v[:a] = 1
    U[a] = v
g = next(gg for gg in range(2, p) if all(pow(gg, (p - 1) // ell, p) != 1 for ell in {2, 3, 5, 7, 11, 13, 17, 17, 51} if (p - 1) % ell == 0))
# robust primitive root
def is_prim(gg):
    n, f, fs = p - 1, 2, set()
    while f * f <= n:
        while n % f == 0:
            fs.add(f); n //= f
        f += 1
    if n > 1: fs.add(n)
    return all(pow(gg, (p - 1) // ell, p) != 1 for ell in fs)
g = next(gg for gg in range(2, p) if is_prim(gg))
assert g <= m
UGP = cpow(U[g], p)                   # u_g^p == g (mod pR)

def build_unit(expo, k, const_target):
    P = np.zeros(p, dtype=np.int64); P[0] = 1
    for a, n in expo.items():
        if n:
            P = cmul(P, cpow(U[a], n))
    P = np.roll(P, k)                 # times zeta^k
    c = 1
    for a, n in expo.items():
        c = c * pow(a, n, p) % p
    # need c * g^j == const_target (mod p)
    dl = {pow(g, j, p): j for j in range(p - 1)}
    j = dl[const_target * pow(c, p - 2, p) % p]
    P = cmul(P, cpow(UGP, j))
    return P

def to_ray(F, s_mod_q):
    """F in (Z/p^2)[x]/(x^p-1) representing f in R. Return ray mod p with sum s."""
    a = (F[:p - 1] - F[p - 1]) % q          # representative mod Phi_p
    tot = int(a.sum()) % q
    assert tot % p == 0, "element not in the prime above p"
    k = ((s_mod_q - tot) % q) // p % p
    ray = np.zeros(p, dtype=np.int64)
    ray[:p - 1] = (a + k) % p
    ray[p - 1] = k
    return ray

def corrected_P(r, corrected):
    # P(X) = 1 + X + sum_{m>=2, r*m in I_p} c_m X^m with [X^m] log P = 0 for those m (if corrected)
    Pc = {1: 1}
    if not corrected:
        return Pc
    deg = M // r + 1
    for mm in range(2, deg + 1):
        if r * mm in I_p:
            # compute [X^mm] log P with current coefficients; adding c X^mm changes it by c
            ser = [0] * (deg + 1); ser[0] = 1
            for k, c in Pc.items():
                if k <= deg: ser[k] = c % p
            L = ser_log1p([0] + ser[1:], deg + 1) if deg + 1 <= p - 1 else None
            Pc[mm] = (-L[mm]) % p
    if r in I_p:
        return None
    return Pc

def attempt(r, aa=1, corrected=True):
    e = 2 * r + 1
    global Pcoef
    Pcoef = corrected_P(r, corrected)
    if Pcoef is None:
        return {"r": r, "e": e, "solvable": False, "reason": "r itself is an irregular half-index"}
    inv = lambda z: pow(z % p, p - 2, p)
    logE = [0] * N                      # log((exp T - 1)/T) through degree p-2
    Eser = [invfact[k + 1] for k in range(N)]
    logE = normalized_log(Eser)
    result = {"r": r, "e": e}
    # generic targets
    expos = {}
    for b in range(p):
        G = [0] * N; G[0] = 1
        X = aa * pow(b, 2 * r + 1, p) % p
        for mm, cm in Pcoef.items():
            if 2 * r * mm < N:
                G[2 * r * mm] = (G[2 * r * mm] + cm * pow(X, mm, p)) % p
        logG = normalized_log(G)
        target = [(x - y) % p for x, y in zip(logG, logE)]   # log(u_b/u_b(0))
        sol = lift(target, p - 2)
        if sol is None:
            bad = [i for i in I_p if target[2 * i] % p]
            result.update(solvable=False, first_bad_b=b, irregular_rows_nonzero=bad)
            return result
        expos[b] = sol
    # exceptional target: u_inf = -E^{-e} * a  mod T^(p-1-e) ... normalized log = -e log E
    maxdeg = p - 2 - e
    target = [0] * N
    for j in range(1, maxdeg + 1):
        target[j] = (-e * logE[j]) % p
    sol_inf = lift(target, maxdeg)
    if sol_inf is None:
        result.update(solvable=False, exceptional=True)
        return result
    # build actual units mod p^2 and rays
    s_mod_q = 0                        # s = p^2
    one_minus_x = np.zeros(p, dtype=np.int64); one_minus_x[0] = 1; one_minus_x[1] = q - 1
    rays = {}
    for b in range(p):
        expo, k = expos[b]
        ub = build_unit(expo, k, p - 1)          # constant -1
        fb = cmul(one_minus_x, ub)
        rays[b] = to_ray(fb, s_mod_q)
    expo, k = sol_inf
    const = (-1) ** e * aa % p                    # (1-e^T)^e u = (-T E)^e u ; want a T^e
    uinf = build_unit(expo, k, const % p)
    finf = cmul(cpow(one_minus_x, e), uinf)
    rays["inf"] = to_ray(finf, s_mod_q)
    # integral reconstruction congruences
    bad_cells = set()
    for i in range(p):
        for j in range(p):
            tot = int(rays["inf"][j]) + sum(int(rays[b][(i + b * j) % p]) for b in range(p))
            if tot % p:
                bad_cells.add((i, j))
    bad = len(bad_cells)
    ij_pattern = bad_cells == {(i, j) for i in range(1, p) for j in range(1, p)}
    # moment diagnostics: h_n(b) = [T^n] sum_t c_b(t) exp(tT), checked against deg h_n <= n
    # and h_{inf,n} = [b^n] h_n for 1 <= n <= p-2 (the nonconstant moment conditions)
    powt = [[pow(t, n, p) for n in range(p)] for t in range(p)]
    def transform_coeffs(ray):
        return [sum(int(ray[t]) * powt[t][n] for t in range(p)) * invfact[n] % p for n in range(p)]
    hb = [transform_coeffs(rays[b]) for b in range(p)]
    hinf = transform_coeffs(rays["inf"])
    degree_failures, infinity_failures = [], []
    moment_polynomials = {}
    for n in range(1, p - 1):
        vals = [hb[b][n] for b in range(p)]
        coeffs = [vals[0]] + [(-sum(vals[b] * powt[b][p - 1 - k] for b in range(p))) % p
                              for k in range(1, p - 1)] + [(-sum(vals)) % p]
        deg = max((k for k in range(p) if coeffs[k]), default=-1)
        if witness_dir is not None:
            moment_polynomials[str(n)] = coeffs
        if deg > n:
            degree_failures.append([n, deg])
        if hinf[n] != coeffs[n]:
            infinity_failures.append(n)
    digest = hashlib.sha256(b"".join(rays[b].tobytes() for b in list(range(p)) + ["inf"])).hexdigest()
    result.update(solvable=True, congruence_failures=bad, cells=p * p,
                  failures_exactly_at_cells_with_ij_nonzero=ij_pattern,
                  failing_degree_conditions=degree_failures,
                  failing_infinity_conditions=infinity_failures,
                  rays_sha256=digest, determinant_exponent=p + 3 + 2 * r)
    if witness_dir is not None:
        def full_unit_data(exponents, zeta_exponent, target_constant):
            initial_constant = 1
            for a, n in exponents.items():
                initial_constant = initial_constant * pow(a, n, p) % p
            ratio = target_constant * pow(initial_constant, -1, p) % p
            j = next(j for j in range(p-1) if pow(g, j, p) == ratio)
            final_exponents = dict(exponents)
            final_exponents[g] = final_exponents.get(g, 0) + p*j
            return {'circular_unit_exponents': {str(a): int(n) for a,n in exponents.items()},
                    'zeta_exponent': int(zeta_exponent),
                    'target_constant_mod_p': int(target_constant),
                    'constant_correction_j': j, 'constant_correction_exponent_of_u_g': p*j,
                    'final_circular_unit_exponents': {str(a): int(n) for a,n in final_exponents.items()}}
        ray_data = {str(b): [int(x) for x in rays[b]] for b in list(range(p)) + ['inf']}
        witness = {
            'p': p, 'r': r, 'e': e, 'construction_a': aa, 'augmentation': p*p,
            'arithmetic_modulus': q, 'P_coefficients': Pcoef,
            'primitive_root': g, 'circular_unit_indices': avec,
            'generic_unit_data': {
                str(b): full_unit_data(expos[b][0], expos[b][1], p-1)
                for b in range(p)},
            'exceptional_unit_data': full_unit_data(sol_inf[0], sol_inf[1], int(const)),
            'rays_mod_p': ray_data,
            'rays_canonical_json_sha256': hashlib.sha256(json.dumps(ray_data, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            'generic_transform_coefficients': hb,
            'exceptional_transform_coefficients': hinf,
            'moment_polynomials_mod_p': moment_polynomials,
            'failing_cells': [list(cell) for cell in sorted(bad_cells)],
            'result': result.copy(),
            'scope': 'Complete finite reconstruction residues and unit-exponent inputs; no huge full integer vector is asserted.'}
        path = witness_dir / f'p{p}_r{r}.json'
        path.write_text(json.dumps(witness, sort_keys=True, separators=(',', ':')) + '\n')
        result['witness_file'] = path.name
        result['witness_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result

out = {"p": p, "irregular_half_indices": I_p, "bernoulli_formula_matches_direct_unit_logs": formula_ok,
       "primitive_root": g}
admissible = lambda r: all((2 * r + 1) * i > r * (p - 1) for i in I_p if i % r == 0)
out["admissible_r_le_10"] = [r for r in range(1, 11) if admissible(r)]
out["attempts"] = []
rmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
for r in range(1, rmax + 1):
    res = attempt(r, corrected=True)
    res["admissible"] = admissible(r)
    res["P_coefficients"] = {str(k): v for k, v in (Pcoef or {}).items()}
    out["attempts"].append(res)
    print(res, flush=True)
    if res.get("solvable") and res.get("congruence_failures") == 0:
        break
outpath = sys.argv[3] if len(sys.argv) > 3 else None
if outpath:
    json.dump(out, open(outpath, "w"), indent=1)
print(json.dumps({k: out[k] for k in ("p", "irregular_half_indices", "bernoulli_formula_matches_direct_unit_logs", "admissible_r_le_10")}))
