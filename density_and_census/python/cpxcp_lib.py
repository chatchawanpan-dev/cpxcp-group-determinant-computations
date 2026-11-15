"""Exact arithmetic for integer group determinants of C_p x C_p.

Elements of R = Z[zeta_p] are integer lists of length p (coefficients of
zeta^0..zeta^(p-1)); canon() reduces to the basis 1, zeta, ..., zeta^(p-2).
Vectors v in Z^G are dicts {(i, j): int} with (i, j) in F_p^2.
Notation follows the manuscript: f_b = sum v_ij zeta^(i+bj), f_inf = sum
v_ij zeta^j, s = augmentation, (2.12) = lift of a character value to a
projection, (2.5) = integral reconstruction.
"""
from flint import fmpz_mat


def canon(a, p):
    a = list(a) + [0] * (p - len(a))
    top = a[p - 1]
    return [a[t] - top for t in range(p - 1)]


def full(a, p):
    a = list(a) + [0] * (p - len(a))
    return a[:p]


def add(a, b, p):
    a, b = full(a, p), full(b, p)
    return [x + y for x, y in zip(a, b)]


def scal(c, a, p):
    return [c * x for x in full(a, p)]


def mul(a, b, p):
    a, b = full(a, p), full(b, p)
    out = [0] * p
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[(i + j) % p] += x * y
    return out


def one(p):
    return [1] + [0] * (p - 1)


def varpi(p):
    return [1, -1] + [0] * (p - 2)


def varpi_pow(e, p):
    r = one(p)
    for _ in range(e):
        r = mul(r, varpi(p), p)
    return r


def mult_matrix(a, p):
    """Matrix of multiplication by a on the basis 1, ..., zeta^(p-2)."""
    cols = []
    for k in range(p - 1):
        e = [0] * p
        e[k] = 1
        cols.append(canon(mul(a, e, p), p))
    return fmpz_mat([[cols[k][r] for k in range(p - 1)] for r in range(p - 1)])


def norm(a, p):
    return int(mult_matrix(a, p).det())


def is_zero(a, p):
    return all(x == 0 for x in canon(a, p))


def div_by_varpi(a, p):
    """a/(1-zeta) in R; raises if a is not in the prime ideal above p.

    Uses p/(zeta - 1) = sum_k k zeta^k, so 1/(1-zeta) = -(1/p) sum_k k zeta^k.
    """
    w = [k for k in range(p)]
    c = canon(mul(a, w, p), p)
    if any(x % p for x in c):
        raise ValueError("not divisible by varpi")
    return [-x // p for x in c]


def vp_norm(a, p):
    n = abs(norm(a, p))
    if n == 0:
        return None
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def directions(p):
    return list(range(p)) + ["inf"]


def projection(v, d, p):
    pi = [0] * p
    for (i, j), x in v.items():
        t = j % p if d == "inf" else (i + d * j) % p
        pi[t] += x
    return pi


def char_value(v, d, p):
    return projection(v, d, p)  # evaluation sum_t pi(t) zeta^t, as a length-p list


def lift(f, s, p):
    """(2.12): integer vector on F_p with sum s and evaluation f."""
    a = canon(f, p)
    num = s - sum(a)
    if num % p:
        raise ValueError("f not in the prime above p, or p does not divide s")
    k = num // p
    return [x + k for x in a] + [k]


def reconstruct(pis, s, p):
    """(2.5): returns (v, ok) where ok says every entry is integral."""
    v, ok = {}, True
    for i in range(p):
        for j in range(p):
            num = pis["inf"][j] + sum(pis[b][(i + b * j) % p] for b in range(p)) - s
            if num % p:
                ok = False
            v[(i, j)] = num // p if num % p == 0 else None
    return v, ok


def group_det(v, p):
    els = [(i, j) for i in range(p) for j in range(p)]
    M = [[v[((g[0] - h[0]) % p, (g[1] - h[1]) % p)] for h in els] for g in els]
    return int(fmpz_mat(M).det())


def augmentation(v):
    return sum(v.values())


def vp(n, p):
    if n == 0:
        return None
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e
