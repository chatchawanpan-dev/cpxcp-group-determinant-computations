# Census of r_p, the independence model, and checks for Theorem D

These programs produce Table 2 of the manuscript and finite checks of the
construction in the proof of Theorem D (Section 8). They use only Python,
python-flint, NumPy, mpmath and SymPy; no Magma is needed.

## Install and run

From this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_density.py rerun_outputs 2
```

`run_density.py OUTPUT_DIR JOBS` refuses an existing `OUTPUT_DIR`, runs the five
steps below in order, writes each step's output and a `run_log.json` with the
interpreter, package versions, commands, return codes and times. The census
dominates the running time (about 15 CPU-minutes with python-flint 0.9.0).

| Step | Program and arguments | Output |
|---|---|---|
| census | `python/rp_census.py 200000 census_200000.json census_200000.csv JOBS` | irregular pairs and `r_p` for all primes `7 <= p < 200000` |
| model | `python/rp_model.py 12 model_k12.json` | limiting probabilities `P(r_p = k)`, `P(r_p > k)` for `k <= 12` |
| mechanism | `python/mechanism_test.py 20 20261115 mechanism_test.json` | 60 random cases at `p = 5, 7, 11` |
| construction | `python/theorem_d_construction.py 2026 theorem_d_construction.json` | 11 explicit determinants `p^N m` |
| verify | `python/check_density_outputs.py . independent_check.json` | independent checks of all of the above |

`outputs/` holds the packaged run.

## What each program computes

**rp_census.py.** For each prime `p`, the coefficients of
`(x/2) coth(x/2) = cosh(x/2) / (sinh(x/2)/(x/2))` are computed as a power series
in `y = x^2/4` modulo `p` (FLINT series inversion); the coefficient of `y^k`
is a unit multiple of `B_(2k)`, so its vanishing for `1 <= k <= (p-3)/2`
identifies the irregular pairs `(p, 2k)`. Inverse factorials come from Wilson's
theorem. Then `r_p` is the least `r >= 1` with `(2r+1) i > r (p-1)` for every
irregular half-index `i` divisible by `r`. The CSV has one row per prime:
`p; irregular indices 2i; r_p; e(p) - p` (the last column equals `3 + 2 r_p`,
which is `e(p) - p` when `p` does not divide `h_p^+`, in particular for all
`p < 2^31`). This method differs from that of `../python/tail_exponent_census.py`,
which inverts `(e^x - 1)/x`.

**rp_model.py.** In the limit of the independence model the irregular
half-indices divided by `p` form a Poisson process of intensity one on `(0, 1/2)`
with residues equidistributed modulo every integer; `r` is inadmissible when some
point `i` has `r | i` and `i/p <= r/(2r+1)`. Then
`P(r_p > k) = sum over S subset {1..k} of (-1)^|S| exp(-mu(S))`, where `mu(S)` is
the integral over `0 < x < 1/2` of the density of integers divisible by some
`r in S` with `r/(2r+1) >= x`. The `mu(S)` are computed exactly as rationals
(residues modulo `lcm(1..k)`); the exponentials are evaluated with mpmath at 60
digits. Predicted counts in Table 2 are `17981 * P(r_p = k)` and
`17981 * P(r_p > 9)`.

**mechanism_test.py.** Draws `v0` with augmentation `p`, exceptional direction
`infinity` and `v_p(Theta(v0)) = p + 3`; replaces each character value `f_d` by
`f_d (1 + p beta_d)` with random `beta_d` in `Z[zeta_p]`, adds `c varpi^(p-1)` to
`f_inf`, and takes `s'` with `s'/p = tau0 - c (mod p)`. The vector given by (2.12)
and (2.5) must be integral, have the new character values, and satisfy
`Theta(v') = s' prod N(f'_d)` as an exact `p^2 x p^2` determinant (FLINT). As a
control, `s'/p = tau0 - c + 1` must break integrality.

**theorem_d_construction.py.** Carries out the proof of Theorem D(a): from `v0`,
puts `y_b = f_b/varpi` and `y_inf = g_t/varpi^2` with
`g_t = f_inf + (tau0 - t) varpi^(p-1)`, searches `y'_d = y_d + p w` with prime norm
`l_d` (so `y'_d R` is a prime ideal in the ray class of `y_d R` modulo `pR`), and
builds `v` with augmentation `p^a m0` and `Theta(v) = p^N m`, `m = m0 prod l_d`.
Cases: `(p, N, m0)` = (5,8,1), (5,9,2), (7,10,1), (7,10,3), (7,10,-2), (7,11,1),
(7,11,5), (11,14,1), (11,14,7), (11,15,-1), (11,16,2).

**check_density_outputs.py.** Recomputes `r_p` for every census row from the
listed indices and the histogram; recomputes the irregular indices of 726 sampled
primes by inverting `(e^x - 1)/x`, and of all primes below 1500 from exact
Bernoulli numbers; checks `P(r_p = 1) = exp(-1/3)` and
`P(r_p = 2) = exp(-1/5) - exp(-11/30)`; rebuilds the rows of Table 2; and for every
stored vector recomputes the projections, the reconstruction congruences and
`Theta(v) = s prod Res(Phi_p, F_d)` with SymPy resultants, independently of the
FLINT determinant used by the producers.

These are finite checks of stated ranges and cases. They do not prove Theorem D,
whose proof is in Section 8 of the manuscript, and the model is a heuristic.
