# Computational guide

This guide describes the finite inputs and certificate formats of the computations
for the manuscript. Run commands from the repository root.
`inputs/run_cases.json` is the case manifest; its 13 case identifiers also name the
subdirectories under the selected output root. No private reference-book file is
needed for execution or included in this repository.

## Installation and observed environment

The fresh Python run records identify Python 3.13.5 on macOS ARM64,
`python-flint==0.9.0`, and `numpy==2.4.4`. `requirements.txt` pins the two packages.
Python 3.13 is the tested interpreter series; this guide does not assert tests on
other operating systems or Python versions. Magma is installed and licensed
separately; the fresh Magma logs identify version 2.29-6. Place its `magma`
executable on `PATH` before running Magma cases.

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -c "import sys, flint, numpy; print(sys.version); print(flint.__version__); print(numpy.__version__)"
command -v magma
```

The original November 5 README reported Python 3.13.16, python-flint 0.9.0 and
NumPy 2.5.3 on Linux, and Magma V2.29-6 on macOS. These are historical declarations;
the old logs lack session-version headers. They are not the observed versions of
the fresh reproduction. The separate compact-generator certificate originally
saved on November 8 records Python 3.9.6; its new run records the new interpreter.

The mathematical identity and decision computations use integers, rationals, or
finite rings and fields. Decimal heuristic fractions, predicted counts, printed
bounds, and elapsed times are summaries and are not exact arithmetic certificates.
The tail producer uses NumPy `int64` for modular products and rejects a prime when
its implemented bound `p**4 - p*p <= int64_max` does not hold.

## Run and resume

For a clean reproduction of every manifest case, choose an unused output root:

```sh
.venv/bin/python run_all.py --jobs 2 --output-root rerun_outputs
.venv/bin/python python/verify_complete_outputs.py --root rerun_outputs --output rerun_outputs/independent_verification.json
```

The runner's default is `--jobs 2 --output-root reproduced_outputs`. To reuse
completed matching cases in that default directory, run:

```sh
.venv/bin/python run_all.py --jobs 2 --output-root reproduced_outputs --resume
.venv/bin/python python/verify_complete_outputs.py --root reproduced_outputs --output independent_verification.json
```

`--case ID` may be repeated to select cases; omitting it selects all 13. For example:

```sh
.venv/bin/python run_all.py --jobs 2 --output-root selected_rerun --case python_generator_certificate --case magma_vector_p11_77419
.venv/bin/python python/verify_complete_outputs.py --root selected_rerun --case python_generator_certificate --case magma_vector_p11_77419 --output selected_rerun/independent_verification.json
```

Each producer runs in its own case directory. The Python executable is the one
used to launch `run_all.py`; Magma is resolved from `PATH`. Every Magma process
receives `MAGMA_STARTUP_FILE` pointing to the supplied empty
`inputs/magma_startup.m`. Each Magma program calls `SetSeed(1)` and prints its
version and numerical inputs. The Python algorithms have no random sampling.
This fixes the recorded choices in the tested environment; it is not a promise
of identical Magma generator choices across software versions or platforms.

The runner writes `stdout.log`, `stderr.log`, and `run.json` for each case. It
accepts a run only when the process exits zero, its log contains none of the
implemented error markers, and all manifest `expected_files` exist. A batch
`run_summary.json` is written after all selected futures finish. A batch with a
failed case exits with code 1. Inspect the full logs when diagnosing a failure;
the error-marker scan is an execution check, not a mathematical verifier.

`--resume` requires a prior `passed` record, a matching script hash, exact matching
case input, matching startup-file hash, and matching hashes for the recorded
outputs. It reuses that recorded run; it does not rerun it under the current
interpreter or certify that the current software environment is identical. An
existing nonempty case directory that fails these checks is refused. Preserve
it and choose a new output root for a clean rerun. Do not use the same output root
for two concurrent runner invocations.

## The 13 finite cases

Bounds on primes are inclusive in the producer manifest. The census cutoff
20000 is composite, so `7 <= p <= 20000` and `7 <= p < 20000` select the same primes.
Each tail case tries `r=1,...,RMAX`, stopping at the first solvable attempt with
zero reconstruction failures. It records failed attempts before that success.

| Case identifier | Producer and exact input | Certificate files beyond logs and run record |
|---|---|---|
| `magma_p13_B40000000` | `magma/defect_distribution_p13.m`, `B:=40000000 records:=decisions.csv`; every prime `ell <= 40000000` with `ell mod 13 = 1` | `decisions.csv` |
| `magma_p11_B30000000` | `magma/prime_cofactor_counts.m`, `p:=11 B:=30000000 records:=decisions.csv`; every prime `ell <= 30000000` with `ell mod 11 = 1` | `decisions.csv` |
| `magma_p13_B12000000` | `magma/defect_distribution_p13.m`, `B:=12000000 records:=decisions.csv`; every prime `ell <= 12000000` with `ell mod 13 = 1` | `decisions.csv` |
| `python_census_20000` | `python/tail_exponent_census.py 20000 2000 census.json census_all_primes.csv`; all primes `7 <= p <= 20000`, producer exact-Bernoulli cross-check for primes through 2000 | `census.json`, `census_all_primes.csv` |
| `python_tail_p103_r5` | `python/check_bernoulli_tail.py 103 5 result.json witnesses` | `result.json`; `witnesses/p103_r1.json` through `p103_r5.json` |
| `python_tail_p157_r2` | `python/check_bernoulli_tail.py 157 2 result.json witnesses` | `result.json`; `witnesses/p157_r1.json`, `p157_r2.json` |
| `python_tail_p131_r2` | `python/check_bernoulli_tail.py 131 2 result.json witnesses` | `result.json`; `witnesses/p131_r1.json`, `p131_r2.json` |
| `python_tail_p37_r1` | `python/check_bernoulli_tail.py 37 1 result.json witnesses` | `result.json`; `witnesses/p37_r1.json` |
| `python_tail_p59_r1` | `python/check_bernoulli_tail.py 59 1 result.json witnesses` | `result.json`; `witnesses/p59_r1.json` |
| `python_tail_p67_r1` | `python/check_bernoulli_tail.py 67 1 result.json witnesses` | `result.json`; `witnesses/p67_r1.json` |
| `python_tail_p101_r1` | `python/check_bernoulli_tail.py 101 1 result.json witnesses` | `result.json`; `witnesses/p101_r1.json` |
| `magma_vector_p11_77419` | `magma/determinant_example_p11.m`, `ell:=77419`; `p=11` fixed in source | `determinant_example_p11_vector.txt`, `determinant_example_p11_matrix_entries.csv` |
| `python_generator_certificate` | `python/check_generator_certificate.py certificate.json`; `p=11` and coefficients `[-2,-2,-3,-2,-1,-2,-1,-3,0,0]` fixed in source | `certificate.json` |

The manifest records the engine, script, argument list and expected files for
every case. Tail cases additionally record `prime` and `rmax`; prime scans record
`prime`, `bound`, seed, ideal-selection rule and class-number proof setting. The
vector case records the group-matrix convention; the compact-generator case
records its coefficient list. Arbitrary changes to these ranges are new cases;
the independent verifier is designed for the listed 13 cases and their bounds.

## Prime-decision records

Both scan programs work in `K = Q(zeta_p)` with its ring of integers `O`. Before
the scan they call
`ClassNumber(O : Proof := "Full", Bound := MinkowskiBound(O))`, assert the result
is 1, and print it. This is an explicit full class-number calculation in the
producer, separate from assertions that each scanned ideal is principal. The
Python verifier does not independently recalculate a class group.

For every split rational prime in the specified range, the producer chooses
`Decomposition(O,ell)[1][1]`, obtains a generator with `IsPrincipal`, and checks
`Abs(Norm(y)) == ell`. Every such prime, including every rejected prime, has one
row in `decisions.csv`, in increasing order of `ell`. The delimiter is a semicolon:

| Column | Meaning |
|---|---|
| `ell` | Rational prime being tested |
| `norm` | Exact absolute norm of the recorded generator, asserted equal to `ell` |
| `power_basis_coefficients` | List `[a_0,...,a_(p-2)]` representing `y = sum a_j zeta_p^j` |
| `higher_odd_coordinates` | List `[c_3,c_5,...,c_(p-2)]`, each represented by an integer from 0 to `p-1` |
| `admitted` | Integer 1 if all coordinates except the last two vanish and at least one of the last two is nonzero; otherwise 0 |

List-valued fields use bracketed integer-list syntax. To recompute the
coordinates, substitute `zeta_p = exp(T)` through degree `p-2`, divide the resulting
series by its constant term, and take its formal logarithm in
`F_p[T]/(T^(p-1))`. The verifier uses differentiation, series inversion and
integration rather than the producer's logarithm power sum. It checks all rows
against a separate prime sieve and an integer resultant norm.

The p=11 log gives membership counts at `10^5`, `10^6`, `10^7`, and the final
bound, a complete admitted-prime list and the number of split primes processed.
The p=13 logs also give the full 2197-cell histogram of `(c_3,c_5,c_7)`, encoded
by `c_3 + 13*c_5 + 169*c_7`, and membership checkpoints at `10^7` and `B`.
Individual nonzero tuples and their frequencies depend on the selected prime
ideal and the fixed coordinate/embedding conventions. The manuscript retains the invariant zero-triple
and membership counts. The zero-triple count alone omits the last-two-coordinates
condition; equality of those two sets in these finite datasets is checked, not
assumed as a general identity. The verifier compares admission sets on the overlap
of the two p=13 ranges.

## Census schema and Vandiver scope

`census_all_primes.csv` is semicolon-delimited and contains every prime in the
census, including regular primes. Its columns are:

| Column | Meaning |
|---|---|
| `p` | Prime with `7 <= p <= 20000` |
| `irregular_indices` | JSON list of even indices `k=2,4,...,p-3` for which `p` divides the numerator of `B_k`; empty for a regular prime |
| `r_p` | Least positive `r` satisfying `(2r+1)*i > r*(p-1)` for every irregular half-index `i=k/2` divisible by `r` |
| `p_plus_3_plus_2r_p` | Exact integer `p+3+2*r_p` |
| `within_cited_vandiver_range` | Integer flag `int(p < 2**31)`; metadata about the manuscript's cited external range |

`census.json` has `summary` and `irregular_primes`. Each irregular-prime entry
contains `p`, `irregular_indices`, `r_p`, and `tail_exponent`. The summary records
`PMAX`, prime and irregular-prime counts, the `r_p_histogram`, maximum `r_p`,
smallest prime with `r_p > 1`, maximum index of irregularity, the boolean
`max_r_p_le_(p-5)/2_for_all`, and two decimal comparisons
(`fraction_r_p_gt_1`, `heuristic_1_minus_exp(-1/3)`).

The producer obtains irregular indices by modular power-series inversion and
cross-checks them against exact Bernoulli numbers only through its `CHECK=2000`
argument. The independent verifier instead uses exact rational Bernoulli
numerators for the entire listed census, checking every row and the summaries.
The census has 2259 prime rows and records 900 irregular primes; these are
reference values to test, not a substitute for inspecting verification status.

The `within_cited_vandiver_range` flag does not test Vandiver's conjecture. The
manuscript cites Hart–Harvey–Ong for primes below `2^31`; this repository does not
reproduce that external computation. Interpretation of `p+3+2*r_p` as the stated
least full-tail exponent uses the manuscript's hypotheses, including its
Vandiver condition. Outside the externally supported range, the unconditional
construction bound and the conditional equality must retain their separate
claim strengths. Finite census observations do not prove a distribution law.

## Tail attempts and residue witnesses

Each `result.json` contains `p`, `irregular_half_indices`, `primitive_root`,
`bernoulli_formula_matches_direct_unit_logs`, `admissible_r_le_10`, and `attempts`.
An attempt records `r`, `e=2*r+1`, whether its unit system is `solvable`, the
admissibility flag and the correction polynomial `P_coefficients`. A solvable
attempt also records the number of failed reconstruction cells out of `p*p`,
whether its failure set is exactly the cells with both coordinates nonzero,
`failing_degree_conditions` as `[n,degree]` pairs,
`failing_infinity_conditions` as indices, a ray digest, the proposed
`determinant_exponent=p+3+2*r`, and witness filename and hash. A nonsolvable
attempt instead records its reason or failed unit equations. The manifest's
listed witness files are expected for the fixed cases above; arbitrary parameters
may fail before producing such a witness.

The complete finite residue certificate `witnesses/pP_rR.json` has these groups
of fields:

| Fields | Representation and meaning |
|---|---|
| `p`, `r`, `e`, `construction_a`, `augmentation`, `arithmetic_modulus` | Construction parameters; the supplied cases use `construction_a=1`, augmentation `p*p` and modulus `p*p` |
| `P_coefficients` | JSON object mapping degree strings to coefficients modulo `p` |
| `primitive_root`, `circular_unit_indices` | Chosen primitive root and circular-unit indices used in the lift |
| `generic_unit_data` | Object keyed by strings `"0"`,...,`"p-1"`, giving the lift data for each finite direction |
| `exceptional_unit_data` | Lift data for the direction at infinity |
| `rays_mod_p` | Object keyed by the finite direction strings and `"inf"`; each value lists the `p` residues in order `t=0,...,p-1` |
| `rays_canonical_json_sha256` | SHA-256 of the ray object serialized with sorted keys and separators `(',', ':')`, encoded as UTF-8 |
| `generic_transform_coefficients` | `p` rows, indexed by finite direction, each containing coefficients of degrees `0,...,p-1` |
| `exceptional_transform_coefficients` | Coefficients of degrees `0,...,p-1` for the infinity direction |
| `moment_polynomials_mod_p` | Object keyed by `n=1,...,p-2`; each length-`p` list gives the coefficients in degrees `0,...,p-1` of the canonical polynomial in the direction variable |
| `failing_cells` | Sorted list of all zero-based pairs `[i,j]` failing the reconstruction congruence |
| `result`, `scope` | A snapshot of the attempt's residue diagnostics and an explicit scope statement |

Each unit-data object records `circular_unit_exponents`, `zeta_exponent`,
`target_constant_mod_p`, `constant_correction_j`,
`constant_correction_exponent_of_u_g`, and `final_circular_unit_exponents`.
The final exponent of the primitive-root circular unit includes the correction
`p*j`; all circular-unit maps use stringified integer indices. These inputs permit
the verifier to reconstruct the units modulo `p^2`, derive the rays, and check
every one of the `p^2` reconstruction cells. It independently derives all moment
transforms and canonical moment polynomials from those rays.

`result.rays_sha256` is the older digest of concatenated NumPy array bytes and
can depend on byte order. Use `rays_canonical_json_sha256` for portable comparison;
`witness_sha256` checks the exact saved witness bytes. The witness's nested
`result` is captured before the outer attempt acquires its witness-link fields.
These are complete residue certificates for the tested constructions, not huge
integer coefficient vectors or direct full-matrix determinants for those primes.

## Explicit matrix and compact generator

`determinant_example_p11_vector.txt` is a nested list of 11 rows of 11 integers,
representing `v[i][j]` for `0 <= i,j < 11`. Magma may wrap long integers with a
backslash followed by a newline; join such continuations before parsing.
The group order is `G=[(0,0),(0,1),...,(10,10)]`, with the first coordinate major.

`determinant_example_p11_matrix_entries.csv` is semicolon-delimited, with header
`row_index;column_index;value`. Indices are **zero-based**, from 0 through 120,
and all `121^2=14641` entries are written. If a row is indexed by `w=(i,j)` and
a column by `x=(k,l)`, the entry is `v[(k-i) mod 11][(l-j) mod 11]`.
The producer computes the exact integer determinant and prints the generator
coefficients, norm, higher odd coordinates, augmentation and character norms.
The independent verifier reconstructs every matrix entry from the vector and
compares a direct FLINT integer determinant with a separate product of exact
character norms. The target is `11^15*77419 = 323398376027990284769`.

`certificate.json` has a `python` version field and a `certificate` object with
`p`, `power_basis_coefficients`, `norm`, `constant_residue`,
`odd_logarithmic_coordinates`, `historical_example_galois_multiplier`, and
`checks`. The compact generator is represented by
`[-2,-2,-3,-2,-1,-2,-1,-3,0,0]`, in ascending powers of `zeta_11`. The expected
norm is 77419, constant residue is 6, and coordinates in degrees 3,5,7,9 are
`[0,0,10,9]`. The recorded multiplier 3 relates them to the historical
`[0,0,6,5]`. This small certificate checks that generator; the separate matrix
case supplies the complete integer vector and matrix.

## Independent verification and status records

`python/verify_complete_outputs.py` reads the recorded outputs without importing
or executing a producer. It independently checks prime coverage and all decision
rows, the full census, every tail witness cell and moment polynomial, every
matrix entry and its determinant, and the compact generator. It shares FLINT
with some producer calculations, so algorithmic independence is not an
independent arithmetic-library implementation. The current version and the
precise checks performed are recorded in its output report.

The verifier reports `passed_range` (exit 0), `failed` (exit 1), or `incomplete`
(exit 2). Missing files and cases whose `run.json` is not `passed` are incomplete.
A mismatching datum is failed. Its report contains its own hash, selected cases,
input-root path, command, environment, case-manifest hash, optional manuscript
snapshot hash, source locators, limitations and per-case input hashes, mismatches,
missing items, timing and computed results. A partial selection verifies only
that selection. Consult the report's actual case list and statuses.

Each producer `run.json` contains `id`, `input`, `script_sha256`,
`actual_command`, start/end UTC timestamps, environment information,
`returncode`, elapsed seconds, detected errors, missing expected files and
`outputs`. The latter maps relative output paths to byte lengths and SHA-256
hashes; `run.json` excludes itself from that map. A running case lacks final
completion fields. `run_summary.json` identifies the manifest hash, selected
cases, per-case execution outcomes and `all_selected_passed`. A summary is
specific to that invocation and may cover only a subset.

The revision's current independent report is
[independent_verification.json](independent_verification.json).
Use it and the per-case records for current status. This guide specifies how to
reproduce and interpret the data; it does not assert completion of a run still
in progress, resolve source-support gaps, or certify a general proof.

## Historical evidence and comparison

`outputs/python/` preserves the seven original tail JSON/log pairs and census
JSON/output pair, plus the separately identified compact-generator certificate.
`outputs/magma/` preserves the original p=11 scan, two p=13 scans, and explicit
example log/vector: 21 original output files in total before the added compact
certificate. The current programs add certificates and execution records, so
the old files are not expected to match fresh outputs byte for byte. Exact copies of the five original programs are included in
`historical_scripts/python/` and `historical_scripts/magma/`, extracted from the
November 5 version of the supplement.

Compare exact mathematical records using their stated conventions and invariant
counts. Timing, environment paths, JSON formatting, generator choices and
nonzero-tuple frequencies can differ without establishing a mathematical
contradiction. Conversely, a retained finite claim with a genuine mismatch
requires diagnosis; a successful build or a previous favorable report does not
resolve it. `MANIFEST.json` lists the SHA-256 fingerprint of every file in this
repository; the README shows how to check them.

## Additional statement-review controls

A separate
[review-control directory](review_controls/README.md) holds additional finite checks. It contains the five
reviewer-derived programs, the numerical statement addendum, execution metadata
and all saved standard output and error logs. Its requirements are
`python-flint==0.9.0`, `sympy==1.14.0` and `mpmath==1.3.0`; the Magma controls
use separately licensed Magma 2.29-6. From the repository root, a clean
run into an unused output root is:

```sh
python3.13 -m venv .venv-review
.venv-review/bin/python -m pip install -r review_controls/requirements.txt
command -v magma
.venv-review/bin/python review_controls/run_controls.py --output-root review_controls/rerun_outputs
```

This runner executes all five cases sequentially and refuses a nonempty output
root. It has no resume option. Its `execution_passed` status checks process
completion and known error markers; it does not automatically adjudicate all
printed mathematical results. The linked README gives each exact finite scope,
the observed results, the original-source provenance and the statement-review
interpretation. In particular, the quick class-number calls use Magma defaults;
they are not independent unconditional certificates. The full Minkowski-bound
producer records described above and the cited Washington tables provide the
separate unconditional support.

The five controls provide additional bounded checks; they do not extend the
tested ranges of the 13 cases.

## Editorial notation and provenance

The manuscript uses `I_p^u` and `r^u` for the indices missing from the global-unit image and their least admissible integer; older programs, review addenda and saved logs may use `I_p^E` and `r^E` for the same objects. The exceptional character-value valuation is now written `e_infty`, distinct from the tail exponent `e(p)`. These changes do not alter any algorithm, input, computed coordinate or output. Original filenames containing `revision`, dated execution records and historical review statements retain their original wording.

## Census to 200000, model and Theorem D checks

The directory `density_and_census/` is self-contained and separate from the 13
cases above; `run_all.py` and `python/verify_complete_outputs.py` do not use it.
It has its own pinned requirements (python-flint 0.9.0, NumPy 2.5.3, mpmath
1.3.0, SymPy 1.14.0):

```sh
python3.13 -m venv .venv-density
.venv-density/bin/python -m pip install -r density_and_census/requirements.txt
.venv-density/bin/python density_and_census/run_density.py rerun_density 2
```

The runner refuses an existing output directory and writes, for each of its five
steps, the output files, a log and `run_log.json`. The packaged run is
`density_and_census/outputs/`; its `independent_check.json` records the checks
listed in `density_and_census/README.md`. The census CSV has one semicolon-separated
row per prime `p` with `7 <= p < 200000`: `p`, the irregular indices `2i`
separated by spaces, `r_p`, and `3 + 2 r_p`.
