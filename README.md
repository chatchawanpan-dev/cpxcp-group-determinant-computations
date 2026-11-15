# Computations for "Bernoulli numbers and the integer group determinants of C_p × C_p"

Chatchawan Panraksa, Applied Mathematics Program, Mahidol University International College

This repository contains the programs, inputs, outputs and verification records for
the computations reported in the manuscript *Bernoulli numbers and the integer
group determinants of C_p × C_p* (2026). The theorems of the manuscript are proved
there; the finite computations here illustrate them and check the numbers quoted in
its tables and examples. They do not prove the general results, and the conditional
use of Vandiver's conjecture is as stated in the manuscript.

## Where each computation appears in the manuscript

| Manuscript | Directory or case |
|---|---|
| Table 1 (irregular primes `p < 400`, `e(p) = p + 3 + 2 r_p`) and the census below 20000 | `python_census_20000`, `python_tail_p37_r1`, …, `python_tail_p157_r2` in `inputs/run_cases.json` |
| Section 9, the obstruction at `p = 103` (degrees 36, 30, 28, 27; the 10404 failing cells) and at `p = 131, 157` | `python_tail_p103_r5`, `python_tail_p131_r2`, `python_tail_p157_r2` |
| Table 2 (`r_p` for the 17981 primes `7 <= p < 200000`, observed and predicted by the independence model) | [`density_and_census/`](density_and_census/README.md) |
| Section 8, finite checks of the construction in the proof of Theorem D | [`density_and_census/`](density_and_census/README.md) (`mechanism_test.py`, `theorem_d_construction.py`) |
| Table 3 (primes `ℓ` with `p^{p+4} ℓ` a determinant, `p = 11, 13`) and the counts 202714 and 31 | `magma_p11_B30000000`, `magma_p13_B12000000`, `magma_p13_B40000000` |
| Equation (9.1), the generator of norm 77419 and its coordinates | `python_generator_certificate` |
| Section 9, the vector `v` with `Θ(v) = 11^15 · 77419` | `magma_vector_p11_77419` |

The 13 cases of `inputs/run_cases.json` are produced by `run_all.py`; their saved runs
are in `reproduced_outputs/`, one directory per case, each with a `run.json` record.
[`python/verify_complete_outputs.py`](python/verify_complete_outputs.py) is an
independent second implementation that checks every saved case; its report is
[`independent_verification.json`](independent_verification.json).
[`COMPUTATIONAL_GUIDE.md`](COMPUTATIONAL_GUIDE.md) describes all cases, bounds and
data formats. `review_controls/` holds five additional finite checks, and
`outputs/` and `historical_scripts/` keep the programs and outputs of the first
version of the computations.

## Quick checks (Python only, a few minutes)

Python 3.13 with the pinned packages:

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
# all 13 saved cases (about two minutes)
.venv/bin/python python/verify_complete_outputs.py --root reproduced_outputs --output my_verification.json
# only the p = 11 vector and its determinant (seconds)
.venv/bin/python python/verify_complete_outputs.py --root reproduced_outputs --case magma_vector_p11_77419 --output vector_check.json
```

The census to 200000 and the model have their own requirements:

```sh
python3.13 -m venv .venv-density
.venv-density/bin/python -m pip install -r density_and_census/requirements.txt
cd density_and_census
../.venv-density/bin/python python/check_density_outputs.py outputs my_check.json
```

## Full reproduction

`run_all.py` reruns the 13 cases into a new directory; four of them need a licensed
Magma (version 2.29-6 was used):

```sh
.venv/bin/python run_all.py --jobs 2 --output-root rerun_outputs
.venv/bin/python python/verify_complete_outputs.py --root rerun_outputs --output rerun_outputs/independent_verification.json
```

`density_and_census/run_density.py` reruns the census, the model and the Theorem D
checks (about 15 CPU-minutes); see its [README](density_and_census/README.md).
Rerunning the nine Python cases and the census in a clean environment on 15 November
2026 reproduced every saved output file byte for byte, apart from timing fields.

## File integrity

[`MANIFEST.json`](MANIFEST.json) lists the size and SHA-256 hash of every other file.
To check them:

```sh
python3 -c "import json,hashlib; m=json.load(open('MANIFEST.json')); bad=[f['path'] for f in m['files'] if hashlib.sha256(open(f['path'],'rb').read()).hexdigest()!=f['sha256']]; print('all files match' if not bad else bad)"
```

In the saved run records (`run.json`, `run_log.json`, `independent_verification.json`)
and in `review_controls/historical_review/statement_core_executed.py`, absolute local
paths have been replaced by the placeholders `<repository>`, `<venv>` and
`<review-packets>`; nothing else in those files was changed. Docstrings of the
older programs refer to section numbers of an earlier version of the manuscript.

## License

MIT; see [LICENSE](LICENSE).
