# Additional computational controls

These five programs provide additional finite checks of reconstruction, affine
conditions, prime coordinates, class numbers and exact arithmetic. Their programs,
explicit inputs and complete saved run logs are included below.

Original program snapshots are preserved in
[historical_review](historical_review/). The [provenance record](PROVENANCE.json) describes
the small changes needed for these portable copies. The core program now locates
its [numerical addendum](numerical_statement_addendum.txt) beside the script. The
Magma programs prepend seed 1, a column width of 1024 and a version header. The
affine and arithmetic programs are unchanged copies. All other numerical inputs
are explicit in the programs; Python sampling uses seed `20261108`.

## Reproduction

Run these commands from the repository root. Python 3.13.5 on
macOS ARM64 was observed in the saved run records. The separate requirements pin
`python-flint==0.9.0`, `sympy==1.14.0` and `mpmath==1.3.0`. Magma 2.29-6
(version 2.29.6) must be installed and licensed separately, with `magma` on `PATH`.
No source PDF or private reference book is required.

```sh
python3.13 -m venv .venv-review
.venv-review/bin/python -m pip install -r review_controls/requirements.txt
.venv-review/bin/python -c "import sys, flint, sympy, mpmath; print(sys.version); print(flint.__version__, sympy.__version__, mpmath.__version__)"
command -v magma
.venv-review/bin/python review_controls/run_controls.py --output-root review_controls/rerun_outputs
```

Choose a new output directory. The runner refuses a nonempty output root and has
no resume option. It runs all five programs sequentially, using the Python
interpreter that launches it and the supplied `../inputs/magma_startup.m` for
Magma startup. Preserve the packaged `reproduced_outputs` directory; use another
unused name for each subsequent run.

## Scope and saved results

| Case | Program and finite scope | Saved result |
|---|---|---|
| `core` | `core_controls.py`: 80 transform controls, 48 projection/recovery controls and 96 positive/perturbed moment controls at p=5,7,11; the supplied 121 by 121 integer matrix; selected reconstruction obstructions at p=37,59,67,103,131,157; all 2259 primes with 7 <= p < 20000 | No foundation-control failures; determinant 323398376027990284769; 900 irregular primes and 657 primes with r > 1; the obstruction-cell sets match the reported controls |
| `affine` | `affine_and_prime_controls.py`: 48 constructed moment packets at p=7,11,13, 46 prime-coordinate vectors at p=11,13 with permitted conjugates, and zero-unit controls | No recorded affine or prime-orbit failures |
| `prime` | `prime_count_control.m`: all 945 split primes for p=11 through 100000, plus the individual p=13 candidate 6285241 | p=11 hits 77419 and 82237; the specified p=13 candidate has coordinates [0,0,0,11,8] |
| `classnumber` | `classnumber_control.m`: Magma’s default `ClassNumber` calculation for the cyclotomic fields at p=11,13 | Both calls return 1 |
| `arithmetic` | `arithmetic_controls.py`: an independent SymPy resultant and primality check for 77419, six printed prime-count comparisons, and the two branch bounds at p=11,13 | Resultant 77419, primality true, and the recorded numerical comparisons |

Each case directory under [reproduced_outputs](reproduced_outputs/) contains
`stdout.log`, `stderr.log` and `run.json`. The run records contain input hashes,
actual commands, environment, timestamps, elapsed time, return code, detected
errors and output hashes. The [batch summary](reproduced_outputs/run_summary.json)
records `execution_passed` for all five cases. The complete saved outputs are included in each case directory.

`execution_passed` records a zero exit code, no detected error marker, and all
required output markers. Inspect the complete diagnostic output for the finite
mathematical results; these checks do not establish general theorems.

The class-number quick control uses Magma defaults and is not an independent
unconditional certificate. The separate manuscript producers explicitly use
`ClassNumber(O : Proof := "Full", Bound := MinkowskiBound(O))`; their preserved
records and the cited Washington tables supply the unconditional support described
in the computational guide. This small prime control does not independently repeat the
large scans or prove p=13 minimality.

The historical snapshots include an obsolete initial affine diagnostic that fixes
one ideal rather than allowing the conjugates in the manuscript. It is retained
for provenance and is not used by `run_controls.py`. Use the five programs listed
in the table and the supplied runner for the current controls.

The original 13 manuscript experiments, their complete data and their independent
verification are in the main directories of this repository. Their dated execution records are preserved. See
[the computational guide](../COMPUTATIONAL_GUIDE.md) and
[the independent verification report](../independent_verification.json) for their
larger ranges and distinct verification status. These controls supplement the larger recorded experiments.

The runner also requires the expected opening/result/completion markers in each standard-output log. A zero exit code with an empty or incomplete log is rejected. These completeness checks supplement the mathematical interpretation described above.
