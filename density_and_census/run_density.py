#!/usr/bin/env python3
"""Run the Revision 4 computations into OUTPUT_DIR and record a run log.

Usage (from supplementary/density_and_census):
    python run_density.py OUTPUT_DIR [JOBS]
The census to 200000 takes about 15 CPU-minutes with python-flint 0.9.0.
"""
import json
import os
import platform
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = [
    ("census", ["rp_census.py", "200000", "census_200000.json", "census_200000.csv", "{jobs}"]),
    ("model", ["rp_model.py", "12", "model_k12.json"]),
    ("mechanism", ["mechanism_test.py", "20", "20261115", "mechanism_test.json"]),
    ("construction", ["theorem_d_construction.py", "2026", "theorem_d_construction.json"]),
    ("verify", ["check_density_outputs.py", ".", "independent_check.json"]),
]


def main():
    out = os.path.abspath(sys.argv[1])
    jobs = sys.argv[2] if len(sys.argv) > 2 else "2"
    os.makedirs(out, exist_ok=False)
    log = {"python": sys.version, "platform": platform.platform(), "steps": []}
    try:
        import flint, numpy, mpmath, sympy  # noqa: E401
        log["packages"] = {"python-flint": flint.__version__, "numpy": numpy.__version__,
                           "mpmath": mpmath.__version__, "sympy": sympy.__version__}
    except ImportError as exc:
        log["packages"] = f"import failed: {exc}"
    for name, cmd in STEPS:
        args = [sys.executable, os.path.join(HERE, "python", cmd[0])] + [c.format(jobs=jobs) for c in cmd[1:]]
        t0 = time.time()
        r = subprocess.run(args, cwd=out, capture_output=True, text=True)
        open(os.path.join(out, f"{name}.log"), "w").write(r.stdout + r.stderr)
        log["steps"].append({"step": name, "command": args[1:], "returncode": r.returncode,
                             "seconds": round(time.time() - t0, 1)})
        if r.returncode:
            break
    json.dump(log, open(os.path.join(out, "run_log.json"), "w"), indent=1)
    print(json.dumps(log, indent=1))
    sys.exit(max(s["returncode"] for s in log["steps"]))


if __name__ == "__main__":
    main()
