#!/usr/bin/env python3
"""Reproduce the five independent statement-review controls and save all I/O.

This runner checks execution only. Interpret the finite controls using the
statement-review report; they do not prove the general theorems.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
CASES = [('core', 'core_controls.py'),
         ('affine', 'affine_and_prime_controls.py'),
         ('prime', 'prime_count_control.m'),
         ('classnumber', 'classnumber_control.m'),
         ('arithmetic', 'arithmetic_controls.py')]
REQUIRED_OUTPUT = {
    'core': ['foundation_controls', "'failures': []", 'GEN1', 'DET1',
             'IR4 37', 'IR4 59', 'IR4 67', 'OBS1 157 1', 'CEN1_below400'],
    'affine': ['affine_constructed_moment_packets 48 failures []',
               'prime_orbit_controls 11 22', 'prime_orbit_controls 13 24',
               'affine_and_prime_control_failures []', 'zero_unit_tau_controls'],
    'prime': ['p11_candidates 945', '77419', '82237', 'p13_first_asserted true'],
    'classnumber': ['11 1', '13 1'],
    'arithmetic': ['GEN1 resultant 77419 prime True',
                   'density comparison 13 40000000', 'bounds 13'],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, default=ROOT/'rerun_outputs')
    args = parser.parse_args()
    output = args.output_root.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error('Choose an unused output directory; existing results are preserved.')
    output.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment['MAGMA_STARTUP_FILE'] = str(ROOT.parent/'inputs/magma_startup.m')
    runtime = {'python': sys.version, 'platform': platform.platform(),
               'python_flint': importlib.metadata.version('python-flint'),
               'sympy': importlib.metadata.version('sympy'),
               'mpmath': importlib.metadata.version('mpmath')}
    results = {}
    for case, script in CASES:
        folder = output/case
        folder.mkdir()
        source = ROOT/script
        if source.suffix == '.m':
            magma = shutil.which('magma')
            if not magma:
                raise SystemExit('Magma executable not found.')
            command = [magma, '-b', str(source)]
        else:
            command = [sys.executable, '-u', str(source)]
        inputs = {script: sha(source)}
        if case == 'core':
            inputs['numerical_statement_addendum.txt'] = sha(ROOT/'numerical_statement_addendum.txt')
        record = {'id': case, 'script': script, 'input_sha256': inputs,
                  'actual_command': command, 'environment': runtime,
                  'magma_startup_sha256': sha(ROOT.parent/'inputs/magma_startup.m'),
                  'started_utc': datetime.now(timezone.utc).isoformat()}
        start = time.monotonic()
        print('START', case, flush=True)
        with (folder/'stdout.log').open('w') as stdout, (folder/'stderr.log').open('w') as stderr:
            result = subprocess.run(command, cwd=folder, env=environment,
                                    stdout=stdout, stderr=stderr)
        logs = (folder/'stdout.log').read_text()+(folder/'stderr.log').read_text()
        errors = re.findall(r'Runtime error[^\n]*|User error[^\n]*|Assertion failed[^\n]*|Syntax error[^\n]*', logs)
        stdout_text = (folder/'stdout.log').read_text()
        missing_output = [marker for marker in REQUIRED_OUTPUT[case] if marker not in stdout_text]
        record.update(returncode=result.returncode, errors=errors,
                      missing_output_markers=missing_output,
                      elapsed_seconds=time.monotonic()-start,
                      finished_utc=datetime.now(timezone.utc).isoformat(),
                      status='execution_passed' if result.returncode == 0 and not errors and not missing_output else 'failed',
                      interpretation='Execution check only; mathematical interpretation is in the independent statement review.')
        record['outputs'] = {name: {'bytes':(folder/name).stat().st_size,
                                    'sha256':sha(folder/name)} for name in ['stdout.log','stderr.log']}
        (folder/'run.json').write_text(json.dumps(record, indent=2)+'\n')
        results[case] = record['status']
        print(record['status'].upper(), case, flush=True)
    summary = {'runner_sha256': sha(Path(__file__)), 'cases': results,
               'all_execution_passed': all(v == 'execution_passed' for v in results.values())}
    (output/'run_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    if not summary['all_execution_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
