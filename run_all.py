#!/usr/bin/env python3
"""Run every stated experiment into isolated, hash-identified result directories.

Run with the Python interpreter containing requirements.txt dependencies.
Magma must be installed separately. No historical output is overwritten.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_case(case, output_root, resume):
    folder = output_root / case['id']
    source = (ROOT / case['script']).resolve()
    source.relative_to(ROOT)
    metadata_path = folder / 'run.json'
    if folder.exists() and any(folder.iterdir()):
        if resume and metadata_path.exists():
            old = json.loads(metadata_path.read_text())
            if (old.get('status') == 'passed' and old.get('script_sha256') == sha(source)
                    and old.get('input') == case
                    and old.get('environment', {}).get('magma_startup_sha256')
                        == sha(ROOT/'inputs/magma_startup.m')):
                if all((folder/name).is_file() and sha(folder/name) == row['sha256']
                       for name,row in old['outputs'].items()):
                    print('REUSED', case['id'], flush=True)
                    return old
        raise RuntimeError(f'{folder} exists without a matching completed run; choose a fresh output root.')
    folder.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment['MAGMA_STARTUP_FILE'] = str(ROOT/'inputs/magma_startup.m')
    if case['engine'] == 'python':
        command = [sys.executable, '-u', str(source)] + case['arguments']
    else:
        magma = shutil.which('magma')
        if not magma:
            raise RuntimeError('Magma executable was not found.')
        command = [magma, '-b'] + case['arguments'] + [str(source)]
    metadata = {'id': case['id'], 'status':'running', 'input':case,
                'script_sha256':sha(source), 'actual_command':command,
                'started_utc':datetime.now(timezone.utc).isoformat(),
                'environment': {'python':sys.version, 'platform':platform.platform(),
                    'machine':platform.machine(), 'python_flint':importlib.metadata.version('python-flint'),
                    'numpy':importlib.metadata.version('numpy'),
                    'magma_startup_sha256':sha(ROOT/'inputs/magma_startup.m')}}
    metadata_path.write_text(json.dumps(metadata,indent=2)+'\n')
    print('START',case['id'],flush=True)
    start = time.monotonic()
    with (folder/'stdout.log').open('w') as stdout, (folder/'stderr.log').open('w') as stderr:
        process = subprocess.run(command,cwd=folder,env=environment,stdout=stdout,stderr=stderr)
    log = (folder/'stdout.log').read_text() + (folder/'stderr.log').read_text()
    errors = re.findall(r"Runtime error[^\n]*|User error[^\n]*|Assertion failed[^\n]*|Can't open file[^\n]*|Syntax error[^\n]*",log)
    missing = [name for name in case['expected_files'] if not (folder/name).is_file()]
    metadata.update({'returncode':process.returncode,'elapsed_seconds':time.monotonic()-start,
                     'finished_utc':datetime.now(timezone.utc).isoformat(),
                     'errors':errors,'missing_files':missing,
                     'status':'passed' if process.returncode==0 and not errors and not missing else 'failed'})
    metadata['outputs'] = {str(f.relative_to(folder)):{'bytes':f.stat().st_size,'sha256':sha(f)}
                           for f in sorted(folder.rglob('*')) if f.is_file() and f!=metadata_path}
    metadata_path.write_text(json.dumps(metadata,indent=2)+'\n')
    print(metadata['status'].upper(),case['id'],f"{metadata['elapsed_seconds']:.1f}s",flush=True)
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path,default=ROOT/'reproduced_outputs')
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--case',action='append',default=[])
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    config=json.loads((ROOT/'inputs/run_cases.json').read_text())
    cases=[case for case in config['cases'] if not args.case or case['id'] in args.case]
    unknown=set(args.case)-{c['id'] for c in config['cases']}
    if unknown:
        parser.error('Unknown case IDs: '+', '.join(sorted(unknown)))
    if not cases or args.jobs<1:
        parser.error('Select at least one valid case and a positive job count.')
    output=args.output_root.resolve();output.mkdir(parents=True,exist_ok=True)
    results=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures={pool.submit(run_case,c,output,args.resume):c['id'] for c in cases}
        for future in as_completed(futures):
            try:results.append(future.result())
            except Exception as error:
                results.append({'id':futures[future],'status':'failed','exception':str(error)})
                print('FAILED',futures[future],str(error),flush=True)
    summary={'manifest_sha256':sha(ROOT/'inputs/run_cases.json'),
             'selected_cases':[c['id'] for c in cases],
             'results':{r['id']:{k:v for k,v in r.items() if k in ['status','elapsed_seconds','errors','missing_files','exception']} for r in results},
             'all_selected_passed':all(r['status']=='passed' for r in results)}
    (output/'run_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not summary['all_selected_passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
