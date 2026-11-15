#!/usr/bin/env python3
"""Independent, exact verifier of the complete finite reproduction certificates.

Run from any directory with Python and python-flint installed::

    python python/verify_complete_outputs.py --output independent_verification.json

The default input root is this supplement's ``reproduced_outputs`` directory.
No author program is imported or executed. Missing or unfinished cases give
``incomplete`` and exit 2; any mathematical/data mismatch gives ``failed`` and
exit 1. ``passed_range`` means only the explicitly enumerated finite checks.
"""
import argparse
import ast
from collections import Counter
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import re
import sys
import time

try:
    import flint
    from flint import fmpq, fmpz_mat, fmpz_poly, nmod_mat, nmod_poly
except ImportError:
    flint = None

SUPPLEMENT = Path(__file__).resolve().parents[1]
PRIME_CASES = [(11, 30000000), (13, 12000000), (13, 40000000)]
TAIL_CASES = [(37, 1), (59, 1), (67, 1), (101, 1), (103, 5), (131, 2), (157, 2)]
CASE_IDS = ([f'magma_p{p}_B{b}' for p, b in PRIME_CASES]
            + ['python_census_20000']
            + [f'python_tail_p{p}_r{r}' for p, r in TAIL_CASES]
            + ['magma_vector_p11_77419', 'python_generator_certificate'])
SOURCE_LOCATORS = {
    'PC1': 'manuscript/main.tex: Theorem C(a), Table tab:primes, par:invariant-counts',
    'CEN1': 'manuscript/main.tex: definition of admissible r, Table tab:tail, census paragraph',
    'OBS1': 'manuscript/main.tex: eq:reconstruction-congruence, eq:degree, eq:constant, eq:tail-data',
    'DET1': 'manuscript/main.tex: recorded p=11 determinant example',
    'GEN1': 'manuscript/main.tex: eq:generator-certificate and eq:defects',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def integer_list(text):
    value = ast.literal_eval(text.replace('\\\n', '').replace('\\\r\n', ''))
    if not isinstance(value, list) or any(type(x) is not int for x in value):
        raise ValueError('expected a list of integers')
    return value


def sieve(bound):
    flags = bytearray(b'\1') * (bound + 1)
    flags[:2] = b'\0\0'
    for p in range(2, isqrt(bound) + 1):
        if flags[p]:
            flags[p*p::p] = b'\0' * ((bound - p*p) // p + 1)
    return flags


def poly_coefficients(poly, length):
    return [int(poly[k]) for k in range(length)]


@lru_cache(maxsize=None)
def exponential_weights(p, length):
    fact = 1
    invfacts = [1]
    for n in range(1, p - 1):
        fact = fact * n % p
        invfacts.append(pow(fact, -1, p))
    return tuple(tuple(pow(j, n, p) * invfacts[n] % p for j in range(length)) for n in range(p-1))


def full_logarithm(coefficients, p):
    """Substitute exp(T), then integrate F'/F, rather than log(1+u)."""
    expansion = [sum(c*w for c,w in zip(coefficients, weights)) % p
                 for weights in exponential_weights(p, len(coefficients))]
    if not expansion[0]:
        raise ValueError('generator is not a unit modulo (1-zeta)')
    series = nmod_poly(expansion, p)
    logarithm = (series.derivative() * series.inverse_series_trunc(p - 2)).truncate(p - 2).integral()
    return poly_coefficients(logarithm, p - 1)


def log_coordinates(coefficients, p):
    return full_logarithm(coefficients, p)[3:p-1:2]


class Check:
    def __init__(self, case_id, root):
        self.case_id, self.directory = case_id, root / case_id
        self.errors, self.missing, self.files, self.data = [], [], {}, {}
        self.start = time.monotonic()

    def expect(self, condition, check, **detail):
        if not condition:
            self.errors.append({'check': check, **detail})
        return condition

    def file(self, name):
        path = self.directory / name
        if not path.is_file():
            self.missing.append(name)
            return None
        self.files[name] = digest(path)
        return path

    def json(self, name):
        path = self.file(name)
        return None if path is None else json.loads(path.read_text())

    def ready(self):
        record = self.json('run.json')
        if record is None:
            return False
        if record.get('status') != 'passed':
            self.missing.append('completed run.json with status passed')
            self.data['reproduction_status'] = record.get('status')
            return False
        return True

    def result(self):
        return {'case_id': self.case_id,
                'status': 'failed' if self.errors else 'incomplete' if self.missing else 'passed_range',
                'elapsed_seconds': round(time.monotonic() - self.start, 6),
                'input_sha256': self.files, 'missing': self.missing,
                'mismatches': self.errors, **self.data}


def check_stdout_hits(check, hits, p, bound, rows, zeros, flags, histogram):
    candidates = [p for p in check.directory.iterdir() if p.is_file() and
                  (p.name.startswith('stdout') or p.name.endswith('.stdout.log'))]
    if not candidates:
        check.missing.append('stdout log containing complete all admitted list')
        return
    found = False
    for path in candidates:
        text = path.read_text().replace('\\\n', '')
        match = re.search(r'all admitted:\s*(\[[\s\S]*?\])', text)
        if match:
            check.file(path.name)
            found = True
            check.expect(integer_list(match.group(1)) == hits, 'complete_stdout_hit_list')
        for prefix in re.finditer(r'^first admitted:\s*(\[[\s\S]*?\])', text, re.M):
            values = integer_list(prefix.group(1))
            check.expect(values == hits[:len(values)], 'stdout_first_admitted_prefix')
        for prefix in re.finditer(r'^admitted:\s*(\d+)\s*(\[[\s\S]*?\])', text, re.M):
            values = integer_list(prefix.group(2))
            check.expect(int(prefix.group(1)) == len(hits) and values == hits[:len(values)], 'stdout_admitted_count_and_prefix')
        for checkpoint in re.finditer(r'^x=(\d+) pi\(x\)=(\d+) admitted=(\d+)', text, re.M):
            x, prime_count, admitted_count = map(int, checkpoint.groups())
            check.expect(x <= bound and prime_count == flags.count(1, 0, x+1)
                         and admitted_count == sum(h <= x for h in hits), 'stdout_prime_count_and_checkpoint', x=x)
        for checkpoint in re.finditer(r'^membership checkpoint:\s*(\d+)\s+(\d+)', text, re.M):
            x, admitted_count = map(int, checkpoint.groups())
            check.expect(x <= bound and admitted_count == sum(h <= x for h in hits), 'stdout_membership_checkpoint', x=x)
        split = re.search(r'^split primes processed:\s*(\d+)', text, re.M)
        if split:
            check.expect(int(split.group(1)) == rows, 'stdout_split_prime_count')
        split = re.search(r'^split (\d+) cells (\d+) count\(0,0,0\)=\s*(\d+)', text, re.M)
        if split:
            check.expect(list(map(int, split.groups())) == [rows, p**3, len(zeros)], 'stdout_split_and_zero_triple_counts')
        if p == 13:
            supplied_histogram = re.search(r'^full tuple histogram by key x\+13\*y\+169\*z:\s*(\[[\s\S]*?\])', text, re.M)
            computed_histogram = [histogram[k] for k in range(13**3)]
            if supplied_histogram:
                check.expect(integer_list(supplied_histogram.group(1)) == computed_histogram, 'within_run_complete_tuple_histogram')
            else:
                check.missing.append('complete tuple histogram in stdout log')
            for label, stride in [('x,0,0', 1), ('0,y,0', 13), ('0,0,w', 169)]:
                match_line = re.search(r'^counts of \('+re.escape(label)+r'\)[^:]*:\s*(\[[\s\S]*?\])', text, re.M)
                if match_line:
                    check.expect(integer_list(match_line.group(1)) == [histogram[stride*k] for k in range(13)], 'within_run_histogram_axis', axis=label)
            distribution = re.search(r'count distribution \(count:frequency\):\s*(\[[\s\S]*?\])', text)
            if distribution:
                supplied = ast.literal_eval(distribution.group(1).replace('<','(').replace('>',')'))
                frequencies = Counter(computed_histogram)
                check.expect(supplied == [(k,frequencies[k]) for k in range(max(computed_histogram)+1)], 'within_run_histogram_count_distribution')
    if not found:
        check.missing.append('complete all admitted list in stdout log')


def verify_prime_case(check, p, bound, flags):
    path = check.file('decisions.csv')
    if path is None:
        return
    expected = (ell for ell in range(p + 1, bound + 1, p) if flags[ell])
    phi = fmpz_poly([1] * p)
    hits, zeros, rows, seen, bad_norms, bad_logs = [], [], 0, 0, 0, 0
    histogram = Counter()
    coordinates_sha = hashlib.sha256()
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter=';')
        check.expect(reader.fieldnames == ['ell', 'norm', 'power_basis_coefficients', 'higher_odd_coordinates', 'admitted'], 'decision_header', actual=reader.fieldnames)
        for rows, row in enumerate(reader, 1):
            expected_ell = next(expected, None)
            try:
                ell, claimed_norm, admitted = int(row['ell']), int(row['norm']), int(row['admitted'])
                coefficients, supplied = integer_list(row['power_basis_coefficients']), integer_list(row['higher_odd_coordinates'])
                check.expect(ell == expected_ell, 'ordered_complete_prime_coverage', row=rows, actual=ell, expected=expected_ell)
                check.expect(2 <= ell <= bound and bool(flags[ell]) and ell % p == 1, 'prime_and_split', row=rows, ell=ell)
                check.expect(len(coefficients) == p - 1, 'power_basis_length', row=rows, actual=len(coefficients))
                if len(coefficients) != p - 1:
                    continue
                norm = int(phi.resultant(fmpz_poly(coefficients)))
                if not check.expect(norm == ell == claimed_norm and norm > 0, 'exact_integer_norm', row=rows, ell=ell, computed=norm, supplied=claimed_norm):
                    bad_norms += 1
                coordinates = log_coordinates(coefficients, p)
                if p == 13:
                    histogram[coordinates[0] + 13*coordinates[1] + 169*coordinates[2]] += 1
                if not check.expect(supplied == coordinates, 'higher_odd_log_coordinates', row=rows, ell=ell, computed=coordinates, supplied=supplied):
                    bad_logs += 1
                criterion = all(c == 0 for c in coordinates[:-2]) and any(coordinates[-2:])
                check.expect(admitted in (0, 1) and bool(admitted) == criterion, 'admission_iff_criterion', row=rows, ell=ell, computed=int(criterion), supplied=admitted)
                if criterion:
                    hits.append(ell)
                    check.expect(ell * 2 ** (p - 1) > p ** (p - 4), 'strict_prime_size_bound', ell=ell)
                    if p == 13:
                        check.expect(coordinates[-2] != 0, 'p13_exceptional_branch_in_range', ell=ell)
                if all(c == 0 for c in coordinates[:-2]):
                    zeros.append(ell)
                coordinates_sha.update(f'{ell};{coordinates};{int(criterion)}\n'.encode())
                seen += 1
            except (ValueError, KeyError, TypeError, SyntaxError) as exc:
                check.errors.append({'check': 'parse_or_arithmetic', 'row': rows, 'error': str(exc)})
    remainder = list(expected)
    check.expect(not remainder, 'no_missing_terminal_primes', missing_count=len(remainder), first_missing=remainder[:10])
    checkpoints = {100000: 2, 1000000: 54, 10000000: 537, 30000000: 1510} if p == 11 else {10000000: 3, 12000000: 4, 40000000: 31}
    counts = {str(x): sum(h <= x for h in hits) for x in checkpoints if x <= bound}
    for x, want in checkpoints.items():
        if x <= bound:
            check.expect(counts[str(x)] == want, 'paper_checkpoint_count', x=x, expected=want, actual=counts[str(x)])
    check.expect(bool(hits) and hits[0] == (77419 if p == 11 else 6285241), 'least_admitted_prime', actual=hits[:1])
    if p == 13:
        check.expect(len(zeros) == (4 if bound == 12000000 else 31), 'p13_zero_triple_count', actual=len(zeros))
        check.expect(zeros == hits, 'p13_zero_triples_equal_admitted_set')
        check.expect(rows == (65618 if bound == 12000000 else 202714), 'p13_split_count', actual=rows)
    check_stdout_hits(check, hits, p, bound, rows, zeros, flags, histogram)
    check.data.update(claim_ids=['PC1'], prime=p, inclusive_bound=bound, rows=rows,
                      exact_norms_tested=seen, independently_computed_logarithms=seen,
                      norm_mismatches=bad_norms, logarithm_mismatches=bad_logs,
                      checkpoint_counts=counts, admitted_primes=hits, zero_low_coordinate_primes=zeros,
                      canonical_decisions_sha256=coordinates_sha.hexdigest(),
                      choice_dependent_histogram_checked_within_same_run=(p == 13),
                      method='Independent Eratosthenes sieve; integer resultant; formal derivative/inverse/integration; exact finite-field admission criterion.',
                      scope_limit='Finite inclusive range only. Neither convergence rate nor asymptotic density follows from these counts. Nonzero tuple frequencies are checked only against the generators in the same run; they are not asserted independent of prime-ideal choices.')


def irregular_indices_exact(p, bernoulli_numerators):
    return [n for n in range(2, p - 1, 2) if bernoulli_numerators[n] % p == 0]


def least_r(p, indices):
    halves = [n // 2 for n in indices]
    return next(r for r in range(1, p) if all((2*r + 1)*i > r*(p - 1) for i in halves if i % r == 0))


def verify_census(check):
    data, path = check.json('census.json'), check.file('census_all_primes.csv')
    if data is None or path is None:
        return
    flags = sieve(19999)
    primes = [p for p in range(7, 20000) if flags[p]]
    # Exact rational Bernoulli numbers are independent of modular series inversion.
    fmpq.bernoulli(19996, cache=True)
    bernoulli = {n: fmpq.bernoulli(n).numerator for n in range(2, 19997, 2)}
    histogram, irregular, computed_rows = Counter(), [], []
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter=';'))
    check.expect(len(rows) == len(primes) == 2259, 'census_complete_row_count', actual=len(rows))
    for index, p in enumerate(primes):
        indices = irregular_indices_exact(p, bernoulli)
        r = least_r(p, indices)
        histogram[r] += 1
        item = {'p': p, 'irregular_indices': indices, 'r_p': r, 'tail_exponent': p + 3 + 2*r}
        computed_rows.append(item)
        if indices:
            irregular.append(item)
        check.expect(r <= (p - 5)//2, 'census_uniform_r_bound', p=p, r=r)
        if index >= len(rows):
            continue
        row = rows[index]
        actual = {'p': int(row['p']), 'irregular_indices': integer_list(row['irregular_indices']), 'r_p': int(row['r_p']), 'tail_exponent': int(row['p_plus_3_plus_2r_p'])}
        check.expect(actual == item, 'exact_bernoulli_census_row', row=index+1, expected=item, actual=actual)
        check.expect(str(row['within_cited_vandiver_range']).lower() in ('true', '1'), 'census_reported_vandiver_range_flag', p=p)
    expected_hist = {1: 1602, 2: 283, 3: 224, 4: 67, 5: 59, 6: 6, 7: 16, 8: 2}
    check.expect(dict(histogram) == expected_hist, 'paper_r_histogram', actual=dict(histogram))
    check.expect(len(irregular) == 900, 'paper_irregular_count', actual=len(irregular))
    check.expect(data['irregular_primes'] == irregular, 'census_json_all_irregular_rows')
    summary = data['summary']
    for key, expected in [('PMAX', 20000), ('primes', 2259), ('irregular', 900), ('r_p_histogram', {str(k):v for k,v in expected_hist.items()}), ('max_r_p', 8), ('smallest_with_r_p_gt_1', 103), ('max_index_of_irregularity', 4)]:
        check.expect(summary.get(key) == expected, 'census_summary', field=key, expected=expected, actual=summary.get(key))
    check.expect(summary.get('max_r_p_le_(p-5)/2_for_all') is True, 'census_reported_uniform_bound')
    check.expect(summary.get('fraction_r_p_gt_1') == 657/2259, 'census_fraction_above_first_tail')
    archival = SUPPLEMENT / 'outputs/python/tail_exponent_census_2e4.json'
    if archival.is_file():
        check.files['archived:outputs/python/tail_exponent_census_2e4.json'] = digest(archival)
        check.expect(json.loads(archival.read_text())['irregular_primes'] == irregular, 'archived_irregular_rows_against_exact_bernoulli')
    check.data.update(claim_ids=['CEN1'], prime_interval='7 <= p < 20000', rows=len(primes), irregular_primes=len(irregular),
                      r_p_histogram=dict(histogram), complete_rows_sha256=canonical_digest(computed_rows),
                      exact_Bernoulli_even_indices='2,4,...,19996',
                      method='Exact rational FLINT Bernoulli numerators, separately reduced for every prime; direct admissibility inequality; independent sieve.',
                      scope_limit='All finite census rows checked. No independent verification of the cited external Vandiver theorem/range is claimed; recorded flags are checked as metadata only.')


def unit_ray(data, p, e, augmentation, primitive_root, check, label):
    q = p*p
    phi = nmod_poly([1] * p, q)
    exponents = {int(a): int(n) for a, n in data['circular_unit_exponents'].items()}
    final = {int(a): int(n) for a, n in data['final_circular_unit_exponents'].items()}
    correction = int(data['constant_correction_exponent_of_u_g'])
    expected_final = dict(exponents)
    expected_final[primitive_root] = expected_final.get(primitive_root, 0) + correction
    check.expect(final == expected_final and correction == p * data['constant_correction_j'], 'unit_constant_correction_exponents', ray=label)
    unit = nmod_poly([1], q)
    constant = 1
    for a, n in final.items():
        if not (1 <= a < p and n >= 0):
            raise ValueError(f'invalid circular unit exponent {a}:{n}')
        unit = (unit * nmod_poly([1] * a, q).pow_mod(n, phi)) % phi
        constant = constant * pow(a, n, p) % p
    check.expect(constant == data['target_constant_mod_p'] % p, 'unit_constant_residue', ray=label, actual=constant)
    zeta_power = nmod_poly([0, 1], q).pow_mod(int(data['zeta_exponent']) % p, phi)
    value = (unit * zeta_power * nmod_poly([1, -1], q).pow_mod(e, phi)) % phi
    coefficients = poly_coefficients(value, p)
    check.expect((augmentation - sum(coefficients)) % p == 0, 'ray_integral_augmentation_correction', ray=label)
    c = (augmentation - sum(coefficients)) // p
    return [(a + c) % p for a in coefficients]


def verify_tail_witness(check, witness, attempt, p, r):
    e, s, g = 2*r + 1, p*p, witness['primitive_root']
    for field, want in [('p', p), ('r', r), ('e', e), ('augmentation', s), ('arithmetic_modulus', p*p)]:
        check.expect(witness.get(field) == want, 'witness_parameter', r=r, field=field, expected=want, actual=witness.get(field))
    check.expect(len({pow(g, k, p) for k in range(p - 1)}) == p - 1, 'primitive_root', p=p, g=g)
    rays = witness['rays_mod_p']
    labels = [str(b) for b in range(p)] + ['inf']
    check.expect(set(rays) == set(labels), 'all_ray_labels', r=r)
    check.expect(canonical_digest(rays) == witness['rays_canonical_json_sha256'], 'canonical_ray_digest', r=r)
    for label in labels:
        ray = rays[label]
        if not check.expect(len(ray) == p and all(type(x) is int and 0 <= x < p for x in ray), 'canonical_ray_shape', r=r, ray=label):
            raise ValueError('malformed ray')
        check.expect(sum(ray) % p == 0, 'ray_common_augmentation_mod_p', r=r, ray=label)
        unit_data = witness['exceptional_unit_data'] if label == 'inf' else witness['generic_unit_data'][label]
        actual = unit_ray(unit_data, p, e if label == 'inf' else 1, s, g, check, label)
        check.expect(actual == ray, 'unit_exponents_reproduce_complete_ray_mod_p', r=r, ray=label, **({'computed':actual,'supplied':ray} if actual != ray else {}))
    failures = [[i, j] for i in range(p) for j in range(p)
                if (rays['inf'][j] + sum(rays[str(b)][(i + b*j) % p] for b in range(p)) - s) % p]
    check.expect(failures == witness['failing_cells'], 'all_reconstruction_cells', r=r)
    fact, invfacts = 1, [1]
    for n in range(1, p):
        fact = fact * n % p
        invfacts.append(pow(fact, -1, p))
    power_matrix = nmod_mat([[pow(t, n, p) * invfacts[n] % p for n in range(p)] for t in range(p)], p)
    transforms = nmod_mat([rays[label] for label in labels], p) * power_matrix
    generic = [[int(transforms[b, n]) for n in range(p)] for b in range(p)]
    infinity = [int(transforms[p, n]) for n in range(p)]
    tau = (s // p) % p
    for row in generic:
        row[p-1] = (row[p-1] - tau) % p
    infinity[p-1] = (infinity[p-1] - tau) % p
    check.expect(generic == witness['generic_transform_coefficients'], 'all_generic_transforms_from_rays', r=r)
    check.expect(infinity == witness['exceptional_transform_coefficients'], 'exceptional_transform_from_ray', r=r)
    # Lagrange basis: delta_b(X)=1-(X-b)^(p-1). This produces the
    # unique degree <= p-1 representative, including b=0 correctly.
    interpolation = nmod_mat([[int(b == 0)] + [(-pow(b, p-1-k, p)) % p for k in range(1, p)] for b in range(p)], p)
    polys = nmod_mat(generic, p).transpose() * interpolation
    degree_failures, infinity_failures = [], []
    for n in range(1, p - 1):
        coefficients = [int(polys[n, k]) for k in range(p)]
        degree = max((k for k, c in enumerate(coefficients) if c), default=-1)
        if degree > n:
            degree_failures.append([n, degree])
        if coefficients[n] != infinity[n]:
            infinity_failures.append(n)
        check.expect(coefficients == witness['moment_polynomials_mod_p'][str(n)], 'canonical_moment_polynomial', r=r, n=n)
    constant = (sum(row[p-1] for row in generic) + infinity[p-1] + tau) % p
    check.expect(constant == 0, 'constant_moment', r=r, actual=constant)
    ij_pattern = failures == [[i,j] for i in range(1,p) for j in range(1,p)]
    computed = {'congruence_failures': len(failures), 'cells': p*p,
                'failures_exactly_at_cells_with_ij_nonzero': ij_pattern,
                'failing_degree_conditions': degree_failures,
                'failing_infinity_conditions': infinity_failures,
                'determinant_exponent': p + 3 + 2*r}
    for key, value in computed.items():
        check.expect(attempt.get(key) == value, 'result_matches_independent_reconstruction', r=r, field=key, computed=value, supplied=attempt.get(key))
        check.expect(witness['result'].get(key) == value, 'witness_result_matches_independent_reconstruction', r=r, field=key)
    expected_bad = (p == 103 and r < 5) or (p in (131,157) and r == 1)
    expected_degrees = {(103,1):[[25,36]], (103,2):[[25,30]], (103,3):[[25,28]], (103,4):[[25,27]], (131,1):[[23,33]], (157,1):[[63,93]]}
    check.expect(len(failures) == ((p-1)**2 if expected_bad else 0), 'paper_failure_count', r=r)
    check.expect(ij_pattern == expected_bad, 'paper_failure_pattern', r=r)
    check.expect(degree_failures == expected_degrees.get((p,r), []) and not infinity_failures, 'paper_moment_failures', r=r, degrees=degree_failures, infinity=infinity_failures)
    # Compare the entire prescribed H_b=T P(a b^e T^(2r)), not just tests.
    P = {int(k):int(v) for k,v in witness['P_coefficients'].items()}
    a = witness.get('construction_a', 1)
    for b in range(p):
        prescribed = [0]*p
        prescribed[1] = 1
        for k, c in P.items():
            if 2*r*k + 1 < p:
                prescribed[2*r*k+1] = c * pow(a * pow(b,e,p), k, p) % p
        check.expect(generic[b] == prescribed, 'prescribed_generic_tail_expansion', r=r, b=b)
    prescribed_inf = [0]*p
    prescribed_inf[e] = a % p
    check.expect(infinity == prescribed_inf, 'prescribed_exceptional_tail_expansion', r=r)
    return {'r':r, 'e':e, 'cells_checked':p*p, 'unit_rays_checked':p+1,
            'moment_polynomials_checked':p-2, **computed}


def verify_tail_case(check, p, rmax):
    data = check.json('result.json')
    if data is None:
        return
    exact_irregular = [n//2 for n in range(2, p-1, 2) if fmpq.bernoulli(n).numerator % p == 0]
    check.expect(data['p'] == p and data['irregular_half_indices'] == exact_irregular,
                 'exact_bernoulli_irregular_half_indices', expected=exact_irregular)
    admissible = [r for r in range(1, 11) if all((2*r+1)*i > r*(p-1) for i in exact_irregular if i % r == 0)]
    check.expect(data['admissible_r_le_10'] == admissible, 'all_reported_admissibility_flags')
    attempts = data['attempts']
    check.expect([a['r'] for a in attempts] == list(range(1, rmax+1)), 'complete_attempt_sequence')
    summaries = []
    for r in range(1, rmax+1):
        witness = check.json(f'witnesses/p{p}_r{r}.json')
        if witness is None or r > len(attempts):
            continue
        attempt = attempts[r-1]
        check.expect(attempt.get('witness_sha256') == check.files[f'witnesses/p{p}_r{r}.json'], 'witness_file_digest', r=r)
        check.expect(attempt.get('admissible') == (r in admissible), 'attempt_admissibility', r=r)
        check.expect(witness['P_coefficients'] == attempt['P_coefficients'], 'P_coefficients_result_witness_agree', r=r)
        summaries.append(verify_tail_witness(check, witness, attempt, p, r))
    # Direct logarithmic derivative of every circular unit u_a in the chosen
    # generating set, compared with exact rational Bernoulli coefficients.
    units_checked = 0
    if summaries:
        invfacts, factorial = [1], 1
        for n in range(1, p-1):
            factorial = factorial*n % p
            invfacts.append(pow(factorial, -1, p))
        for a in witness['circular_unit_indices']:
            actual = full_logarithm([1]*a, p)
            expected = [0]*(p-1)
            expected[1] = (a-1)*pow(2,-1,p) % p
            for n in range(2, p-1, 2):
                bn = fmpq.bernoulli(n)
                expected[n] = (int(bn.numerator % p)*pow(int(bn.denominator % p),-1,p)
                               *pow(n,-1,p)*invfacts[n]*(pow(a,n,p)-1)) % p
            check.expect(actual == expected, 'exact_bernoulli_formula_vs_unit_log_derivative', a=a)
            units_checked += 1
        check.expect(data.get('bernoulli_formula_matches_direct_unit_logs') is True,
                     'producer_bernoulli_formula_flag')
    check.data.update(claim_ids=['OBS1','F2','F3','F4'], prime=p, attempts=summaries,
                      exact_bernoulli_half_indices=exact_irregular, circular_unit_logarithms_checked=units_checked,
                      method='Independent cyclotomic quotient arithmetic modulo p² from all recorded circular-unit exponents; every p² cell; moment transforms from rays; Lagrange interpolation in F_p[b]/(b^p-b).',
                      scope_limit='Finite residue certificates only; no huge full integer determinant vectors or general theorem proof are asserted by these reconstruction cases.')


def verify_vector(check):
    vector_path = check.file('determinant_example_p11_vector.txt')
    matrix_path = check.file('determinant_example_p11_matrix_entries.csv')
    if vector_path is None or matrix_path is None:
        return
    p, ell = 11, 77419
    vector = ast.literal_eval(vector_path.read_text().replace('\\\n', '').replace('\\\r\n', ''))
    if not check.expect(len(vector)==p and all(len(row)==p and all(type(v) is int for v in row) for row in vector), 'vector_shape'):
        return
    augmentation = sum(map(sum, vector))
    check.expect(augmentation == p, 'vector_augmentation', actual=augmentation)
    matrix = [[vector[(x//p-w//p)%p][(x%p-w%p)%p] for x in range(p*p)] for w in range(p*p)]
    seen = set()
    with matrix_path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter=';')
        check.expect(reader.fieldnames == ['row_index','column_index','value'], 'matrix_header', actual=reader.fieldnames)
        for record in reader:
            i,j,value = int(record['row_index']),int(record['column_index']),int(record['value'])
            # The public CSV deliberately uses zero-based row/column indices.
            check.expect((i,j) not in seen and 0 <= i < p*p and 0 <= j < p*p, 'matrix_index_unique_and_valid', row=i, column=j)
            seen.add((i,j))
            if 0 <= i < p*p and 0 <= j < p*p:
                check.expect(value == matrix[i][j], 'every_group_matrix_entry', row=i, column=j, supplied=value, expected=matrix[i][j])
    check.expect(len(seen)==p**4, 'complete_group_matrix_entries', actual=len(seen))
    direct = int(fmpz_mat(matrix).det())
    expected = p**15 * ell
    check.expect(direct == expected, 'direct_exact_121_by_121_determinant', actual=direct, expected=expected)
    phi = fmpz_poly([1]*p)
    norms = []
    for b in range(p):
        ray = [sum(vector[i][j] for i in range(p) for j in range(p) if (i+b*j)%p==t) for t in range(p)]
        norm = int(phi.resultant(fmpz_poly(ray)))
        norms.append(norm)
        check.expect(norm == p, 'generic_projection_norm', b=b, actual=norm)
    exceptional = [sum(vector[i][j] for i in range(p)) for j in range(p)]
    exceptional_norm = int(phi.resultant(fmpz_poly(exceptional)))
    check.expect(exceptional_norm == p**3*ell, 'exceptional_projection_norm', actual=exceptional_norm)
    product = augmentation * exceptional_norm
    for norm in norms:
        product *= norm
    check.expect(product == direct, 'norm_factorization_matches_direct_matrix_determinant')
    check.expect(max(abs(x) for row in vector for x in row) < 10**36, 'paper_coefficient_size')
    check.data.update(claim_ids=['DET1','F2'], matrix_dimension=121, matrix_entries_checked=len(seen), augmentation=augmentation,
                      direct_exact_determinant=direct, expected_determinant=expected, generic_norms=norms,
                      exceptional_norm=exceptional_norm, maximum_absolute_coefficient=max(abs(x) for row in vector for x in row),
                      method='Every entry independently reconstructed from v[x-w], then FLINT exact integer matrix determinant; separate exact character-norm factorization.')


def verify_generator(check):
    data = check.json('certificate.json')
    if data is None:
        return
    certificate = data['certificate']
    coeff = [-2,-2,-3,-2,-1,-2,-1,-3,0,0]
    norm = int(fmpz_poly([1]*11).resultant(fmpz_poly(coeff)))
    coordinates = log_coordinates(coeff,11)
    check.expect(certificate['power_basis_coefficients']==coeff, 'compact_generator_coefficients')
    check.expect(certificate['norm']==norm==77419, 'compact_generator_exact_norm')
    check.expect(certificate['odd_logarithmic_coordinates']=={str(k):v for k,v in zip([3,5,7,9],coordinates)} and coordinates==[0,0,10,9], 'compact_generator_log_coordinates')
    check.expect(certificate['constant_residue']==sum(coeff)%11==6, 'compact_generator_constant')
    check.expect([c*pow(3,j,11)%11 for c,j in zip([0,0,6,5],[3,5,7,9])]==coordinates, 'historical_coordinate_galois_multiplier')
    check.data.update(claim_ids=['GEN1'], norm=norm, higher_odd_coordinates=coordinates,
                      method='Exact integer resultant and independent formal derivative/inverse/integration.')


def main():
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root',type=Path,default=SUPPLEMENT/'reproduced_outputs')
    parser.add_argument('--output',type=Path,default=SUPPLEMENT/'independent_verification.json')
    parser.add_argument('--case',action='append',choices=CASE_IDS,help='Verify only this case (repeatable); omitted means all cases.')
    args=parser.parse_args()
    started=time.monotonic()
    selected=args.case or CASE_IDS
    report={'schema_version':1,'generated_at_utc':datetime.now(timezone.utc).isoformat(),
            'verifier_sha256':digest(Path(__file__)), 'input_root':str(args.root.resolve()),
            'command_argv':[sys.executable]+sys.argv,
            'environment':{'python':platform.python_version(),'python_implementation':platform.python_implementation(),
                           'platform':platform.platform(),'python_flint':getattr(flint,'__version__',None)},
            'independence':'Second implementation from manuscript definitions. No author program imported or executed; author outputs are claims to check, not correctness evidence.',
            'source_locators':SOURCE_LOCATORS,'selected_cases':selected,'cases':[],
            'limits':['Only listed finite cases and ranges are checked. No convergence rate, infinite assertion, proof correctness, novelty, or author approval is certified.',
                      'The verifier shares the FLINT exact arithmetic library with some producer computations; independent algorithms and representations are used as stated per case.']}
    manifest=SUPPLEMENT/'inputs/run_cases.json'
    if manifest.is_file():
        report['case_manifest_sha256']=digest(manifest)
    manuscript=SUPPLEMENT.parent/'manuscript/main.tex'
    if manuscript.is_file():
        report['manuscript_snapshot_sha256']=digest(manuscript)
    if flint is None:
        report.update(status='incomplete',missing=['python-flint is required'])
    else:
        flags=None
        for case_id in selected:
            check=Check(case_id,args.root)
            try:
                if check.ready():
                    print(f'Verifying {case_id}',flush=True)
                    if case_id.startswith('magma_p'):
                        p,b=map(int,re.fullmatch(r'magma_p(\d+)_B(\d+)',case_id).groups())
                        if flags is None:
                            flags=sieve(40000000)
                        verify_prime_case(check,p,b,flags)
                    elif case_id=='python_census_20000':
                        verify_census(check)
                    elif case_id.startswith('python_tail_'):
                        p,r=map(int,re.fullmatch(r'python_tail_p(\d+)_r(\d+)',case_id).groups())
                        verify_tail_case(check,p,r)
                    elif case_id=='magma_vector_p11_77419':
                        verify_vector(check)
                    else:
                        verify_generator(check)
            except Exception as exc:
                check.errors.append({'check':'verifier_exception','type':type(exc).__name__,'error':str(exc)})
            result=check.result()
            report['cases'].append(result)
            print(f"{case_id}: {result['status']} ({len(check.errors)} mismatches, {len(check.missing)} missing)",flush=True)
        statuses=[c['status'] for c in report['cases']]
        report['status']='failed' if 'failed' in statuses else 'incomplete' if 'incomplete' in statuses else 'passed_range'
        # The overlap comparison is invariant under choosing another prime ideal.
        p13={c['case_id']:c for c in report['cases'] if c.get('prime')==13}
        short=p13.get('magma_p13_B12000000'); long=p13.get('magma_p13_B40000000')
        if short and long and 'admitted_primes' in short and 'admitted_primes' in long:
            good=short['admitted_primes']==[x for x in long['admitted_primes'] if x<=12000000]
            report['p13_overlap_admission_sets_agree']=good
            if not good:
                report['status']='failed'
    report['elapsed_seconds']=round(time.monotonic()-started,6)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(f"{report['status']}: {args.output}",flush=True)
    return {'passed_range':0,'failed':1,'incomplete':2}[report['status']]


if __name__=='__main__':
    sys.exit(main())
