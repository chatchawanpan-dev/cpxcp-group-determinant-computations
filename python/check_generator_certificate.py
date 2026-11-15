#!/usr/bin/env python3
"""Exact certificate for the p=11 generator printed in Section 8.

Uses Python's standard library only. The norm is a Bareiss determinant of
multiplication in Z[X]/(1+X+...+X^10). Logarithmic coordinates are calculated
independently by substituting X=exp(T) modulo (11,T^10).
The certificate concerns this generator, not the historical 121x121 vector.
"""
import json
import platform
import sys


def bareiss(matrix):
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a) - 1):
        if a[k][k] == 0:
            pivot_row = next((i for i in range(k + 1, len(a)) if a[i][k]), None)
            if pivot_row is None:
                return 0
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        for i in range(k + 1, len(a)):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def certificate():
    p = 11
    coefficients = [-2, -2, -3, -2, -1, -2, -1, -3, 0, 0]
    n = p - 1
    columns = []
    for j in range(n):
        product = [0] * (2 * n - 1)
        for i, value in enumerate(coefficients):
            product[i + j] += value
        for degree in range(len(product) - 1, n - 1, -1):
            value = product[degree]
            for k in range(n):
                product[degree - n + k] -= value
            product[degree] = 0
        columns.append(product[:n])
    matrix = [[columns[j][i] for j in range(n)] for i in range(n)]
    norm = bareiss(matrix)
    factorial = [1] * n
    for j in range(1, n):
        factorial[j] = factorial[j - 1] * j % p
    series = [sum(a * pow(i, j, p) for i, a in enumerate(coefficients))
              * pow(factorial[j], -1, p) % p for j in range(n)]
    constant = series[0]
    u = [a * pow(constant, -1, p) % p for a in series]
    u[0] = 0

    def multiply(a, b):
        return [sum(a[i] * b[j - i] for i in range(j + 1)) % p
                for j in range(n)]

    logarithm, power = [0] * n, u
    for k in range(1, n):
        factor = (1 if k % 2 else -1) * pow(k, -1, p)
        logarithm = [(a + factor * b) % p for a, b in zip(logarithm, power)]
        power = multiply(power, u)
    odd = {str(j): logarithm[j] for j in range(3, p - 1, 2)}
    assert norm == 77419
    assert list(odd.values()) == [0, 0, 10, 9]
    # These are conjugate to the historical author's coordinates (0,0,6,5).
    assert [(pow(3, j, p) * a) % p
            for j, a in zip(range(3, p - 1, 2), [0, 0, 6, 5])] == list(odd.values())
    return {"p": p, "power_basis_coefficients": coefficients,
            "norm": norm, "constant_residue": constant,
            "odd_logarithmic_coordinates": odd,
            "historical_example_galois_multiplier": 3,
            "checks": "exact integer norm and finite-field logarithm; no full vector rebuild"}


if __name__ == "__main__":
    result = {"python": platform.python_version(), "certificate": certificate()}
    output = json.dumps(result, indent=2) + "\n"
    if len(sys.argv) == 2:
        with open(sys.argv[1], "w") as stream:
            stream.write(output)
    print(output, end="")
