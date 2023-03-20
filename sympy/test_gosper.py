from sympy import Symbol, binomial, Integer, Rational, factorial
from gosper import compute_bc, compute_pqr, degree_bound_f, compute_f, gosper_sum
import pytest


def test_compute_bc_poly():
    k = Symbol('k', integer=True)
    a = k ** 2
    b, c = compute_bc(a, k)
    assert b.equals(k ** 2 + 2 * k + 1)
    assert c.equals(k ** 2)


def test_compute_bc_rat():
    k = Symbol('k', integer=True)
    a = (3 * (k ** 2)) / k
    b, c = compute_bc(a, k)
    assert b.equals(k + 1)
    assert c.equals(k)


def test_compute_bc_exp():
    k = Symbol('k', integer=True)
    a = k * (2 ** k)
    b, c = compute_bc(a, k)
    assert b.equals(2 * k + 2)
    assert c.equals(k)


def test_compute_bc_binom():
    k = Symbol('k', integer=True)
    n = Symbol('n')
    a = ((-1) ** k) * binomial(n, k)
    b, c = compute_bc(a, k)
    assert b.equals(k - n)
    assert c.equals(k + 1)


def test_compute_bc_fact():
    pass


def test_compute_pqr():
    k = Symbol('k', integer=True)
    b = k * (k + 3)
    c = k + 1
    p, q, r = compute_pqr(b, c, k)
    assert p.equals((k + 2) * (k + 3))
    assert q.equals(k)
    assert r.equals(1)


def test_degree_bound_f_case_1():
    k = Symbol('k', integer=True)
    p = (k + 2) * (k + 3)
    q = k
    r = Integer(1)
    d = degree_bound_f(p, q, r, k)
    assert d == 1


def test_degree_bound_f_case_2_z_0():
    k = Symbol('k', integer=True)
    p = k ** 4
    q = k
    r = k + 1
    d = degree_bound_f(p, q, r, k)
    assert d == 4


def test_degree_bound_f_case_2_z_float():
    k = Symbol('k', integer=True)
    p = k ** 4
    q = k ** 2 + 3 * k
    r = k ** 2
    d = degree_bound_f(p, q, r, k)
    assert d == 3


def test_degree_bound_f_case_2_z_neg():
    k = Symbol('k', integer=True)
    p = k ** 4
    q = k ** 2 + 2 * k
    r = k ** 2
    d = degree_bound_f(p, q, r, k)
    assert d == 3


def test_degree_bound_f_case_3_z_less():
    k = Symbol('k', integer=True)
    p = k ** 4
    q = k ** 2 - k
    r = k ** 2
    d = degree_bound_f(p, q, r, k)
    assert d == 3


def test_degree_bound_f_case3_z_more():
    k = Symbol('k', integer=True)
    p = k
    q = k ** 2 - Rational(1, 2) * k
    r = k ** 2
    d = degree_bound_f(p, q, r, k)
    assert d == 1


def test_compute_f():
    pass


def gosper_sum_data():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)
    return [
        (1, k, k + 1),
        (k, k, (k * (k + 1)) / 2),
        (k ** 2, k, (k * (k + 1) * (2 * k + 1)) / 6),
        (k ** 3, k, ((k * (k + 1)) / 2) ** 2),
        (2 ** k, k, 2 ** (k + 1)),
        (k * (2 ** k), k, (k - 1) * (2 ** (k + 1))),
        ((k ** 2) * (2 ** k), k, 2 * (k ** 2 - 2 * k + 3) * (2 ** k)),
        ((1 / (k + 1)) - (1 / k), k, (1 / (k + 1))),
        (k * factorial(k), k, factorial(k + 1)),
        (((-1) ** k) * binomial(n, k), k, -(k * (-1) ** k * binomial(n, k)) / n),
    ]


@pytest.mark.parametrize("a,k,s_exp", gosper_sum_data())
def test_gosper_sum(a, k, s_exp):
    s_act = gosper_sum(a, k)
    assert s_exp.equals(s_act)