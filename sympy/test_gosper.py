from sympy import Symbol, binomial, Integer
from gosper import compute_bc, compute_pqr, degree_bound_f


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
    n = Symbol('n', integer=True)
    a = ((-1) ** k) * binomial(n, k)
    b, c = compute_bc(a, k)
    assert b.equals(k - n)
    assert c.equals(k + 1)


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
    q = k * (k + 1)
    r = 2 * (k + 2) * (k + 3)
    d = degree_bound_f(p, q, r, k)
    assert d == 3


def test_degree_bound_f_case_2_z_neg():
    pass


def test_degree_bound_f_case_3():
    pass