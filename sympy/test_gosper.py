from sympy import Symbol, binomial
from gosper import compute_bc, compute_pqr


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


def test_compute_pqr_poly():
    k = Symbol('k', integer=True)
    b = k * (k + 3)
    c = k + 1
    p, q, r = compute_pqr(b, c)
    assert p.equals()