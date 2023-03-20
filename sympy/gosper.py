from sympy import fraction, combsimp, cancel, resultant, roots, Symbol, Dummy
from sympy import gcd, degree


def compute_bc(a, k):
    x = a.subs(k, k + 1) / a
    x = combsimp(x)
    x = cancel(x)
    b, c = fraction(x)
    # TODO: ensure b, c are polys in k
    return b, c


def compute_pqr(b, c, k):
    # TODO: precons
    j = Dummy("j")
    rp = resultant(b, c.subs(k, k+j))
    rz = roots(rp, multiple=True, filter="Z")

    p = 1
    q = b
    r = c

    for z in rz:
        if z.is_negative: continue
        g = gcd(q, r.subs(k, k + z))
        q /= g
        r /= g.subs(k, k - z)
        for i in range(0, -z, -1):
            p *= g.subs(k, k + i)

    # TODO: postconditions
    return p, q, r


def degree_bound_f(p, q, r, k):
    x = q + r.subs(k, k-1)
    y = q - r.subs(k, k-1)

    n = degree(x, k)
    print(x, y)
    if n <= degree(y, k):
        return degree(p.subs(k, k-1), k) - degree(y, k)
    
    a = x.coeff(k, n)
    b = y.coeff(k, n - 1)
    assert not a.equals(0)
    print(a, b, n)

    z = (-2 * b) / a
    if not z.is_integer or z < 0:
        return degree(p.subs(k, k-1), k) - n + 1
    
    print(z, degree(p.subs(k, k-1), k) - n + 1)
    return max(z, degree(p.subs(k, k-1), k) - n + 1)


def gosper_sum(a, k):
    b, c = compute_bc(a, k)
    p, q, r = compute_pqr(b, c, k)
    d = degree_bound_f(q, p, r, k)
