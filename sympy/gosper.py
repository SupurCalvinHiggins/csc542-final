from sympy import fraction, combsimp, cancel, resultant, roots, Symbol, Dummy
from sympy import gcd, degree, linsolve, simplify


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
    p = cancel(p)
    q = cancel(q)
    r = cancel(r)
    return p, q, r


def degree_bound_f(p, q, r, k):
    x = q + r.subs(k, k-1)
    y = q - r.subs(k, k-1)

    n = degree(x, k)
    if n <= degree(y, k):
        return degree(p.subs(k, k-1), k) - degree(y, k)
    
    a = x.coeff(k, n)
    b = y.coeff(k, n - 1)
    assert not a.equals(0)

    z = (-2 * b) / a
    if not z.is_integer or z < 0:
        return degree(p.subs(k, k-1), k) - n + 1
    
    return max(z, degree(p.subs(k, k-1), k) - n + 1)


def compute_f(p, q, r, d, k):
    f = 0
    cs = []
    for i in range(0, d + 1):
        c = Dummy(f"c{i}")
        f += c * (k ** i)
        cs.append(c)
    
    print(f)
    p.subs(k, k-1)
    e = q * f - r.subs(k, k-1) * f.subs(k, k-1) - p.subs(k, k-1)
    print("e: ", e)
    print("e: ", simplify(e))
    es = []
    for i in range(0, 2 * (d + 1)): # TODO: idk why we need * 2
        es.append(e.subs(k, i))
    
    print(es)
    ss = linsolve(es, cs)
    if len(ss) != 1:
        raise ValueError()

    s = next(iter(ss))
    f = f.subs(list(zip(cs, s)))
    f = f.subs(list(zip(cs, [0] * len(cs))))
    # print(f)
    return f


def compute_s(a, p, r, f, k):
    s = (r.subs(k, k-1) / p.subs(k, k-1)) * f.subs(k, k-1) * a
    s = cancel(s)
    return s

k = Symbol('k', integer=True)
# p = k
# q = k ** 2 - k
# r = k ** 2
# d = degree_bound_f(p, q, r, k)
# print(d)
# compute_f(p, q, r, d, k)


def gosper_sum(a, k):
    b, c = compute_bc(a, k)
    p, q, r = compute_pqr(b, c, k)
    print(p, q, r)
    d = degree_bound_f(p, q, r, k)
    print(d)
    if d < 0:
        raise ValueError()
    f = compute_f(p, q, r, d, k)
    s = compute_s(a, p, r, f, k)
    return s

from sympy import Integer
k = Symbol('k', integer=True)
a = (2 ** k) * (k ** 8) * (3 ** k)
print(gosper_sum(a, k))