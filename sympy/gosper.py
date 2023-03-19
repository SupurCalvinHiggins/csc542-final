from sympy import fraction, combsimp, cancel, resultant, roots, Symbol, Dummy
from sympy import gcd

def compute_bc(a, k):
    x = a.subs(k, k + 1) / a
    x = combsimp(x)
    x = cancel(x)
    b, c = fraction(x)
    # TODO: ensure b, c are polys in k
    return b, c


def compute_pqr(b, c, k):
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
        for e in range(0, -z, -1):
            p *= g.subs(k, k + e)

    return p, q, r


def gosper_sum(a, k):
    b, c = compute_bc(a, k)
    p, q, r = compute_pqr(b, c, k)