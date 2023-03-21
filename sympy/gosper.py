from sympy import fraction, combsimp, cancel, resultant, roots, Symbol, Dummy
from sympy import gcd, degree, linsolve, simplify


DEBUG_MODE = 1


def debug(s):
    if DEBUG_MODE:
        print(s)


def compute_bc(a, k):
    debug(f"called compute_bc(a := {a}, k := {k})")
    x = a.subs(k, k + 1) / a
    x = combsimp(x)
    x = cancel(x)
    b, c = fraction(x)
    # TODO: ensure b, c are polys in k
    debug(f"returned (b := {b}, c := {c}) from compute_bc")
    return b, c


def compute_pqr(b, c, k):
    # TODO: precons
    debug(f"called compute_pqr(b := {b}, c := {c}, k := {k})")
    j = Dummy("j")
    rp = resultant(b, c.subs(k, k+j), k)
    debug(f"resultant is {rp}")
    rz = roots(rp, j, multiple=True, filter='Z')
    debug(f"resultant roots are {rz}")

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
    debug(f"returned (p := {p}, q := {q}, r := {r}) from compute_pqr")
    return p, q, r


def degree_bound_f(p, q, r, k):
    debug(f"called degree_bound_f(p := {p}, q := {q}, r := {r}, k := {k})")
    x = q + r.subs(k, k-1)
    y = q - r.subs(k, k-1)
    n = degree(x, k)
    if n <= degree(y, k):
        d = degree(p.subs(k, k-1), k) - degree(y, k)
        debug(f"returned d := {d} from degree_bound_f case 1")
        return d

    a = x.coeff(k, n)
    b = y.coeff(k, n - 1)
    assert not a.equals(0)

    z = (-2 * b) / a
    debug(f"found z := {z}")
    # TODO: needs to check for int literal, not just integer
    if not z.is_integer or z < 0 or z.free_symbols:
        d = degree(p.subs(k, k-1), k) - n + 1
        debug(f"returned d := {d} from degree_bound_f case 2")
        return d
    
    d = max(z, degree(p.subs(k, k-1), k) - n + 1)
    debug(f"returned d := {d} from degree_bound_f case 3")
    return d


def compute_f(p, q, r, d, k):
    debug(f"called compute_f(p := {p}, q := {q}, r := {r}, d := {d}, k := {k})")
    f = 0
    cs = []
    for i in range(0, d + 1):
        c = Dummy(f"c{i}")
        f += c * (k ** i)
        cs.append(c)
    
    e = q * f - r.subs(k, k-1) * f.subs(k, k-1) - p.subs(k, k-1)

    es = []
    for i in range(0, 2 * (d + 1)): # TODO: idk why we need * 2
        es.append(e.subs(k, i))
    
    debug(f"built system of equations es := {es}")
    
    ss = linsolve(es, cs)
    debug(f"found solution set ss := {ss}")

    if len(ss) != 1:
        raise ValueError()

    s = next(iter(ss))
    f = f.subs(list(zip(cs, s)))
    f = f.subs(list(zip(cs, [0] * len(cs))))
    debug(f"returned f := {f} from compute_f")
    return f


def compute_s(a, p, r, f, k):
    debug(f"called compute_s(a := {a}, p := {p}, r := {r}, f := {f}, k := {k})")
    s = (r.subs(k, k-1) / p.subs(k, k-1)) * f.subs(k, k-1) * a
    debug(s)
    # s = cancel(s.subs(k, k+1))
    debug(f"returned s := {s} from compute_s")
    return s


def gosper_sum(a, k):
    b, c = compute_bc(a, k)
    p, q, r = compute_pqr(b, c, k)
    d = degree_bound_f(p, q, r, k)
    if d < 0:
        raise ValueError()
    f = compute_f(p, q, r, d, k)
    s = compute_s(a, p, r, f, k)
    # NOTE: s is indefinite (s - s.subs(k, 0) is the real answer)
    return s
