from sympy import fraction, combsimp, cancel, resultant, roots, Symbol, Dummy
from sympy import gcd, degree, linsolve, simplify


DEBUG_MODE = 1


def debug(s: str) -> None:
    """
    Prints debug information to stdout if DEBUG_MODE is set.

    Parameters
    ----------
        s: The string to print.
    """
    if DEBUG_MODE:
        print(s)


def compute_bc(a, k):
    """
    Computes polynomials $b_k$ and $c_k$ such that $a_{k+1} / a_k = b_k / c_k$. If no
    such polynomials can be found, then $b_k = a_{k+1} / a_k$ and $c_k = 1$.

    Parameters
    ----------
        a: The summand $a_k$.
        k: The index variable $k$.

    Returns
    -------
        The polynomials $b_k$ and $c_k$. If the procedure failed, then $b_k$ will be 
        non-polynomial and $c_k$ will be $1$.
    """
    debug(f"called compute_bc(a := {a}, k := {k})")
    x = a.subs(k, k + 1) / a
    x = combsimp(x)
    x = cancel(x)
    b, c = fraction(x)
    # TODO: Ensure that b and c are polynomials in k.
    debug(f"returned (b := {b}, c := {c}) from compute_bc")
    return b, c


def compute_pqr(b, c, k):
    """
    Computes polynomials $p_{k+1}, q_{k+1}$ and $r_{k+1}$ such that $b_k / c_k = 
    (p_{k+1} / p_k) * (q_{k+1} * r_{k+1})$ and $\\mathrm{gcd}(q_k, r_{k+j}) = 1$ for all
    integers $j \\geq 0$.

    Parameters
    ----------
        b: The polynomial $b_k$.
        c: The polynomial $c_k$.
        k: The index variable $k$.
    
    Returns
    -------
        The polynomials $p_{k+1}, q_{k+1}$ and $r_{k+1}$. 
    """
    # TODO: Ensure that b and c are polynomials in k.
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

    p = cancel(p)
    q = cancel(q)
    r = cancel(r)

    # TODO: Ensure that p, q and r satisfy the conditions. 
    debug(f"returned (p := {p}, q := {q}, r := {r}) from compute_pqr")
    return p, q, r


def degree_bound_f(p, q, r, k):
    """
    Computes an upper bound for the degree of $f_k$.

    Parameters
    ----------
        p: The polynomial $p_{k+1}$.
        q: The polynomial $q_{k+1}$.
        r: The polynomial $r_{k+1}$.
        k: The index variable $k$.

    Returns
    -------
        The upper bound on the degree of $f_k$. The upper bound will be $0$ if $f_k$ is
        $0$ and negative if $f_k$ cannot be computed by Gosper's algorithm.
    """
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
    if not z.is_integer or z < 0 or z.free_symbols:
        d = degree(p.subs(k, k-1), k) - n + 1
        debug(f"returned d := {d} from degree_bound_f case 2")
        return d
    
    d = max(z, degree(p.subs(k, k-1), k) - n + 1)
    debug(f"returned d := {d} from degree_bound_f case 3")
    return d


def compute_f(p, q, r, d, k):
    """
    Computes the polynomial $f_k$ such that $p_k = q_{k+1} f_k - r_k f_{k-1}$. If $f_k$
    cannot be found, a ValueError will be raised.

    Parameters
    ----------
        p: The polynomial $p_{k+1}$.
        q: The polynomial $q_{k+1}$.
        r: The polynomial $r_{k+1}$.
        d: The degree bound of $f_k$.
        k: The index variable $k$.
    
    Returns
    -------
        The polynomial $f_k$. 
    """
    debug(f"called compute_f(p := {p}, q := {q}, r := {r}, d := {d}, k := {k})")
    f = 0
    cs = []
    for i in range(0, d + 1):
        c = Dummy(f"c{i}")
        f += c * (k ** i)
        cs.append(c)
    
    e = q * f - r.subs(k, k-1) * f.subs(k, k-1) - p.subs(k, k-1)

    es = []
    for i in range(0, 2 * (d + 1)): # TODO: We shouldn't need to multiply by 2.
        es.append(e.subs(k, i))
    
    debug(f"built system of equations es := {es}")
    
    ss = linsolve(es, cs)
    debug(f"found solution set ss := {ss}")

    if len(ss) != 1:
        raise ValueError("Summand is not Gosper summable.")

    s = next(iter(ss))
    f = f.subs(list(zip(cs, s)))
    f = f.subs(list(zip(cs, [0] * len(cs))))
    debug(f"returned f := {f} from compute_f")
    return f


def compute_s(a, p, r, f, k):
    """
    Computes the indefinite sum $s_k$.

    Parameters
    ----------
        a: The summand $a_k$.
        p: The polynomial $p_{k+1}$.
        r: The polynomial $r_{k+1}$.
        f: The polynomial $f_k$.
        k: The index variable $k$.

    Returns
    -------
        The indefinite sum $s_k$.
    """
    debug(f"called compute_s(a := {a}, p := {p}, r := {r}, f := {f}, k := {k})")
    s = (r.subs(k, k-1) / p.subs(k, k-1)) * f.subs(k, k-1) * a
    debug(f"returned s := {s} from compute_s")
    return s


def gosper_sum(a, k):
    """
    Computes the indefinite sum of $a_k$ with Gosper's algorithm. If the summand is not
    Gosper summable, a ValueError will be raised.

    Parameters
    ----------
        a: The summand $a_k$.
        k: The index variable $k$.

    Returns
    -------
        The indefinite sum $s_k$.
    """
    b, c = compute_bc(a, k)
    p, q, r = compute_pqr(b, c, k)
    d = degree_bound_f(p, q, r, k)
    if d < 0:
        raise ValueError("Summand is not Gosper summable.")
    f = compute_f(p, q, r, d, k)
    s = compute_s(a, p, r, f, k)
    return s
