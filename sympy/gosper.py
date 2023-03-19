from sympy import fraction, combsimp, cancel, Symbol


def compute_bc(a, k):
    x = a.subs(k, k + 1) / a
    x = combsimp(x)
    x = cancel(x)
    b, c = fraction(x)
    # TODO: ensure b, c are polys in k
    return b, c


def gosper_sum(a, k):
    b, c = compute_bc(a, k)