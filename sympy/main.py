from sympy import Symbol, Integer, Rational, factorial, binomial, simplify, expand
from gosper import gosper_sum


def exam():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)

    print("1a.")
    expr = Integer(1)
    s = gosper_sum(expr, k)
    expr = s.subs(k, n + 1) - s.subs(k, 1)
    s = gosper_sum(expr, n)
    expr = s.subs(n, n + 1) - s.subs(n, 1)
    expr = simplify(expand(expr))
    print(expr)
    print()

    print("1b.")
    expr = n - k + 1
    s = gosper_sum(expr, k)
    expr = s.subs(k, n + 1) - s.subs(k, 1)
    s = gosper_sum(expr, n)
    expr = s.subs(n, n + 1) - s.subs(n, 1)
    expr = simplify(expand(expr))
    print(expr)
    print()

    print("2.")
    expr = k * (k + 1)
    s = gosper_sum(expr, k)
    expr = s.subs(k, n + 1) - s.subs(k, 1)
    expr = simplify(expand(expr))
    print(expr)
    print()

    print("3.")
    expr = k * factorial(k)
    s = gosper_sum(expr, k)
    expr = s.subs(k, n + 1) - s.subs(k, 0)
    expr = simplify(expand(expr))
    print(expr)
    print()

    print("4.")
    expr = k * 2 ** k
    s = gosper_sum(expr, k)
    expr = s.subs(k, n) - s.subs(k, 0)
    expr = simplify(expand(expr))
    print(expr)
    print()


def main():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)
    expr = Integer(1)
    s = gosper_sum(expr, k)
    print(s.subs(k, n) - s.subs(k, 0))


if __name__ == "__main__":
    SHOW_EXAM = False
    
    if SHOW_EXAM:
        exam()
    else:
        main()
