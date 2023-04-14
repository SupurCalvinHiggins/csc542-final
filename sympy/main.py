from sympy import Symbol, Integer, Rational, factorial, binomial, simplify
from gosper import gosper_sum


def main():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)
    expr = k
    print(gosper_sum(expr, k))


if __name__ == "__main__":
    main()