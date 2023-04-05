from sympy import Symbol, Integer, factorial, binomial
from gosper import gosper_sum


def main():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)
    expr = Integer(1)
    print(gosper_sum(expr, k))


if __name__ == "__main__":
    main()