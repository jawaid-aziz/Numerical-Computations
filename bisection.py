import math


def bisection_method(func, a, b, tolerance=1e-6, max_iterations=100):

    fa = func(a)
    fb = func(b)

    if fa * fb >= 0:
        raise ValueError(
            "Bisection method requires f(a) and f(b) to have opposite signs."
        )

    print(
        f"{'Iter':<6}{'a':<15}{'b':<15}{'Midpoint':<15}"
        f"{'f(midpoint)':<20}{'Abs Error':<15}"
    )
    print("-" * 86)

    prev_mid = None

    for iteration in range(1, max_iterations + 1):

        c = (a + b) / 2
        fc = func(c)

        if prev_mid is None:
            error = float("nan")
        else:
            error = abs(c - prev_mid)

        print(
            f"{iteration:<6}"
            f"{a:<15.6f}"
            f"{b:<15.6f}"
            f"{c:<15.6f}"
            f"{fc:<20.6f}"
            f"{error:<15.6f}"
        )

        if abs(fc) < tolerance:
            return c, iteration

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

        if prev_mid is not None and error < tolerance:
            return c, iteration

        prev_mid = c

    return c, max_iterations


if __name__ == "__main__":

    expression = input("Enter f(x): ")

    def f(x):
        return eval(
            expression,
            {
                "__builtins__": None,
                "x": x,
                "sin": math.sin,
                "cos": math.cos,
                "tan": math.tan,
                "log": math.log,
                "exp": math.exp,
                "sqrt": math.sqrt,
                "pi": math.pi,
                "e": math.e,
            },
        )

    a = float(input("Enter lower bound a: "))
    b = float(input("Enter upper bound b: "))
    tolerance = float(input("Enter tolerance ε: "))
    max_iterations = int(input("Enter maximum iterations: "))

    root, steps = bisection_method(
        f,
        a,
        b,
        tolerance,
        max_iterations,
    )

    print("\nResult")
    print(f"Root = {root:.10f}")
    print(f"Iterations = {steps}")