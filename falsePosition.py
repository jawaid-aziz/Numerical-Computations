def false_position_method(func, a, b, tolerance=1e-6, max_iterations=1000):
    fa = func(a)
    fb = func(b)

    if fa * fb >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs.")

    c = a
    for i in range(max_iterations):
        c = (a * fb - b * fa) / (fb - fa)
        fc = func(c)

        if abs(fc) < tolerance or abs(b - a) < tolerance:
            return c, i + 1

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    return c, max_iterations


if __name__ == "__main__":
    import math

    func_str = input("Enter f(x): ")
    func = lambda x: eval(func_str, {"x": x, "math": math, **vars(math)})

    a = float(input("Enter a: "))
    b = float(input("Enter b: "))
    tol = float(input("Tolerance (default 1e-6): ") or 1e-6)
    max_iter = int(input("Max iterations (default 1000): ") or 1000)

    root, steps = false_position_method(func, a, b, tol, max_iter)
    print(f"\nRoot: {root}")
    print(f"Iterations: {steps}")