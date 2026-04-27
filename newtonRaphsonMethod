import math
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application


def newton_raphson(f, f_prime, x0, tolerance=1e-6, max_iterations=1000):
    x_prev = x0

    print("\nIteration\tx_old\t\t\tf(x_old)\t\tf'(x_old)\t\tx_new\t\t\t|x_new - x_old|")

    for i in range(max_iterations):
        fx = f(x_prev)
        fpx = f_prime(x_prev)

        if abs(fpx) < 1e-12:
            raise ValueError("Derivative is zero or near-zero. Method fails.")

        x_next = x_prev - fx / fpx
        error = abs(x_next - x_prev)

        print(f"{i + 1}\t\t{x_prev:.10f}\t{fx:.10f}\t{fpx:.10f}\t{x_next:.10f}\t{error:.10f}")

        if abs(x_next) > 1e10:
            raise ValueError("Method is diverging. Try a different x0.")

        if error < tolerance:
            return x_next, i + 1

        x_prev = x_next

    return x_prev, max_iterations


if __name__ == "__main__":
    func_str = input("Enter f(x): ")

    # Auto-differentiate using sympy
    x = sp.Symbol('x')
    f_expr = parse_expr(func_str, local_dict={'x': x}, transformations=(standard_transformations + (implicit_multiplication_application,)))
    f_prime_expr = sp.diff(f_expr, x)

    print(f"f'(x) = {f_prime_expr}")

    f       = sp.lambdify(x, f_expr, modules=["math"])
    f_prime = sp.lambdify(x, f_prime_expr, modules=["math"])

    x0       = float(input("Enter initial guess x0: "))
    tol      = float(input("Tolerance (default 1e-6): ") or 1e-6)
    max_iter = int(input("Max iterations (default 1000): ") or 1000)

    try:
        root, steps = newton_raphson(f, f_prime, x0, tol, max_iter)
        print(f"\nRoot: {root}")
        print(f"Iterations: {steps}")
    except ValueError as e:
        print(f"\nError: {e}")