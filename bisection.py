def bisection_method(func, a, b, tolerance=1e-6, max_iterations=1000):
	fa = func(a)
	fb = func(b)

	if fa * fb >= 0:
		raise ValueError("Bisection method requires f(a) and f(b) to have opposite signs.")

	iterations = 0

	while abs(b - a) > tolerance and iterations < max_iterations:
		c = (a + b) / 2
		fc = func(c)

		if fc == 0:
			return c, iterations + 1

		if fa * fc < 0:
			b = c
			fb = fc
		else:
			a = c
			fa = fc

		iterations += 1

	return (a + b) / 2, iterations


if __name__ == "__main__":
	def f(x):
		return x**3 - x - 2


	root, steps = bisection_method(f, 1, 2, tolerance=1e-6)
	print("Root:", root)
	print("Iterations:", steps)