"""Solve a system of linear equations using LU decomposition (Doolittle method).
 
The program accepts equations in natural form like:
    2x - y + z = 2
    3x + 3y + 9z = -1
    3x + 3y + 5z = 4
 
It extracts the matrix A and vector b automatically, then solves using:
    L y = b  (forward substitution)
    U x = y  (backward substitution)
"""
 
import re
 
 
def print_matrix(matrix, title):
    """Display a matrix in a readable format."""
    print(f"\n{title}")
    for row in matrix:
        print("  ", "\t".join(f"{value:12.6f}" for value in row))
 
 
def print_vector(vector, title):
    """Display a vector in a readable format."""
    print(f"\n{title}")
    for index, value in enumerate(vector, start=1):
        print(f"  {index:>2}: {value:12.6f}")
 
 
# ──────────────────────────────────────────────
# Equation Parser
# ──────────────────────────────────────────────
 
def parse_equations(equation_strings):
    """Parse a list of equation strings into matrix A and vector b.
 
    Supports formats like:
        2x - y + z = 2
        -3x + 2.5y - z = 0
        x + 4z = 7           (missing variables get coefficient 0)
 
    Returns:
        variables : list of variable names in order of first appearance
        A         : coefficient matrix (list of lists)
        b         : right-hand side vector (list)
    """
    # ── Step 1: collect every variable name across all equations ──
    variables = []
    seen = set()
    for eq in equation_strings:
        lhs = eq.split("=")[0]
        for match in re.finditer(r"[a-zA-Z_][a-zA-Z0-9_]*", lhs):
            name = match.group()
            if name not in seen:
                variables.append(name)
                seen.add(name)
 
    n = len(variables)
    var_index = {v: i for i, v in enumerate(variables)}
 
    # ── Step 2: parse each equation ──
    A = []
    b = []
 
    for eq in equation_strings:
        eq = eq.strip()
        if "=" not in eq:
            raise ValueError(f"No '=' found in equation: '{eq}'")
 
        lhs, rhs = eq.split("=", 1)
        b.append(float(rhs.strip()))
 
        coeffs = [0.0] * n
 
        # Normalise the LHS so the regex always sees an explicit sign.
        # Insert '+' at the start if the first token isn't signed.
        lhs = lhs.strip()
        if lhs and lhs[0] not in "+-":
            lhs = "+" + lhs
 
        # Pattern: optional sign, optional number, variable name
        # e.g.  +2x   -y   +3.5z   +x
        token_re = re.compile(
            r"([+-])\s*(\d+\.?\d*)?\s*\*?\s*([a-zA-Z_][a-zA-Z0-9_]*)"
        )
 
        for m in token_re.finditer(lhs):
            sign_char, coeff_str, var_name = m.groups()
            sign = 1.0 if sign_char == "+" else -1.0
            coeff = float(coeff_str) if coeff_str else 1.0
            coeff *= sign
 
            if var_name not in var_index:
                raise ValueError(
                    f"Variable '{var_name}' appears only in one equation. "
                    "All variables must appear in all equations (or be 0)."
                )
            coeffs[var_index[var_name]] = coeff
 
        A.append(coeffs)
 
    if len(A) != n:
        raise ValueError(
            f"System is not square: {len(A)} equations but {n} variables {variables}."
        )
 
    return variables, A, b
 
 
# ──────────────────────────────────────────────
# LU Decomposition (Doolittle)
# ──────────────────────────────────────────────
 
def lu_decomposition_doolittle(a):
    """Return L and U such that A = L * U using the Doolittle method."""
    n = len(a)
    l = [[0.0] * n for _ in range(n)]
    u = [[0.0] * n for _ in range(n)]
 
    for i in range(n):
        for k in range(i, n):
            total = sum(l[i][j] * u[j][k] for j in range(i))
            u[i][k] = a[i][k] - total
 
        if abs(u[i][i]) < 1e-12:
            raise ValueError(
                "Zero pivot encountered during LU decomposition. "
                "The matrix may be singular or require pivoting."
            )
 
        l[i][i] = 1.0
        for k in range(i + 1, n):
            total = sum(l[k][j] * u[j][i] for j in range(i))
            l[k][i] = (a[k][i] - total) / u[i][i]
 
    return l, u
 
 
def forward_substitution(l, b):
    """Solve L y = b for y (L is lower triangular with unit diagonal)."""
    n = len(l)
    y = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(l[i][j] * y[j] for j in range(i))
    return y
 
 
def backward_substitution(u, y):
    """Solve U x = y for x (U is upper triangular)."""
    n = len(u)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(u[i][i]) < 1e-12:
            raise ValueError("Division by zero in backward substitution.")
        x[i] = (y[i] - sum(u[i][j] * x[j] for j in range(i + 1, n))) / u[i][i]
    return x
 
 
# ──────────────────────────────────────────────
# Input helpers
# ──────────────────────────────────────────────
 
def read_equations_from_user():
    """Read equations one per line from the keyboard."""
    print("\nEnter your equations one per line.")
    print("Use variable names like x, y, z (or x1, x2, x3, ...).")
    print("Example:  2x - y + z = 2")
    print("Type a blank line when finished.\n")
 
    equations = []
    while True:
        line = input(f"  Eq {len(equations) + 1}: ").strip()
        if line == "":
            if len(equations) == 0:
                print("  Please enter at least one equation.")
                continue
            break
        equations.append(line)
 
    return equations
 
 
def predefined_test_case():
    """A small example for quick verification."""
    return [
        "2x - y + z = 2",
        "3x + 3y + 9z = -1",
        "3x + 3y + 5z = 4",
    ]
 
 
def validate_square_matrix(a):
    n = len(a)
    for row in a:
        if len(row) != n:
            raise ValueError("Parsed matrix is not square — check your equations.")
 
 
# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
 
def main():
    print("=" * 50)
    print("  LU Decomposition — Doolittle Method")
    print("=" * 50)
    print("1. Enter equations manually")
    print("2. Use predefined test case")
 
    choice = input("\nChoose (1 or 2): ").strip()
 
    if choice == "2":
        equation_strings = predefined_test_case()
        print("\nUsing equations:")
        for eq in equation_strings:
            print(f"  {eq}")
    else:
        equation_strings = read_equations_from_user()
 
    variables, A, b = parse_equations(equation_strings)
    validate_square_matrix(A)
 
    print(f"\nDetected variables: {variables}")
 
    l, u = lu_decomposition_doolittle(A)
    y = forward_substitution(l, b)
    x = backward_substitution(u, y)
 
    print_matrix(A, "Coefficient Matrix A")
    print_matrix(l, "Lower Triangular Matrix L")
    print_matrix(u, "Upper Triangular Matrix U")
    print_vector(y, "Intermediate Vector y  (from Ly = b)")
    print_vector(x, f"Solution Vector x  (from Ux = y)")
 
    print("\n── Solution ──")
    for var, val in zip(variables, x):
        print(f"  {var} = {val:.6f}")
 
 
if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print(f"\nError: {error}")
    except ZeroDivisionError:
        print("\nError: Division by zero encountered.")