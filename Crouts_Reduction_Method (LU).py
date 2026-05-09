import numpy as np
import re

# -----------------------------
# Function: Crout Decomposition
# -----------------------------
def crout_decomposition(A):
    n = len(A)

    L = np.zeros((n, n))
    U = np.identity(n)

    for j in range(n):

        # Compute L
        for i in range(j, n):
            sum1 = 0
            for k in range(j):
                sum1 += L[i][k] * U[k][j]

            L[i][j] = A[i][j] - sum1

        # Compute U
        for i in range(j + 1, n):
            sum2 = 0
            for k in range(j):
                sum2 += L[j][k] * U[k][i]

            if L[j][j] == 0:
                raise ZeroDivisionError("Division by zero!")

            U[j][i] = (A[j][i] - sum2) / L[j][j]

    return L, U


# -----------------------------
# Forward Substitution
# Solve Ly = b
# -----------------------------
def forward_substitution(L, b):
    n = len(L)
    y = np.zeros(n)

    for i in range(n):
        sum1 = 0
        for j in range(i):
            sum1 += L[i][j] * y[j]

        y[i] = (b[i] - sum1) / L[i][i]

    return y


# -----------------------------
# Backward Substitution
# Solve Ux = y
# -----------------------------
def backward_substitution(U, y):
    n = len(U)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        sum1 = 0
        for j in range(i + 1, n):
            sum1 += U[i][j] * x[j]

        x[i] = y[i] - sum1

    return x


# -----------------------------
# Input Section
# -----------------------------
n = int(input("Enter number of equations: "))

A = []
B = []

print("\nEnter equations like:")
print("2x + y - z = 8\n")

variables = ['x', 'y', 'z', 'w', 'a', 'b', 'c']

for i in range(n):

    eq = input(f"Equation {i+1}: ")

    left, right = eq.split('=')
    right = float(right.strip())

    coeffs = []

    for var in variables[:n]:

        # Find coefficient of variable
        pattern = r'([+-]?\s*\d*)' + var
        match = re.search(pattern, left.replace(' ', ''))

        if match:
            coeff = match.group(1)

            if coeff in ['', '+']:
                coeff = 1
            elif coeff == '-':
                coeff = -1
            else:
                coeff = int(coeff)

        else:
            coeff = 0

        coeffs.append(coeff)

    A.append(coeffs)
    B.append(right)

A = np.array(A, dtype=float)
B = np.array(B, dtype=float)

# -----------------------------
# Crout Method
# -----------------------------
L, U = crout_decomposition(A)

# Solve Ly = b
Y = forward_substitution(L, B)

# Solve Ux = y
X = backward_substitution(U, Y)

# -----------------------------
# Output
# -----------------------------
print("\nMatrix A:")
print(A)

print("\nMatrix L:")
print(L)

print("\nMatrix U:")
print(U)

print("\nSolution:")

for i in range(n):
    print(f"{variables[i]} = {X[i]}")