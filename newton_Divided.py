# Newton Divided Difference Interpolation

# Input number of data points
n = int(input("Enter the number of data points: "))

# Input x values
x = []
print("\nEnter x values:")
for i in range(n):
    x.append(float(input(f"x[{i}] = ")))

# Input y values
y = []
print("\nEnter y values:")
for i in range(n):
    y.append(float(input(f"y[{i}] = ")))

# Create Divided Difference Table
table = [[0 for _ in range(n)] for _ in range(n)]

# First column contains y values
for i in range(n):
    table[i][0] = y[i]

# Calculate divided differences
for j in range(1, n):
    for i in range(n - j):
        table[i][j] = (
            table[i + 1][j - 1] - table[i][j - 1]
        ) / (x[i + j] - x[i])

# Display Table
print("\nNewton Divided Difference Table:")

width = 15

headers = ["x", "f(x)"] + [f"Order {i}" for i in range(1, n)]

for h in headers:
    print(f"{h:<{width}}", end="")
print()

print("-" * width * len(headers))

for i in range(n):

    print(f"{x[i]:<{width}.4f}", end="")

    for j in range(n):

        if j < n - i:
            print(f"{table[i][j]:<{width}.6f}", end="")
        else:
            print(" " * width, end="")

    print()

# Interpolation point
xp = float(input("\nEnter the value of x for interpolation: "))

# Newton Divided Difference Formula
result = table[0][0]

product_term = 1

for i in range(1, n):

    product_term *= (xp - x[i - 1])

    result += product_term * table[0][i]

# Output
print(f"\nInterpolated value at x = {xp} is:")
print(f"f({xp}) = {result:.6f}")