import math

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

# Create Forward Difference Table
diff_table = [[0 for _ in range(n)] for _ in range(n)]

for i in range(n):
    diff_table[i][0] = y[i]

for j in range(1, n):
    for i in range(n - j):
        diff_table[i][j] = diff_table[i + 1][j - 1] - diff_table[i][j - 1]

# Display Difference Table
print("\nForward Difference Table:")

width = 12

headers = ["x", "y"] + [f"Δ^{i}y" for i in range(1, n)]

for h in headers:
    print(f"{h:<{width}}", end="")
print()

print("-" * width * len(headers))

for i in range(n):
    print(f"{x[i]:<{width}.4f}", end="")

    for j in range(n):
        if j < n - i:
            print(f"{diff_table[i][j]:<{width}.4f}", end="")
        else:
            print(" " * width, end="")

    print()

# Interpolation point
xp = float(input("\nEnter the value of x for interpolation: "))

# Calculate h
h = x[1] - x[0]

# Calculate p
p = (xp - x[0]) / h

# Newton Forward Formula
result = diff_table[0][0]
p_term = 1

for i in range(1, n):
    p_term *= (p - (i - 1))
    result += (p_term * diff_table[0][i]) / math.factorial(i)

# Output result
print(f"\nInterpolated value at x = {xp} is:")
print(f"f({xp}) = {result:.6f}")