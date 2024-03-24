from scipy.optimize import linprog
import numpy as np

# Define the coefficients of the objective function (we set these to 0
# because we're not optimizing any specific function)
c = [0, 0]

# Define the coefficients of the inequality constraints (none in this case)
A_ub = None
b_ub = None

# Define the coefficients of the equality constraints (AX = B)
A_eq = np.array([[1, 2], [3, 4]])
b_eq = np.array([5, 11])

# Define the bounds for each variable
x0_bounds = (0, None) # a >= 0
x1_bounds = (0, None) # b >= 0

# Solve the problem
res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[x0_bounds, x1_bounds], method='highs')

# Print the results
print(f'Status: {res.message}')
if res.success:
    print(f'Solution: a = {res.x[0]}, b = {res.x[1]}')
else:
    print('No solution found that satisfies the constraints.')