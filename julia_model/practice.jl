using ModelingToolkit
using DifferentialEquations

@parameters t
D = Differential(t)

n = 3  # number of variables and parameters

@parameters r[1:n]           # creates r[1], r[2], ..., r[n]
@variables (x(t))[1:n]       # creates x[1](t), x[2](t), ..., x[n](t)

# Differential equations: dxᵢ/dt = rᵢ * xᵢ
eqs = [D(x[i]) ~ r[i] * x[i] for i in 1:n]

@named sys = ODESystem(eqs, t)
sys_simplified = structural_simplify(sys)

# Initial conditions
u0 = [1.0 for _ in 1:n]

# Parameter values (must match order of r)
p_vals = [0.5 + 0.1 * i for i in 1:n]

prob = ODEProblem(sys_simplified, u0, (0.0, 10.0), p_vals)
#sol = solve(prob)