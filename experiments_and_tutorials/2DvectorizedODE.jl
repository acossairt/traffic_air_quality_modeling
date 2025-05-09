using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

#Build the model symbolically

#number of patches,

N=2

@variables p(t)[1:N]
@parameters c[1:N,1:N]

#define equations

eqs = [
    D.(p) ~ collect(c*p)
        ]

#build model symbolically
@mtkbuild model = ODESystem(eqs, t)

# ===================================================
# build numerical problem to solve
# set parameter matrix cm and set initial conditions

cm=[-0.1 0.5;
    -0.5 -0.1]

prob = ODEProblem(model, [p[1]=>0.8,p[2]=>0.5], (0.0, 30), [c => cm])

#show egenvalues

eigvals(cm)

#= numerically solve the problem  "prob" and plot it - state 
variables as a function of time =#

sol = solve(prob,Tsit5())
plot(sol, idxs = (p[1]), ylim=(-1,1))
plot!(sol, idxs = (p[2]), xlim=(0,30))

# Then plot it in state space

plot(sol, idxs = (p[1],p[2]), xlim=(-1,1), ylim=(-1,1))