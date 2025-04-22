using ModelingToolkit
#using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations

@parameters t
D = Differential(t)

pars = @parameters r₁ = 0.038
vars = @variables x1(t) = 1.0 x2(t) = 2.0

eqs = [
    D(x1) ~ r₁ * x1
    D(x2) ~ r₁ * x2
]

@named sys = ODESystem(eqs, t, vars, pars);
simplified_sys = structural_simplify(sys)
prob = ODEProblem(sys, [], (0.0, 600), [])

#@mtkmodel FOL begin
#    @parameters begin
#        C = [1.0 1.0] #[1.0 2.0; 1.0 1.0] #ones(2, 2)
#    end
#    @variables begin
#        x1(t) = 1
#        x2(t) = 2
#    end
#    @equations begin
#        D(x1) ~ x1
#        D(x2) ~ x2
#        #D(x1) ~ C[1, 1] * x1 + C[1, 2] * x2
#        #D(x2) ~ C[2, 2] * x2 + C[2, 1] * x1
#    end
#end

#@mtkbuild fol = FOL()
#prob = ODEProblem(fol, [], (0.0, 100.0), [])