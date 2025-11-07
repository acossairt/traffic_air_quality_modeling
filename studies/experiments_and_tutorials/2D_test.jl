using ModelingToolkit
using DifferentialEquations

@parameters t
D = Differential(t)

my_r = [1.5 2.0 2.5 3.0]
n = length(my_r)

r = @parameters [Symbol("r_$i") for i in 1:n]...  # creates r_1, r_2, ..., r_n
@variables x1(t) = 1.0 x2(t) = 2.0

eqs = [
    D(x1) ~ r[1] * x1
    D(x2) ~ r[2] * x2
]

@named sys = ODESystem(eqs, t);
fol = structural_simplify(sys)
prob = ODEProblem(fol, [], (0.0, 10.0), [])

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