
using ModelingToolkit, DifferentialEquations, Plots

# ------------------------------------------------
# Dynamic Modeling of WES ------------------------
# ------------------------------------------------

#=Step 1 - define time as a parameter and the differential
operator as an ODE.=#

@parameters t
D = Differential(t)

#=Parameters. Note: must define paramater "t" before defining the state variables (thus defined above), otherwise julia will 
complain about undefined parameter 't'. Time, t, can also be defined as a variable, but the distinction is not clear. These
are put into named vectors for convenience. Note that the default parameter values are set in the declaration. This is not 
required, it is just convenient. See comment below on vars=# 

parspop = @parameters r₁=0.038 Kp₁=1.5 r₂=0.042 r₁₂=0 r₂₁=0 Kp₂=9.7 ml₁₂=10 ml₂₁=10 mp₁₂=10 mp₂₁=10  #population/migration dynamics
parsecn = @parameters δ₁=0.05 δ₂=0.05 α₁=0.5 α₂=0.5 s₁=0.25 s₂=0.21 Cₘ₁=0.7 Cₘ₂=0.7 a₁=2.7 a₂=1.7 dmₓ=100 #economic dynamics
parsear = @parameters re₁=0.1 re₂=0.1 eb=0.00004 u=0.0025 α=0.1 σ₁=0.03 σ₂=0.03 G₀=20 G₁=5 Gₘ=20 Tg=50 η=1 #earth system and climate 
pars = [parspop;parsecn;parsear] #combine into a single vector.

#=define the state variables. Assign them to the name "vars" for convenience. Again, initial values are set in the 
declaration. Alternatively, they can be set as follows:

u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]

Then the ODEProblem function can be called with u0.=#

vars = @variables P₁(t)=0.24 P₂(t)=0.24 K₁(t)=0 K₂(t)=0 G(t)=2.8 e₁(t)=0.0004 e₂(t)=0.0004 z(t)=0

#=below is a differentiable threshold function that goes from 0
to 1 as its argument goes from negative to positive. This keeps
things well-behaved and physically meaninful (no negative stocks).
This could be done in many different ways....=#

θ₁(x,p,h) = (max(0,x))^p/(h^p + (max(0,x))^5)
@register θ₁(x,p,h)

θ(x) = θ₁(x,5,0.0001)
@register θ(x)


#=this is a differentiable form of max(x,0) - i.e. approximates it very 
closely with no kink at zero.=#

ϕ(x) = x*θ(x)
@register ϕ(x)

#=economic subsystem. Society may use industrial technology
or a backstop technology that requires only labor inputs.=#

#=industrial output regime - Cobb-Douglass Technology with constant
returns to scale.=#

Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α))
@register Yᵢ(a,α,K,P)

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)


#= background production regime (labor based production)
is given by ab₁P₁ but we assume the multiplier is 1. We then
use the threshold function to switch to industrical production
only when it is more productive than the backdrop technology.=#

Y(P,Y) = P + ϕ(Y-P)
@register Y(P,K)

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)

#=emissions from backdrop or industrialized. Assume eb=0.1*e1 in the 
scenario below, but this is arbitrary. The 'r' in er refers to 
'realized' caron intensity depending on which technology is used=#

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂

#decarbonization. e1 and e2 obey simple odes with exponential decay. z is a Markov
# variable that turns on decarobonization (represents a policy change
# to begin decarbonization) at a certian level of the global
# externality. It is zero until G>Tg (Tg = threshold to begin
# decarbonization. The following function is a fudge - doesn't work because
# after decarbonization begins, when G falls below Tg, decarbonization stops.
# In the XPPAUT model I used a Markov transition matrix to just turn decarbonization
# on and leave it on.  Not sure how to fix in Julia. One way to do it is
# to make z a differential equation - i.e. and integral controller.

# dz/dt= η*θ(G-Tg)(1-Z) where η is the rate at which decarbonization turns on.


#=investment depends on disposable income... which is income
beyond a subsistence level Cₘ=#

Yd₁ = ϕ(Y₁ - P₁*Cₘ₁)
Yd₂ = ϕ(Y₂ - P₂*Cₘ₂)



#= The earth system model. This is the lowest order reasonable 
representation of atmospheric carbon dynamcis given by a single
'externality' - i.e. the carbon stock in the atmoshpere which has a 
component that follows a natural source-sink dynamic that equilbrates
at 280 ppm (pre-industrial level) given by  u*(2.8-G) and 
a nonlinear threshold component, i.e. large scale melting
of ice sheets, rapid release of green house gases from the soil, etc.
given by α*θ(G-G₀). Finally, industrial production emits into
the global carbon stock=#


#= Damage Function(s). We assume that there is the possibility of
a 'press' disturbance and a discontinous disturbance or 'shock'.
For the continous portion we assume  damage = σ₁*exp(G-G₁)*K₁, i.e.
damages grow expoentially beyond threshold G₁.  The next stage of
development is to model the stochastic portion.=#

dm₁ = σ₁*dmₓ*exp(G-G₁)/(dmₓ + exp(G-G₁))
dm₂ = σ₂*dmₓ*exp(G-G₁)/(dmₓ + exp(G-G₁))


#= Popuation dynamics.  Local population is logistic growth, migration 
is discrete choice with logit. First cut - suppose that migration depends on ability 
to migrate (per capita income) and relative welfare perceptions between regions 
(ratio of per capita income.). From the discrete choice framework we can define 
"pressure to migrate" from region i to j as mᵢⱼ as defining the log odds of migrating
(a discrete choice) with mₒ being the half saturation point of the logistic function =#

#define standard logistic

λ(x) = exp(x)/(exp(x)+1)
@register λ(x)

# Then for the migration pressure we get (various attempts and possibilities):

#m₁₂ = m₁*(Y₁/P₁ + Y₂*P₁/(P₂*Y₁) - mₒ)
#m₂₁ = m₁*(Y₂/P₂ + Y₁*P₂/(P₁*Y₂) - mₒ)

#m₁₂ = m₁*(Y₂*P₁/(P₂*Y₁) - mₒ)
#m₂₁ = m₁*(Y₁*P₂/(P₁*Y₂) - mₒ)

#define a function to calculate the percentage difference between to numbers
#for code readability

pd(a,b) = (a-b)/(a+b)
@register pd(a,b)

#not for the migration pressures....

m₁₂ = mp₁₂*pd(Y₂/P₂,Y₁/P₁) - ml₁₂
m₂₁ = mp₂₁*pd(Y₁/P₁,Y₂/P₂) - ml₂₁


# and for the probability of migration, pᵢⱼ, we get

p₁₂ = λ(m₁₂)
p₂₁ = λ(m₂₁)

#and the rate at which these decisions are made is given by rᵢⱼ to
#turn random events into migration rates. 

#differential equations

eqs = [
    D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - r₁₂*P₁*p₁₂ + r₂₁*P₂*p₂₁, #human population dynamics 
    D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + r₁₂*P₁*p₁₂ - r₂₁*P₂*p₂₁, # standard logistic.
    D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - dm₁*K₁, #change in capital stock = 
    D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - dm₂*K₂, # savings-depreciation-damages
    D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*θ(G-G₀)*θ₁(Gₘ-G,3,1), # global externality.
    D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
    D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
    D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
    ]

@named sys = ODESystem(eqs,t,vars,pars);

#= the variable sys is a symbolic object containing a mathematical representation
of the problem. Next, translate the symbolic problem into code to be solved numerically.
using ODEProblem(systemname,initial conditions, integration time span, parameters. 
Note that since we havae set the default values for vars and pars in sys, we send in
blank placeholders in the form of [].  So we create the problem, then sovle it.=#

prob = ODEProblem(sys, [], (0.0,600), [])
sol = solve(prob,saveat=1.0)

#=Now plot the results using Plots with the GR backend.  This is the simplest but does not
produce as high quality or flexible output as do other backends.=#

plot(sol, idxs=[Y₁]);
plot!(sol, idxs=[Y₂]);
p1=plot!(sol, idxs=[10*G],xlabel="Time", ylabel="GNI");
plot(sol,idxs=[P₁]);
p2=plot!(sol,idxs=[P₂],xlabel="Time", ylabel="Population"); 
p3=plot(sol[Y₁/P₁],sol[Y₂/P₂], xlabel="per capita GNI, region 1",
ylabel="per capita GNI, region 2");
p4=plot(sol[G*100],sol[Y₁+Y₂], xlabel="per capita GNI, region 1",
ylabel="per capita GNI, region 2");
plot(p2, p1, p3, p4, layout=(2,2), legend=false, size=(750,500))
