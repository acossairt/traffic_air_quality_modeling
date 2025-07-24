module WESR_std_model_mod

# use exports to let other modules know the name spaces
# against style recommendation they are defined below

using ModelingToolkit, DifferentialEquations

# ----- define the model params -----

#=Step 1 - define time as a parameter and the differential
operator as an ODE.=#

@parameters t
D = Differential(t)
export t

#=Parameters. Note: must define paramater "t" before defining the state variables (thus defined above), otherwise julia will 
complain about undefined parameter 't'. Time, t, can also be defined as a variable, but the distinction is not clear. These
are put into named vectors for convenience. =#

@parameters r₁ Kp₁ r₂ r₁₂ r₂₁ Kp₂ ml₁₂ ml₂₁ mp₁₂ mp₂₁  #population/migration dynamics
export r₁, Kp₁, r₂, r₁₂, r₂₁, Kp₂, ml₁₂, ml₂₁, mp₁₂, mp₂₁
@parameters δ₁ δ₂ α₁ α₂ s₁ s₂ Cₘ₁ Cₘ₂ a₁ a₂ dmₓ #economic dynamics
export δ₁, δ₂, α₁, α₂, s₁, s₂, Cₘ₁, Cₘ₂, a₁, a₂, dmₓ
@parameters re₁ re₂ eb u α σ₁ σ₂ G₀ G₁ Gₘ Tg η #earth system and climate 
export re₁, re₂, eb, u, α, σ₁, σ₂, G₀, G₁, Gₘ, Tg, η

#=define the state variables. Assign them to the name "vars" for convenience. Again, initial values are set in the 
declaration. Alternatively, they can be set as follows:

u0=[P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]

Then the ODEProblem function can be called with u0.=#

@variables P₁(t) P₂(t) K₁(t) K₂(t) G(t) e₁(t) e₂(t) z(t)=0
export P₁, P₂, K₁, K₂, G, e₁, e₂, z

#=below is a differentiable threshold function that goes from 0
to 1 as its argument goes from negative to positive. This keeps
things well-behaved and physically meaninful (no negative stocks).
This could be done in many different ways....=#

θ₁(x,p,h) = (max(0,x))^p/(h^p + (max(0,x))^p)

θ(x) = θ₁(x,5,0.0001)
export θ

#=this is a differentiable form of max(x,0) - i.e. approximates it very 
closely with no kink at zero.=#

ϕ(x) = x*θ(x)
export ϕ

#=economic subsystem. Society may use industrial technology
or a backstop technology that requires only labor inputs.=#

#=industrial output regime - Cobb-Douglass Technology with constant
returns to scale.=#

Yᵢ(a,α,K,P) = a*((ϕ(K))^α)*((ϕ(P))^(1-α))
export Yᵢ

Yᵢ₁ = Yᵢ(a₁,α₁,K₁,P₁)
Yᵢ₂ = Yᵢ(a₂,α₂,K₂,P₂)


#= background production regime (labor based production)
is given by ab₁P₁ but we assume the multiplier is 1. We then
use the threshold function to switch to industrical production
only when it is more productive than the backdrop technology.=#

Y(P,Y) = P + ϕ(Y-P)
export Y

Y₁ = Y(P₁,Yᵢ₁)
Y₂ = Y(P₂,Yᵢ₂)
export Y₁, Y₂

#=emissions from backdrop or industrialized. Assume eb=0.1*e1 in the 
scenario below, but this is arbitrary. The 'r' in er refers to 
'realized' caron intensity depending on which technology is used=#

er₁ = θ(P₁-Yᵢ₁)*eb + θ(Yᵢ₁-P₁)*e₁
er₂ = θ(P₂-Yᵢ₂)*eb + θ(Yᵢ₂-P₂)*e₂
export er₁, er₂

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
export Yd₁, Yd₂


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
export dm₁, dm₂

#= Popuation dynamics.  Local population is logistic growth, migration 
is discrete choice with logit. First cut - suppose that migration depends on ability 
to migrate (per capita income) and relative welfare perceptions between regions 
(ratio of per capita income.). From the discrete choice framework we can define 
"pressure to migrate" from region i to j as mᵢⱼ as defining the log odds of migrating
(a discrete choice) with mₒ being the half saturation point of the logistic function =#

#define standard logistic

λ(x) = exp(x)/(exp(x)+1)
export λ

# Then for the migration pressure we get (various attempts and possibilities):

#m₁₂ = m₁*(Y₁/P₁ + Y₂*P₁/(P₂*Y₁) - mₒ)
#m₂₁ = m₁*(Y₂/P₂ + Y₁*P₂/(P₁*Y₂) - mₒ)

#m₁₂ = m₁*(Y₂*P₁/(P₂*Y₁) - mₒ)
#m₂₁ = m₁*(Y₁*P₂/(P₁*Y₂) - mₒ)

#define a function to calculate the percentage difference between to numbers
#for code readability

pd(a,b) = (a-b)/(a+b)
export pd

#not for the migration pressures....

m₁₂ = mp₁₂*pd(Y₂/P₂,Y₁/P₁) - ml₁₂
m₂₁ = mp₂₁*pd(Y₁/P₁,Y₂/P₂) - ml₂₁
export m₁₂, m₂₁

# and for the probability of migration, pᵢⱼ, we get

p₁₂ = λ(m₁₂)
p₂₁ = λ(m₂₁)
export p₁₂, p₂₁

#and the rate at which these decisions are made is given by rᵢⱼ to
#turn random events into migration rates. 

#climate thresholds

ct=θ(G-G₀)#climate threshold for carbon release
gmax=θ₁(Gₘ-G,3,1)#max carbon in the atmoshpere
export ct, gmax

#differential equations

eqs = [
D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - r₁₂*P₁*p₁₂ + r₂₁*P₂*p₂₁, #human population dynamics 
D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + r₁₂*P₁*p₁₂ - r₂₁*P₂*p₂₁, # standard logistic.
D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - dm₁*K₁, #change in capital stock = 
D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - dm₂*K₂, # savings-depreciation-damages
D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*ct*gmax, # global externality.
D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
]


function construct_ode_system()
    """This function is the main function to construct our model in ModelingToolkit."""

    @named ode = ODESystem(eqs, t);

    return ode

end # end of function


function construct_sde_system(ode, noiseeqs)

    @named sde = SDESystem(ode, noiseeqs)

    return sde

end # end of function


function construct_ode_problem(ode, tspan, u0, p)
    """
    run_time: e.g. (0.0,1000.0)
    initial_conditions: [P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    parameters: p=[r₁=>
    0.038, Kp₁=>1.5, #population growth, carrying capacity
    r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
    ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
    r₁₂=>0, r₂₁=>0, #migration rates
    α₁=>0.5, α₂=>0.5, #capital factor productivity
    a₁=>2.7, a₂=>1.7, #total factor productivity
    δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
    s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
    Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
    σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
    re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.5,  #decarb, #carbon intensities
    η=>1, #rate at which decarbonization is initiated
    u=>0.0025, α=>0.1, #earth system carbon uptake and release
    G₁=>5, G₀=>6.0, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
    ];
    """

    prob = ODEProblem(complete(ode), u0, tspan, p)

    return prob

end


function construct_sde_problem(sde, tspan, u0, p)
    """
    run_time: e.g. (0.0,1000.0)
    initial_conditions: [P₁=>0.24,P₂=>0.24,K₁=>0,K₂=>0,G=>2.8,e₁=>0.0004,e₂=>0.0004,z=>0]
    parameters: p=[r₁=>
    0.038, Kp₁=>1.5, #population growth, carrying capacity
    r₂=>0.042, Kp₂=>9.7, #population growth, carrying capacity
    ml₁₂=>10, ml₂₁=>10, mp₁₂=>10, mp₂₁=>10, #migration preference parameters
    r₁₂=>0, r₂₁=>0, #migration rates
    α₁=>0.5, α₂=>0.5, #capital factor productivity
    a₁=>2.7, a₂=>1.7, #total factor productivity
    δ₁=>0.05, δ₂=>0.05, #capital entropic decay rates  
    s₁=>0.25, s₂=>0.21, #savings rates (of disposable income)
    Cₘ₁=>0.7, Cₘ₂=>0.7, #minimum subsistence consumption
    σ₁=>0.03, σ₂=>0.03, dmₓ=> 100, #climate damages on infrastructure
    re₁=>0.1, re₂=>0.1, eb=>0.00004, Tg=>4.5,  #decarb, #carbon intensities
    η=>1, #rate at which decarbonization is initiated
    u=>0.0025, α=>0.1, #earth system carbon uptake and release
    G₁=>5, G₀=>6.0, Gₘ=>20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
    ];
    """
    global prob=SDEProblem(complete(sde), u0, tspan, p);

    return prob

end # end of function


#differential equations

eqs = [
D(P₁) ~ r₁*P₁*(1 - P₁/Kp₁) - r₁₂*P₁*p₁₂ + r₂₁*P₂*p₂₁, #human population dynamics 
D(P₂) ~ r₂*P₂*(1 - P₂/Kp₂) + r₁₂*P₁*p₁₂ - r₂₁*P₂*p₂₁, # standard logistic.
D(K₁) ~ s₁*Yd₁ - δ₁*K₁ - dm₁*K₁, #change in capital stock = 
D(K₂) ~ s₂*Yd₂ - δ₂*K₂ - dm₂*K₂, # savings-depreciation-damages
D(G)  ~ er₁*Y₁ + er₂*Y₂ - u*(G-2.8) + α*ct*gmax, # global externality.
D(z)  ~ η*θ(G-Tg)*(1-z), #initiate decarbonization
D(e₁) ~ -z*re₁*e₁, # decarbonization in region 1
D(e₂) ~ -z*re₂*e₂, # decarbonization in region 2
]



end # end of module