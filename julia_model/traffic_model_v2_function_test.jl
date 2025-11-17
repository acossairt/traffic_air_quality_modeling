

using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# Hard code for now
NumPatches = 2
NumCors = 1
export NumPatches, NumCors

# Create patch density, corridor density variables for population model with travel demand
@variables np(t)[1:NumPatches]
@variables nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
@variables γ(t)[1:NumPatches]
export np, nc, γ


#=Functions for trip demand and corridor velocities

    Modeling time-based demand for travel
    -------------------------------------
    We expect travel demand (to leave a given patch) varies with 1) time of day
    and 2) current congestion levels.
    To address (1), the following functions define an internal clock and then 
    determine time-based demand.

    - v(t) and u(t)
        - Periodic functions sin(t) and cos(t), respectively.
        - Together, they create a clock for one day. (See clock_eqs)

    - v_shifted(shift, period)
        - Shifts the periodic function v(t) left or right.
        - Trig identity sin(a-b) = sin(a)cos(b) - cos(a)sin(b)
        - Use case: if `my_period=1440` (24 hours, in minutes), then `v` will peak at 
        `t=360` (aka 6am). If I want morning demand to also peak at 6am, then leave 
        `shift=0`. If I instead want morning demand to peak at 9am, then pass `shift=180` 
        (3 hours * 60 minutes) to `v_shifted()`

    - f(x,r,x0,m)
        - Logistic equation (which is smooth, continuous, and well-behaved). Commonly
        used to model discrete decisions (in our case, vehicles are deciding "do I stay 
        or do I go?")
        - If you pass a periodic function (like v(t)) as the argument `x`, then f(x)
        will look like a well-behaved wave, peaking simultaneously with v(t).
        - Parameters
            - r: logistic growth rate (controls steepness of the curve)
            - x0: x-value of the function's midpoint
            - m: carrying capacity ontinuous, and well-behaved). Commonl(supremum of the values of the function)

    Then, to address (2), we assume that as current congestion levels rise, demand to 
    leave a given patch will decrease (from its baseline time-based level.)

=#

f(x, r, x0, m) = m ./ (1 .+ exp.(-r .* (x .- x0)))
# where r = dsharp, x0 =  dur, m = 1, and x is the time of day function v_shifted(shift)    

# Then need a diurnal internal clock (dynamical) equations

@variables u(t) v(t)
export u, v
@parameters period = 24
export period

clock_eqs = [
    D(u) ~ u * (1 - u^2 - v^2) - (2 * pi / period) * (v),
    D(v) ~ v * (1 - u^2 - v^2) + (2 * pi / period) * (u),
]

# Then need function to create shifted version of v(t) to input into f().

v_shifted(shift) = v * cos(2 * π * shift / period) - u * sin(2 * π * shift / period)

#= Modeling congestion-based corridor velocities
    --------------------------------------------
    - g(x,a,b,c)
        - Sigmoid function, also a generic form of the Drake equation relating
        traffic density and average velocities (as congestion rises, avg velocity falls)
        - Note that: 1 - g ~ f. The transformation (1-g(x)) gives a function (demand) which
        ... rises with congestion? That can't be right...
        - Parameters
            - a: supremum of the values of the function (e.g. maximum velocity)
            - b: fixes the half-point (the value of x where g(x) = (1/2)*a)
            - c: controls sharpness of curve drop-off

    Finally, the complete demand is a product of time-based demand (f) and 
    congestion-based demand (1 - g), plus a background demand constant (bgd).
=#

g(x, a, b, c) = a .* exp.(-log(2) .* (x ./ b) .^ c)


# Parameters for time-based demand function (using logistic function)
@parameters dsharp[1:NumPatches] dur[1:NumPatches] shift[1:NumPatches]
export dsharp, dur, shift
# ? Should I give these names different from what is in the function?
# x0 seems to be duration? Unclear why...


# Parameters for congestion-based demand function (using sigmoid function)
# Each patch has a background demand level, bgd

@parameters demhf[1:NumPatches] bgd[1:NumPatches]
export demhf, bgd

# Create complete demand equations (product of time-based and congestion-based functions)
# Do these not require that I pass some arguments?
# Also, will my dimensions be correct for these vectors?
# And... is it fine that np is just a number, not a density?

demand_eqs = [
    γ[1] .~ (1 .- g(np[1], 1, demhf[1], 4)) .* f(v_shifted(shift[1]), dsharp[1], dur[1], 1) .+ bgd[1],
    γ[2] .~ (1 .- g(np[2], 1, demhf[2], 4)) .* f(v_shifted(shift[2]), dsharp[2], dur[2], 1) .+ bgd[1]
]

#=
    Modeling travel velocities
    --------------------------
    The sigmoid function above (g) is the  generic form of the Drake equation commonly
    used in transportation engineering literature. We can use this function to relate
    traffic density (vehicles / area) to average velocities.

    First compute average travel velocity on the corridor
    - g(nc/(ψ*L), vff, nc_half_ff, vsharp) | Units: (vehicles / time)
        - nc / (ψ * L) | (vehicles / area): current traffic density in corridor.
            - Note: (ψ*L = area of unit of road)
        - vff | (velocity): free-flow velocity of corridor
        - nc_half_ff | (vehicles / area): traffic density value at which average 
          velocities on the corridor drop to 1/2 free-flow velocity
        - vsharp | (unitless): parameter controlling sharpness of declining velocities

    Then compute entry-velocity (e.g. the rate at which vehicles can enter a corridor)
    - g(nc / (ψ*L), onff, on_half, onsharp) | Units: (vehicles / time)
        - nc / (ψ*L) | (vehicles / area): current traffic density in corridor. Note
        (ψ*L = area of unit of road)
        - onff | (vehicles / time): free-flow traffic flux of on-ramp
        - on_half | (vehicles / area): traffic density on the on-ramp at which entry 
          flux rate drops to 1/2 of free-flow on-ramp flux-rate
        - onsharp | (unitless): parameter controlling sharpness of declining velocities
=#

# Parameters for road geometry and velocity-density relation on the corridor
@parameters ψ[1:NumCors] L[1:NumCors]
export ψ, L

# Parameters for velocity-density relation on the corridor
@parameters vff[1:NumPatches, 1:NumPatches, 1:NumCors]
@parameters nc_half_ff[1:NumPatches, 1:NumPatches, 1:NumCors]
@parameters vsharp[1:NumPatches, 1:NumPatches, 1:NumCors]
export vff, nc_half_ff, vsharp

# Function to compute average velocity on the corridor
vel_cor = g.(nc ./ (ψ .* L), vff, nc_half_ff, vsharp)
export vel_cor

# Parameters for velocity-density relation on the on-ramp
@parameters onff on_half onsharp
export onff, on_half, onsharp

# When do I need . notation and when do I not?
# Also, do I need to pass the variables every time? The parameters?
ramp_flux = g.(nc ./ (ψ .* L), onff, on_half, onsharp)
export ramp_flux

#=
    Modeling entry and exit fluxes
    ------------------------------
    - EntryFlux()
        - Demand γ gives number of vehicles available to enter the corridor
        - Vel_on_ramp gives rate at which available vehicles do enter the corridor

    - ExitFlux()
        - nc / L gives number of vehicles on the corridor, available to exit
        # ? How to explain this part, these units?
        - vel_cor gives rate at which available vehicles can exit
=#

# Parameters for testing entry and exit flux functions
#@parameters entry_on exit_on
#export entry_on, exit_on

# Rescale time in terms of aperiod=24 day that is 1440 minutes (24 hours)
#time_rescale = 1440 / period
#export time_rescale

# Which arguments do I need to pass vs. not pass? Do I pass vel_on_ramp (the function
# itself) as an argument of EntryFlux?
# Also, we seem to be missing some factors of ψ and L in these...

EntryFlux = ramp_flux .* γ
ExitFlux = (nc ./ L) .* vel_cor
export EntryFlux, ExitFlux
# Create external variables to track fluxes
@variables ϕ_in(t)[1:NumPatches, 1:NumPatches, 1:NumCors] ϕ_out(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
export ϕ_in, ϕ_out

#=
Define dynamical equations for the full system:
    Notice there are only two equations in the population model, regardless of 
    NumPatches and NumCors.
    Because our state variables `np` and `nc` are vectors, we must sometimes use the
    "collect()" function to appropriately sum. 
    Notice: we only used "collect()" in the equation for D.(c).
    For D.(np), we summed over the columns of the EntryFlux matrix and over the columns 
    of the ExitFlux matrix, then took their difference. If you don't do this (summing),
    you will get an error because you will be attempting to save an object of shape 
    (2,2,1) into an object of shape (2,). (Assuming NumPatches=2)
=#

println("\nTesting variables")
ϕ_out_test = collect(ExitFlux)
display(ϕ_out_test)

# Define dynamical equations
println("\nDefining dynamical equations...")

ϕ_in = collect(EntryFlux)
ϕ_out = collect(ExitFlux)

eqs = [
    D.(np) ~ -[sum(ϕ_in[i, :, :]) for i in 1:NumPatches] .+ [sum(ϕ_out[:, i, :]) for i in 1:NumPatches],
    D.(nc) ~ collect(-ExitFlux(nc, L, vff, nc_half_ff, vsharp) .+ EntryFlux(nc, γ, ψ, L, onff, on_half, onsharp)), #ϕ_in .- ϕ_out,
    clock_eqs...,
    demand_eqs...   # <-- include the triple dots to splice the equations into the list
]

#=
# Claude solution 1
eqs = [
    vec(ϕ_in .~ EntryFlux)...,  # Flatten to individual scalar equations
    vec(ϕ_out .~ ExitFlux)...,  # Flatten to individual scalar equations
    D.(np) .~ -[sum(ϕ_in[i, :, :]) for i in 1:NumPatches] .+ [sum(ϕ_out[:, i, :]) for i in 1:NumPatches],
    vec(D.(nc) .~ ϕ_in .- ϕ_out)...,  # Flatten to individual scalar equations
    clock_eqs...,
    demand_eqs...
]
=#

#=
# Claude solution 2
eqs = [
    [ϕ_in[i, j, k] .~ EntryFlux[i, j, k] for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors]...,
    [ϕ_out[i, j, k] .~ ExitFlux[i, j, k] for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors]...,
    [D(np[i]) .~ -sum(ϕ_in[i, :, :]) .+ sum(ϕ_out[:, i, :]) for i in 1:NumPatches]..., # wrong?
    [D(nc[i, j, k]) .~ ϕ_in[i, j, k] .- ϕ_out[i, j, k] for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors]...,
    clock_eqs...,
    demand_eqs...
]
=#

# Create and return ODE problem and solve it
#prob = ODEProblem(complete(ode), u0, tspan, p)
#sol = solve(prob, Tsit5(); saveat=1.0)

##################################################
# Test for parameter values used in XPPAUT study #
##################################################

NumPatches = 2
NumCors = 1

println("NumPatches: $NumPatches, NumCors: $NumCors")

# Set timescale
my_period = 24           # 24 hours
t_end = my_period        # Change if you only want to examine part of a day (< my_period), or multiple days (> my_period)
tspan = (0.0, t_end);    #time span in minutes.

# Set initial conditions
np_init = [10000; 500]
nc_init = zeros(NumPatches, NumPatches, NumCors)
u_init = 1
v_init = 0
γ_init = [0.6; 0.6] # I think this is basically trprp?

u0 = [nc => nc_init, np => np_init, γ => γ_init, u => u_init, v => v_init];

println("Setting up parameters...")

# Parameters for calculating average space-mean speed, and therefore fluxes
my_dsharp = 200 * ones(NumPatches)
my_dur = 0.97 * ones(NumPatches)
my_shift = 0 * ones(NumPatches)
my_demhf = 0.3 * ones(NumPatches)
my_bgd = 0 * ones(NumPatches)
my_ψ = 100 * ones(NumCors)
my_L = 30 * ones(NumCors)
my_vff = 120 * ones(NumPatches, NumPatches, NumCors)
my_vff[[CartesianIndex(i, i, k) for i in 1:NumPatches, k in 1:NumCors]] .= 0 # no movement in loops
my_vsharp = 1 * ones(NumPatches, NumPatches, NumCors)
my_nc_half_ff = 0.3 * ones(NumPatches, NumPatches, NumCors)
my_nc_half_ff[[CartesianIndex(i, i, k) for i in 1:NumPatches, k in 1:NumCors]] .= 1e9 # set diagonal values to almost Inf
my_onff = 5000
my_on_half = 0.5
my_onsharp = 1

println("Free-flow velocities (my_vff):")
display(my_vff)

println("Half ff densities:")
display(my_nc_half_ff)

# Make parameters list
p = [
    dsharp => my_dsharp,
    dur => my_dur,
    shift => my_shift,
    demhf => my_demhf,
    bgd => my_bgd,
    ψ => my_ψ,
    L => my_L,
    vff => my_vff,
    vsharp => my_vsharp,
    nc_half_ff => my_nc_half_ff,
    onff => my_onff,
    on_half => my_on_half,
    onsharp => my_onsharp
];

# Create ODE system
println("Creating ODE system...")
@named ode = ODESystem(eqs, t)

# Simplify structure and set up model
println("Simplifying structure...") # because the system has vector unknowns, we need to simplify it
simple_sys = structural_simplify(ode)

println("Creating problem...")
#prob = ODEProblem(complete(ode), u0, tspan, p)
prob = ODEProblem(simple_sys, u0, tspan, p)

#println("Model setup done!")

#sol = solve(prob, Tsit5())