module traffic_model

using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# Hard code for now
NumPatches = 2
NumCors = 1
export NumPatches, NumCors

# Create variables for population model with travel demand
@variables np(t)[1:NumPatches] nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors] γ(t)[1:NumPatches]
export np, nc, γ

#=
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
            - m: carrying capacity (supremum of the values of the function)

    Then, to address (2), we assume that as current congestion levels rise, demand to 
    leave a given patch will decrease (from its baseline time-based level.)

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

# Define variables and parameters for internal clock
@variables u(t) v(t)
export u, v

@parameters period = 1440
export period

# Create internal clock (dynamical) equations
# I don't think the dot notation should be necessary? But trying to debug...
clock_eqs = [
    D(u) ~ u .* (1 .- u .^ 2 .- v .^ 2) .- (2 .* pi ./ period) .* (v),
    D(v) ~ v .* (1 .- u .^ 2 - v .^ 2) .+ (2 .* pi ./ period) .* (u),
]

# Parameters for time-based demand function (using logistic function)
@parameters dsharp = 100 dur = 0.97 shift1 = 0 shift2 = 0
export dsharp, dur, shift1, shift2
# ? Should I give these names different from what is in the function?
# x0 seems to be duration? Unclear why...

# Function to shift the periodic function v(t) left or right
# ? Do I need to pass period as an argument?
v_shifted(s) = v * cos(2 * π * s / period) - u * sin(2 * π * s / period)

# Logistic function
f(x, r, x0, m) = m ./ (1 .+ exp.(-r .* (x .- x0)))
# where r = dsharp, x0 =  dur, m = 1

# Sigmoid function
g(x, a, b, c) = a .* exp.(-log(2) .* (x ./ b) .^ c)

# Parameters for congestion-based demand function (using sigmoid function)
# Each patch has a background demand level, bgd
@parameters demhf1 = 3000 demhf2 = 500 bgd1 = 0 bgd2 = 0
export demhf1, demhf2, bgd1, bgd2

# Create complete demand equations (product of time-based and congestion-based functions)
# Do these not require that I pass some arguments?
# Also, will my dimensions be correct for these vectors?
# And... is it fine that np is just a number, not a density?
demand_eqs = [
    #γ ~ (1 .- g(np, 1, demhf1, 4)) .* f(v_shifted(shift1), dsharp, dur, 1) .+ bgd1,
    γ[1] ~ (1 .- g(np[1], 1, demhf1, 4)) .* f(v_shifted(shift1), dsharp, dur, 1) .+ bgd1,
    γ[2] ~ (1 .- g(np[2], 1, demhf2, 4)) .* f(-v_shifted(shift2), dsharp, dur, 1) .+ bgd2
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
@parameters ψ L
export ψ, L

# Parameters for velocity-density relation on the corridor
@parameters vff[1:NumPatches, 1:NumPatches, 1:NumCors] nc_half_ff[1:NumPatches, 1:NumPatches, 1:NumCors] vsharp[1:NumPatches, 1:NumPatches, 1:NumCors]
export vff, nc_half_ff, vsharp

# Function to compute average velocity on the corridor
vel_cor(nc, ψ, L, vff, nc_half_ff, vsharp) = g.(nc ./ (ψ * L), vff, nc_half_ff, vsharp)

# Parameters for velocity-density relation on the on-ramp
@parameters onff on_half onsharp
export onff, on_half, onsharp

# When do I need . notation and when do I not?
# Also, do I need to pass the variables every time? The parameters?
vel_on_ramp(nc, ψ, L, onff, on_half, onsharp) = g.(nc ./ (ψ * L), onff, on_half, onsharp)

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
@parameters entry_on exit_on
export entry_on, exit_on

# Rescale time in terms of a day that is 1440 minutes (24 hours)
time_rescale = 1440 / period
export time_rescale

# Which arguments do I need to pass vs. not pass? Do I pass vel_on_ramp (the function
# itself) as an argument of EntryFlux?
# Also, we seem to be missing some factors of ψ and L in these...
EntryFlux(nc, γ, ψ, L, onff, on_half, onsharp) = vel_on_ramp(nc, ψ, L, onff, on_half, onsharp) .* γ .* time_rescale .* entry_on
ExitFlux(nc, L, vff, nc_half_ff, vsharp) = (nc ./ L) .* vel_cor(nc, ψ, L, vff, nc_half_ff, vsharp) .* time_rescale .* exit_on
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
    For D.(p), we summed over the columns of the EntryFlux matrix and over the columns 
    of the ExitFlux matrix, then took their difference. If you don't do this (summing),
    you will get an error because you will be attempting to save an object of shape 
    (2,2,1) into an object of shape (2,). (Assuming NumPatches=2)
=#

# Define dynamical equations
eqs = [
    ϕ_in ~ collect(EntryFlux(nc, γ, ψ, L, onff, on_half, onsharp)),
    ϕ_out ~ collect(ExitFlux(nc, L, vff, nc_half_ff, vsharp)),
    D.(np) ~ [sum(ExitFlux(nc, L, vff, nc_half_ff, vsharp)[:, i, :]) for i in 1:NumPatches] .- [sum(EntryFlux(nc, γ, ψ, L, onff, on_half, onsharp)[i, :, :]) for i in 1:NumPatches],
    D.(nc) ~ collect(-ExitFlux(nc, L, vff, nc_half_ff, vsharp) .+ EntryFlux(nc, γ, ψ, L, onff, on_half, onsharp)),
    clock_eqs...,
    demand_eqs...   # <-- include the triple dots to splice the equations into the list
]

println("extra functions...")

function construct_ode_system()
    """This function is the main function to construct our model in ModelingToolkit."""

    @named ode = ODESystem(eqs, t)

    return ode

end # end of function

function construct_ode_problem(ode, tspan, u0, p)

    prob = ODEProblem(complete(ode), u0, tspan, p)

    return prob

end # end of function

end #end of module

#############
# OLD TESTS #
#############

#=
# Try it out
try_np = [0.9, 0.1]
try_demhf1 = 1
try_dsharp = 1
try_dur = 1
try_shift1 = 0
println((1 .- g(try_np, 1, try_demhf1, 4)))
#println(f(v_shifted(try_shift1), try_dsharp, try_dur, 1))


# Try it out
try_np = [1, 9]
try_nc = zeros(NumPatches, NumPatches, NumCors)
try_nc[1, 2, 1] = 8
try_nc[2, 1, 1] = 4
display(try_nc)
try_ψ = 100
try_L = 30
try_demhf1 = 3000
try_dsharp = 100
try_dur = 0.97
try_shift1 = 0
try_γ = [0.1, 0.9]
try_onff = 5000
try_onhj = 0.5
try_onsharp = 1
#println((1 .- g(try_np, 1, try_demhf1, 4)))
println("TRY")
try_entryflux = EntryFlux(try_nc, try_γ, try_ψ, try_L, try_onff, try_onhj, try_onsharp)
display(try_entryflux)
display(try_entryflux[1, :, :])
display(try_entryflux[:, 1, :])
println("and the sum...")
display([sum(try_entryflux[i, :, :]) for i in 1:NumPatches])
println("DONE TRYING")
=#