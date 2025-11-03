module traffic_model

using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using Interpolations
using IfElse

# Hard code for now
NumPatches = 2
NumCors = 1
export NumPatches, NumCors

# Create parameters for speed-density curve
@parameters v_f[1:NumPatches, 1:NumPatches, 1:NumCors] λ
export v_f, λ

# Create parameters for "demand-to-leave" function
@parameters r x_0 L shift
export r, x_0, L, shift

# Create parameters for road conditions and patch area
@parameters ψ Le
export ψ, Le

# Create parameters for population model
@parameters α[1:NumPatches] nc_half_jam[1:NumPatches, 1:NumPatches, 1:NumCors]
export α, nc_half_jam

# Create parameters for internal "clock" model
@parameters period
export period

# Create parameters for "turning off" entry and exit fluxes (test-purposes)
# Value of 1 keeps flux on, value of 0 turns flux off
@parameters exit_on = 1.0 entry_on = 1.0
export exit_on, entry_on

# Create variables for population model with travel demand
@variables np(t)[1:NumPatches] nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors] γ(t)[1:NumPatches]
export np, nc, γ

# Create variables for internal "clock" model
@variables u(t) v(t)
export u, v

# Create external variables to track fluxes
@variables ϕ_in(t)[1:NumPatches, 1:NumPatches, 1:NumCors] ϕ_out(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
export ϕ_in, ϕ_out

#=
Define demand function (to create "γ" for Flux functions):
    This is the logistic equation, which is smooth, continuous, and well-behaved
    If you pass a periodic function (like v) as the argument `x`, then f(x) will look
    like a (well-behaved) square wave, peaking at the same place as where the
    periodic function peaks.
=#
f(x, r, x_0, L) = L / (1 + exp(-r * (x - x_0)))

#=
Shift demand function left or right: 
    Using trig identity: sin(a-b) = sin(a)cos(b) - cos(a)sin(b)
    i.e. if `my_period=1440` (24 hours, in minutes), then my `v` will peak at `t=360` 
    (aka 6am). If I want morning demand to also peak at 6am, then leave `shift=0`. 
    However, if I instead want morning demand to peak at 9am, then I should pass 
    `shift=180` (3 hours * 60 minutes) to `v_shifted()`. The result will then
    automatically be passed to `f()`` to generate `γ`.
=#
v_shifted(shift) = v * cos(2 * π * shift / period) - u * sin(2 * π * shift / period)

# How to create options for periodic_logistic vs. static?
γ_eqs = [
    γ[1] ~ f(v_shifted(shift), r, x_0, L), # + 0.001 shift up just a hair so that demand is never 0
    γ[2] ~ f(-v_shifted(shift), r, x_0, L) # + 0.001
]

#=
Rescale time in terms of a day that is 1440 minutes
=#
time_rescale = 1440 / period
export time_rescale

# Define speed-density curve -- lots of other versions at bottom of file

#################################################
# NOTES FROM MARTY (attempt to make it cleaner) #
#################################################
#vel(car_density, vel_ff, b, sharpness) = vel_ff * exp(-(car_density / b)^sharpness)
#b(half_jam, sharpness) = half_jam .* (1./log(2)).^sharpness
#c_den(num_cars, ψ, Le) = num_cars ./ (ψ .* Le)
#avg_speed(nc, ψ, Le, vel_ff, sharpness, nc_half_jam) = vel(c_den(nc, ψ, Le), vel_ff, b(c_den(nc_half_jam, ψ, Le), sharpness), sharpness)

vel(car_density, vel_ff, b_factor, sharpness) = vel_ff .* exp.(-1 .* (car_density ./ b_factor) .^ sharpness)
b(half_jam, sharpness) = half_jam .* (1.0 / log(2)) .^ (1 ./ sharpness) # was somehow nicer when I had just sharpness instead of 1/sharpness?
num_to_den(num, ψ, Le) = num ./ (ψ .* Le)
avg_speed(nc, ψ, Le, v_f, sharpness, nc_half_jam) = vel(num_to_den(nc, ψ, Le), v_f, b(num_to_den(nc_half_jam, ψ, Le), sharpness), sharpness)
export avg_speed

# QUESTION
# When do I have to pass parameters and when can I just call upon them? For example:
#b = calc_b(num_to_den(nc_half_jam, ψ, Le), λ)

# Old version using densities (presumably)
#calc_b(kc_half_jam, λ) = (kc_half_jam) .* (1 ./ log(2)) .^ (1 ./ λ)
#avg_speed(kc, ψ, Le, v_f, λ, nc_half_jam) = v_f .* exp.(-1 .* ((kc) ./ calc_b(kc_half_jam, λ)) .^ λ)

# Exportable versions (pass all arguments)
out_calc_b(nc_half_jam, λ, ψ, Le) = (nc_half_jam ./ (ψ .* Le)) .* (1 ./ log(2)) .^ (1 ./ λ)
out_avg_speed(nc, ψ, Le, v_f, λ, nc_half_jam) = v_f .* exp.(-1 .* ((nc ./ (ψ .* Le)) ./ out_calc_b(nc_half_jam, λ, ψ, Le)) .^ λ)

#=
Define the corridor flux matrices:
    Notice the use of `.` notation before each arithmetic operation, but NOT before 
    `=`. This is because `.=` only makes sense if the LHS of the equation is an 
    existing vector (or matrix) which has exactly same shape as the RHS.
=#

println("define flux functions...")
# eventually need to change γ .* np to some function β(N_D) (defined in my notebook)
EntryFlux(nc, np, γ, nc_half_jam, v_f, λ, Le, ψ, entry_on) = ψ .* avg_speed(nc, ψ, Le, v_f, λ, nc_half_jam) .* γ .* np .* time_rescale .* entry_on
ExitFlux(nc, nc_half_jam, v_f, λ, Le, ψ, exit_on) = (nc ./ Le) .* avg_speed(nc, ψ, Le, v_f, λ, nc_half_jam) .* time_rescale .* exit_on
export EntryFlux, ExitFlux

#=
Define equations for model:
    Notice there are only two equations in the population model, regardless of 
    NumPatches and NumCors. The dynamical clock model has 3 or more equations (for u, v, and γ)
    Question: does it matter if u, v, and γ are evaluated before or after np and nc?
    Notice: we only used "collect()" in the equation for D.(c).
    For D.(p), we summed over the columns of the EntryFlux matrix and over the columns 
    of the ExitFlux matrix, then took their difference. If you don't do this (summing),
    you will get an error because you will be attempting to save an object of shape 
    (2,2,1) into an object of shape (2,). (Assuming NumPatches=2)
=#

println("define eqs...")
eqs = [
    ϕ_in ~ collect(EntryFlux(nc, np, γ, nc_half_jam, v_f, λ, Le, ψ, exit_on)),
    ϕ_out ~ collect(ExitFlux(nc, nc_half_jam, v_f, λ, Le, ψ, exit_on)),
    D.(np) ~ [sum(ExitFlux(nc, nc_half_jam, v_f, λ, Le, ψ, exit_on)[:, i, :]) for i in 1:NumPatches] .- [sum(EntryFlux(nc, np, γ, nc_half_jam, v_f, λ, Le, ψ, entry_on)[i, :, :]) for i in 1:NumPatches],
    D.(nc) ~ collect(-ExitFlux(nc, nc_half_jam, v_f, λ, Le, ψ, exit_on) + EntryFlux(nc, np, γ, nc_half_jam, v_f, λ, Le, ψ, entry_on)),
    D(u) ~ u * (1 - u^2 - v^2) - (2 * pi / period) * (v),
    D(v) ~ v * (1 - u^2 - v^2) + (2 * pi / period) * (u),
    γ_eqs...   # <-- include the triple dots to splice the equations into the list
]

#=
Why don't I have to export the following functions?
=#

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

######################################
# Old / previous avg_speed functions #
# Saving for posterity               #
######################################

#=
# Greenshields (1935)
greenshield(nc, v_f, nc_half_jam) = v_f .* (1 .- nc ./ (2 * nc_half_jam))

# Drake (1967)
drake(nc, v_f, a, nc_half_jam) = v_f .* exp.((-1 / 2) .* (nc ./ (2 * nc_half_jam)) .^ 2) .* a # silly way to cheat parameter a

# Drake flexible
drake_2(nc, v_f, a, nc_half_jam) = v_f .* exp.((-1 / (2 * a)) .* (nc ./ (2 * nc_half_jam)) .^ 2)

# Smulders (1990)
d(v_f, a, k_crit, k_jam) = (v_f - a * k_crit) / ((1 / k_crit) - (1 / k_jam))
smulders(nc, v_f, a, nc_crit, nc_half_jam) = nc ≤ nc_crit ? v_f .- a .* nc : d(v_f, a, k_crit, 2 * nc_half_jam) .* (1 ./ nc .- 1 / (2 * nc_half_jam)) # assume nc > 0 and < k_jam always

# Daganzo (1994)
b(v_f, w, k_crit, k_jam) = (v_f + w * k_crit) / (w * k_jam)
daganzo(nc, v_f, w, nc_crit, nc_half_jam) = nc ≤ nc_crit ? v_f : -w + (b(v_f, w, k_crit, 2 * nc_half_jam) * w * 2 * nc_half_jam) ./ nc # assume nc > 0 and < k_jam always

# Original version
β(nc_half_jam) = 1 ./ (nc_half_jam .* 2)
old_exp(nc, nc_half_jam) = exp.(-β.(nc_half_jam) .* nc)

# Custom version (arctan function)
η(nc, a, nc_half_jam) = (pi / 2 .- atan.(a .* (nc .- nc_half_jam))) ./ pi # arctan function with inflection point at nc_half, asymptotically approaches y=0
ϕ(nc, a, nc_half_jam) = η(nc, a, nc_half_jam) ./ η(0.0, a, nc_half_jam)        # rescale so y-range is (0,1)
custom(nc, v_f, a, nc_half_jam) = v_f .* ϕ(nc, a, nc_half_jam)

# Current (favorite) version
# Calculate parameter b so that v(k_hj) = 1/2 v_f # need a better name for b
calc_b(k_half_jam, λ) = k_half_jam .* (1 ./ log(2)) .^ (1 ./ λ)
favorite(nc, v_f, λ, b) = v_f .* exp.(-1 .* (nc ./ b) .^ λ)

# Greenshields (1935)
greenshield(nc, v_f, nc_half_jam) = v_f .* (1 .- nc ./ (2 * nc_half_jam))

# Drake (1967)
drake(nc, v_f, a, nc_half_jam) = v_f .* exp.((-1 / 2) .* (nc ./ (2 * nc_half_jam)) .^ 2) .* a # silly way to cheat parameter a

# Drake flexible
drake_2(nc, v_f, a, nc_half_jam) = v_f .* exp.((-1 / (2 * a)) .* (nc ./ (2 * nc_half_jam)) .^ 2)

# Smulders (1990)
d(v_f, a, k_crit, k_jam) = (v_f - a * k_crit) / ((1 / k_crit) - (1 / k_jam))
smulders(nc, v_f, a, nc_crit, nc_half_jam) = nc ≤ nc_crit ? v_f .- a .* nc : d(v_f, a, k_crit, 2 * nc_half_jam) .* (1 ./ nc .- 1 / (2 * nc_half_jam)) # assume nc > 0 and < k_jam always

# Daganzo (1994)
b(v_f, w, k_crit, k_jam) = (v_f + w * k_crit) / (w * k_jam)
daganzo(nc, v_f, w, nc_crit, nc_half_jam) = nc ≤ nc_crit ? v_f : -w + (b(v_f, w, k_crit, 2 * nc_half_jam) * w * 2 * nc_half_jam) ./ nc # assume nc > 0 and < k_jam always

# Original version
β(nc_half_jam) = 1 ./ (nc_half_jam .* 2)
old_exp(nc, nc_half_jam) = exp.(-β.(nc_half_jam) .* nc)

# Custom version (arctan function)
η(nc, a, nc_half_jam) = (pi / 2 .- atan.(a .* (nc .- nc_half_jam))) ./ pi # arctan function with inflection point at nc_half, asymptotically approaches y=0
ϕ(nc, a, nc_half_jam) = η(nc, a, nc_half_jam) ./ η(0.0, a, nc_half_jam)        # rescale so y-range is (0,1)
custom(nc, v_f, a, nc_half_jam) = v_f .* ϕ(nc, a, nc_half_jam)

# Current (favorite) version
# Calculate parameter b so that v(k_hj) = 1/2 v_f # need a better name for b
calc_b(k_half_jam, λ) = k_half_jam .* (1 ./ log(2)) .^ (1 ./ λ)
favorite(nc, v_f, λ, b) = v_f .* exp.(-1 .* (nc ./ b) .^ λ)
=#