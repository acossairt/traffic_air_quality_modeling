module traffic_model

using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using Interpolations
using IfElse

# Hard code for now
Np = 2
Nc = 1
export Np, Nc

# Create parameters for speed-density curve
@parameters v_f[1:Np, 1:Np, 1:Nc] λ kc_crit = 0.5
export v_f, λ, kc_crit

# Create parameters for "demand-to-leave" function
@parameters r x_0 L shift
export r, x_0, L, shift

# Create parameters for population model
@parameters α[1:Np] kc_half_jam[1:Np, 1:Np, 1:Nc]
export α, kc_half_jam

# Create parameters for internal "clock" model
@parameters period
export period

# Create parameters for "turning off" entry and exit fluxes (test-purposes)
# Value of 1 keeps flux on, value of 0 turns flux off
@parameters exit_on = 1 entry_on = 1
export exit_on, entry_on

# Create variables for population model with travel demand
@variables kp(t)[1:Np] kc(t)[1:Np, 1:Np, 1:Nc] γ(t)[1:Np]
export kp, kc, γ

# Create variables for internal "clock" model
@variables u(t) v(t)
export u, v

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

# Define speed-density curve -- lots of options to consider!

# Greenshields (1935)
greenshield(kc, v_f, kc_half_jam) = v_f .* (1 .- kc ./ (2 * kc_half_jam))

# Drake (1967)
drake(kc, v_f, a, kc_half_jam) = v_f .* exp.((-1 / 2) .* (kc ./ (2 * kc_half_jam)) .^ 2) .* a # silly way to cheat parameter a

# Drake flexible
drake_2(kc, v_f, a, kc_half_jam) = v_f .* exp.((-1 / (2 * a)) .* (kc ./ (2 * kc_half_jam)) .^ 2)

# Smulders (1990)
d(v_f, a, k_crit, k_jam) = (v_f - a * k_crit) / ((1 / k_crit) - (1 / k_jam))
smulders(kc, v_f, a, kc_crit, kc_half_jam) = kc ≤ kc_crit ? v_f .- a .* kc : d(v_f, a, k_crit, 2 * kc_half_jam) .* (1 ./ kc .- 1 / (2 * kc_half_jam)) # assume kc > 0 and < k_jam always

# Daganzo (1994)
b(v_f, w, k_crit, k_jam) = (v_f + w * k_crit) / (w * k_jam)
daganzo(kc, v_f, w, kc_crit, kc_half_jam) = kc ≤ kc_crit ? v_f : -w + (b(v_f, w, k_crit, 2 * kc_half_jam) * w * 2 * kc_half_jam) ./ kc # assume kc > 0 and < k_jam always

# Original version
β(kc_half_jam) = 1 ./ (kc_half_jam .* 2)
old_exp(kc, kc_half_jam) = exp.(-β.(kc_half_jam) .* kc)

# Custom version (arctan function)
η(kc, a, kc_half_jam) = (pi / 2 .- atan.(a .* (kc .- kc_half_jam))) ./ pi # arctan function with inflection point at kc_half, asymptotically approaches y=0
ϕ(kc, a, kc_half_jam) = η(kc, a, kc_half_jam) ./ η(0.0, a, kc_half_jam)        # rescale so y-range is (0,1)
custom(kc, v_f, a, kc_half_jam) = v_f .* ϕ(kc, a, kc_half_jam)

# Current (favorite) version
# Calculate parameter b so that v(k_hj) = 1/2 v_f # need a better name for b
calc_b(k_half_jam, λ) = k_half_jam .* (1 ./ log(2)) .^ (1 ./ λ)
favorite(kc, v_f, λ, b) = v_f .* exp.(-1 .* (kc ./ b) .^ λ)

#=
# Newest custom version (arctan function)
η(kc, a, kc_half_jam) = (atan.(a .* (kc .- kc_half_jam))) ./ pi # arctan function with inflection point at kc_half, symmetric around y=0
ϕ(kc, a, kc_half_jam) = η(kc, a, kc_half_jam) ./ (2 * η(0.0, a, kc_half_jam))       # rescale so y-range is (-0.5,0.5)
custom(kc, v_f, a, kc_half_jam) = v_f .* (ϕ(kc, a, kc_half_jam) .+ 1 / 2)           # shift and rescale so y-range is (0, v_f)
=#

#=
function avg_speed(kc, v_f, a, kc_half_jam, kc_crit, version) # ::String
    if version == "greenshield"
        return greenshield(kc, v_f, kc_half_jam)
    elseif version == "drake"
        return drake(kc, v_f, kc_half_jam)
    elseif version == "smulders"
        return smulders(kc, v_f, a, kc_crit, kc_half_jam)
    elseif version == "daganzo"
        return daganzo(kc, v_f, a, kc_crit, kc_half_jam) # a = w
    elseif version == "arctan"
        return custom(kc, v_f, a, kc_half_jam)
    elseif version == "old"
        return old_exp(kc, kc_half_jam)
    else # Default to arctan function
        return custom(kc, v_f, a, kc_half_jam)
    end
end
=#

#=
if version == "arctan"
    avg_speed(kc, v_f, a, kc_half_jam) = custom(kc, v_f, a, kc_half_jam)
elseif version == "drake"
    avg_speed(kc, v_f, a, kc_half_jam) = drake(kc, v_f, kc_half_jam)
else
    avg_speed(kc, v_f, a, kc_half_jam) = greenshield(kc, v_f, kc_half_jam)
end
=#
#avg_speed(kc, v_f, a, kc_half_jam, λ) = λ == 0 ? custom(kc, v_f, a, kc_half_jam) : drake(kc, v_f, kc_half_jam)
#felse(λ == 0.0, avg_speed(kc, v_f, a, kc_half_jam, λ)=custom(kc, v_f, a, kc_half_jam), avg_speed(kc, v_f, a, kc_half_jam, λ)=drake(kc, v_f, kc_half_jam))

#export new_avg_speed
# Desmos version for convenience
# (v_f / ((pi / 2 - arctan(a * (0 - h))) / pi)) * (pi / 2 - arctan(a * (x - h))) / pi

# Options! Choose your own adventure
#avg_speed(kc, λ, v_f, a, kc_half_jam) = λ .* new_avg_speed(kc, v_f, a, kc_half_jam) + (1 - λ) .* old_avg_speed(kc, kc_half_jam)

#avg_speed(kc, v_f, a, kc_half_jam) = drake_2(kc, v_f, a, kc_half_jam) #custom(kc, v_f, a, kc_half_jam) #drake(kc, v_f, a, kc_half_jam) .* (1 / a) # silly way to cheat parameter a
avg_speed(kc, v_f, λ, kc_half_jam) = favorite(kc, v_f, λ, calc_b(kc_half_jam, λ))
export avg_speed

#=
Define the corridor flux matrices:
    Notice the use of `.` notation before each arithmetic operation, but NOT before 
    `=`. This is because `.=` only makes sense if the LHS of the equation is an 
    existing vector (or matrix) which has exactly same shape as the RHS.
=#

println("define flux functions...")
EntryFlux(kp, kc, γ, α, kc_half_jam, v_f, λ, entry_on) = avg_speed(kc, v_f, λ, kc_half_jam) .* kp .* γ .* time_rescale .* entry_on
ExitFlux(kc, kc_half_jam, v_f, λ, exit_on) = avg_speed(kc, v_f, λ, kc_half_jam) .* kc .* time_rescale .* exit_on

#=
Define equations for model:
    Notice there are only two equations in the population model, regardless of 
    Np and Nc. The dynamical clock model has 3 or more equations (for u, v, and γ)
    Question: does it matter if u, v, and γ are evaluated before or after kp and kc?
    Notice: we only used "collect()" in the equation for D.(c).
    For D.(p), we summed over the columns of the EntryFlux matrix and over the columns 
    of the ExitFlux matrix, then took their difference. If you don't do this (summing),
    you will get an error because you will be attempting to save an object of shape 
    (2,2,1) into an object of shape (2,). (Assuming Np=2)
=#

println("define eqs...")
eqs = [
    D.(kp) ~ [sum(ExitFlux(kc, kc_half_jam, v_f, λ, exit_on)[:, i, :]) for i in 1:Np] .- [sum(EntryFlux(kp, kc, γ, α[1], kc_half_jam, v_f, λ, entry_on)[i, :, :]) for i in 1:Np],
    D.(kc) ~ collect(-ExitFlux(kc, kc_half_jam, v_f, λ, exit_on) + EntryFlux(kp, kc, γ, α[1], kc_half_jam, v_f, λ, entry_on)),
    D(u) ~ u * (1 - u^2 - v^2) - (2 * pi / period) * (v),
    D(v) ~ v * (1 - u^2 - v^2) + (2 * pi / period) * (u),
    γ_eqs...   # <-- include the triple dots to splice the equations into the list
]

#=
Why don't I export the following functions?
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
