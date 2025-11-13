using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# minimal settings (hardcode for now)
NumPatches = 2
NumCors = 1

#variables, patch, corridor populations, trip demand, diurnal clock
@variables np(t)[1:NumPatches] 
@variables nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors] 
@variables γ(t)[1:NumPatches]
@variables u(t) v(t)

# parameters
@parameters dsharp = 100 dur = 0.97 shift1 = 0 shift2 = 0
@parameters demhf1 = 3000 demhf2 = 500 bgd1 = 0 bgd2 = 0
@parameters ψ =100 L = 30
@parameters vff[1:NumPatches,1:NumPatches,1:NumCors] 
@parameters nc_half_ff[1:NumPatches,1:NumPatches,1:NumCors] 
@parameters vsharp[1:NumPatches,1:NumPatches,1:NumCors]
@parameters onff on_half onsharp
@parameters period = 24

# trip demand and corridor velocity helper functions
f(x, r, x0, m) = m ./ (1 .+ exp.(-r .* (x .- x0)))
g(x, a, b, c) = a .* exp.(-log(2) .* (x ./ b) .^ c)

# time shift equation for diurnal clock
v_shifted(shift) = v * cos(2 * π * shift / period) - u * sin(2 * π * shift / period)

#durnal clock equations (symbolic expressions)
clock_eqs = [
    D(u) ~ u .* (1 .- u .^ 2 .- v .^ 2) .- (2 .* pi ./ period) .* (v),
    D(v) ~ v .* (1 .- u .^ 2 - v .^ 2) .+ (2 .* pi ./ period) .* (u),
]

# trip demand eqs (symbolic expressions)
demand_eqs = [
    γ[1] ~ (1 .- g(np[1], 1, demhf1, 4)) .* f(v_shifted(shift1), dsharp, dur, 1) .+ bgd1,
    γ[2] ~ (1 .- g(np[2], 1, demhf2, 4)) .* f(-v_shifted(shift2), dsharp, dur, 1) .+ bgd2
]

# velocity expressions: don't collect; keep as ArrayOps
vel_cor  = g.(nc ./ (ψ * L), vff, nc_half_ff, vsharp)
ramp_flux = g.(nc ./ (ψ * L), onff, on_half, onsharp)

# Corridor Entry/Exit flux expressions (symbolic)
EntryFlux = ramp_flux .* γ     # shape: (NumPatches,NumPatches,NumCors)
ExitFlux  = (nc ./ L) .* vel_cor



# now build dynamics using expressions directly
eqs = [
#D.(np) ~ - [sum(EntryFlux[i, :, :]) for i in 1:NumPatches] .+ [sum(ExitFlux[:, i, :]) for i in 1:NumPatches],
D.(nc) ~ EntryFlux .- ExitFlux,
demand_eqs...,
clock_eqs...
]

@named ode = ODESystem(eqs, t)
# optionally inspect types now:
println("EntryFlux typeof = ", typeof(EntryFlux))
println("EntryFlux eltype = ", Base.eltype(EntryFlux))
