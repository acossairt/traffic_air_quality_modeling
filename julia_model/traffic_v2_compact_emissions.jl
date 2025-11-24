using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

# minimal settings (hardcode for now)
NumPatches = 2
NumCors = 1

#variables, patch, corridor populations, trip demand, diurnal clock
#average velocity and emission flux per corridor and total emissions
@variables np(t)[1:NumPatches] 
@variables nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors] 
@variables γ(t)[1:NumPatches, 1:NumPatches]
@variables v̄(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
@variables ec(t)[1:NumPatches,1:NumPatches, 1:NumCors]
@variables et(t)
@variables u(t) v(t)

# parameters
@parameters dsharp = 100 dur = 0.97 shift1 = 1 shift2 = -2
@parameters demhf12 = 3000 demhf21 = 500 bgd1 = 0 bgd2 = 0
@parameters ψ =100 L = 30
@parameters vff[1:NumPatches,1:NumPatches,1:NumCors] 
@parameters nc_half_ff[1:NumPatches,1:NumPatches,1:NumCors] 
@parameters vsharp[1:NumPatches,1:NumPatches,1:NumCors]
@parameters onff = 5000 on_half = 0.5 onsharp = 1
@parameters period = 24

# trip demand and corridor velocity helper functions
f(x, r, x0, m) = m ./ (1 .+ exp.(-r .* (x .- x0)))
g(x, a, b, c) = a .* exp.(-log(2) .* (x ./ b) .^ c)

# time shift equation for diurnal clock
v_shifted(shift) = v * cos(2 * π * shift / period) - u * sin(2 * π * shift / period)

# diurnal clock equations (symbolic expressions)
clock_eqs = [
    D(u) ~ u .* (1 .- u .^ 2 .- v .^ 2) .- (2 .* pi ./ period) .* (v),
    D(v) ~ v .* (1 .- u .^ 2 - v .^ 2) .+ (2 .* pi ./ period) .* (u),
]

# trip demand eqs (symbolic expressions)
demand_eqs = [
    γ[1,1] ~ 0,
    γ[1,2] ~ (1 .- g(np[1], 1, demhf12, 4)) .* f(v_shifted(shift1), dsharp, dur, 1) .+ bgd1,
    γ[2,1] ~ (1 .- g(np[2], 1, demhf21, 4)) .* f(-v_shifted(shift2), dsharp, dur, 1) .+ bgd2,
    γ[2,2] ~ 0
]

# velocity expressions: don't collect; keep as ArrayOps
vel_cor  = g.(nc ./ (ψ * L), vff, nc_half_ff, vsharp)
ramp_flux = g.(nc ./ (ψ * L), onff, on_half, onsharp)

# Corridor Entry/Exit flux expressions (symbolic)
CorEntryFlux = ramp_flux .* γ     # shape: (NumPatches,NumPatches,NumCors)
CorExitFlux  = (nc ./ L) .* vel_cor

# sum over destination patches originating out of each patch -> length NumPatches
PatchExitFlux  = [sum(collect(CorEntryFlux[i, :, :])) for i in 1:NumPatches] 

# sum over origin patches arriving in each patch -> length NumPatches
PatchEntryFlux = [sum(collect(CorExitFlux[:, i, :])) for i in 1:NumPatches] 

#Calculate emission fluxes per corridor and total
CorEmissionFlux = 0.03 * nc .* vel_cor
#Sum over all corridors
TotEmissionFlux = sum(collect(CorEmissionFlux[:,:,:]))

# now build dynamics using expressions directly and include some 
# intermediate quantities for visualization
eqs = [
D.(np) ~ PatchEntryFlux .- PatchExitFlux,
D.(nc) ~ CorEntryFlux .- CorExitFlux,
demand_eqs...,
clock_eqs...,
v̄ ~ vel_cor,
ec ~ CorEmissionFlux,
D(et) ~ TotEmissionFlux
]

# Assemble the ODE system
@named ode = ODESystem(eqs, t)

# Assemble parameters and initial conditions
p = merge(
    ModelingToolkit.get_defaults(ode),
    Dict(vff[i,j,k] => (i==j ? 0.0 : 120.0) for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors),
    Dict(nc_half_ff[i,j,k] => (i==j ? 0.0 : 0.3) for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors),
    Dict(vsharp[i,j,k] => (i==j ? 0.0 : 1.0) for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors)
)

u0 = [
    np[1] => 10000,
    np[2] => 0,
    nc => zeros(NumPatches, NumPatches, NumCors),
    u => 1.0,
    v => 0.0,
    et => 0.0,
]

tspan = (0.0, 48.0)

# Create the ODE problem
prob = ODEProblem(mtkcompile(ode), merge(Dict(u0),p), tspan)

# Solve the ODE problem
sol = solve(prob, Tsit5(), saveat=0.1)

# Plot results
plot1=plot(sol, idxs=[np[1], np[2], nc[1,2,1], nc[2,1,1]], xlabel="Time (hours)", ylabel="Population", title="Patch Populations Over Time")
plot2=plot(sol, idxs=[nc[1,2,1],10 .* v̄[1,2,1],ec[1,2,1],nc[2,1,1], 10 .* v̄[2,1,1], ec[2,1,1], et])
display(plot1)
display(plot2)