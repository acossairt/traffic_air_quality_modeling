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




# Code for estimating emissions#
# Make a U-shaped curve using data from the California paper
start = 5
my_step = 5
stop = 100
mph_to_kmh = 1.60934
speed_arr = collect(start:my_step:stop) * mph_to_kmh # convert mph to kmh
emissions_arr = [1200, 950, 700, 500, 425, 350, 325, 310, 309, 308, 308, 308, 309, 320, 330, 350, 375, 400, 450, 550] * mph_to_kmh
plot(speed_arr, emissions_arr)

# Interpolate: emissions as a function of speed
interp_fn = linear_interpolation(speed_arr, emissions_arr, extrapolation_bc=Line())

# Calculate emissions rates from a given array of speeds
function calc_emissions_from_speed(vehicle_pop_arr, my_speed_arr, interp_fn)
    #=
        Returns:
            - emissions: array (dim 1) 
                emission rates (g/km) for whole traffic volume (all vehicles) at each
                time step
        Arguments:
            - vehicle_pop_arr: array (dim 1) of vehicle population densities at 
              each time step
            - my_speed_arr: array (dim 1) of avg vehicle speeds at each time step
            - interp_fn: function (interpolated) relating speeds to emissions
    =#
    interpolated_emission_per_vehicle = interp_fn(my_speed_arr)
    emissions = interpolated_emission_per_vehicle .* my_speed_arr .* vehicle_pop_arr
    return emissions
end

# Calculate total emissions from emissions rates
function integrate_emissions(x, y, a, b)
    # Create an interpolation function over data points
    interp_func = LinearInterpolation(x, y, extrapolation_bc=Line())

    # Integrate the interpolated function from a to b
    result, error = quadgk(interp_func, a, b)

    return result, error
end

# Example implementation (variable names and sol variables will need to be updated)

# Calculate emission rates for Corridor 1
pop_C1 = sol[3, :] # needs updated
time = sol.t
C1_speeds = calc_space_mean_speed_alternative_greenshields.(v_f, pop_C1, C1_half)
C1_emissions = calc_emissions_from_speed(pop_C1, C1_speeds, interp_fn)
C1_flow = calc_flow(pop_C1, v_f, C1_jam)

# Same for Corridor 2
pop_C2 = sol[5, :] # needs updated
time = sol.t
C2_speeds = calc_space_mean_speed_alternative_greenshields.(v_f, pop_C2, C2_half)
C2_emissions = calc_emissions_from_speed(pop_C2, C2_speeds, interp_fn)
C2_flow = calc_flow(pop_C2, v_f, C2_jam)

# Total emissions
C1_total_emissions = integrate_emissions(time, C1_emissions, 0.0, 100.0)[1] # [1] for value, [2] for error
formatted_C1_em = @sprintf("%.3f", C1_total_emissions) # Nicer format for printing on plots

C2_total_emissions = integrate_emissions(time, C2_emissions, 0.0, 100.0)[1]
formatted_C2_em = @sprintf("%.3f", C2_total_emissions)

Total_emissions = C1_total_emissions + C2_total_emissions
formatted_T_em = @sprintf("%.3f", Total_emissions)

C1_fraction_emissions = 100 * C1_total_emissions / Total_emissions
formatted_C1_fraction = @sprintf("%.2f", C1_fraction_emissions)

C2_fraction_emissions = 100 * C2_total_emissions / Total_emissions
formatted_C2_fraction = @sprintf("%.2f", C2_fraction_emissions)