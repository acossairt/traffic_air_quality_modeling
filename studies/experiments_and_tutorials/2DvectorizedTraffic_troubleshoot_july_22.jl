# Fast call for helper functions to run model
aurora_path_to_help_funcs = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/help_funcs.jl"
marty_path_to_help_funcs = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/help_funcs.jl"
include(marty_path_to_help_funcs) # CHANGE FOR YOURSELF
using .HelpFuncs
using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra, Printf

#=
Set initial conditions, parameters, and other arguments to pass to model solver
=#

println("Setting up model parameters...")

# Set number of patches and corridors

Np = 2 # Number of patches
Nc = 1 # Number of corridors

@printf "Number of patches: %d\n" Np
@printf "Number of corridors: %d\n" Nc

# Set initial conditions
kp_init = [1.0; 0.0]
kc_init = zeros(Np, Np, Nc)
u_init = 1
v_init = 0

println("Initial patch densities (kp_init):")
display(kp_init)
println("Initial corridor densities (kc_init):")
display(kc_init)


# Set timescale
my_period = 1440    # 24 hours, in minutes
my_t_end = my_period     # Change if you only want to examine part of a day (< my_period), or multiple days (> my_period)
scale = 1440 / my_period # if using a period other than 1440, this will rescale the other parameters

# Set parameters
my_α = [1]   # Tolerance for congestion (for now, assume same for patch 1 and patch 2)
kc_jam_base = 10
my_kc_jam = kc_jam_base * ones(Np, Np, Nc)  # Inverse road capacity (for now, assume same for all corridors)
my_kc_jam[[CartesianIndex(i, i, k) for i in 1:Np, k in 1:Nc]] .= 1e9 # set diagonal values to almost Inf
display(my_kc_jam)

# TEST: set C_init to C_jam
#kc_init[1,2,1] = 1 ./ my_kc_jam[1,2,1] + 0.05 # Should be 1/20
display(kc_init)

# Arguments for demand function (which is f() = L / (1 + exp(-r * (x-x0))) )
demand_function = "periodic_logistic" # Default option is "static"
my_L = 0.6
my_r = 100
my_x0 = 0.99
my_shift = 0 * 60    # in units of minutes, set to 0 for no shift

# Parameters for calculating space-mean speed
my_v_f = 90
my_a = 20

#=
Build model, solve problem, plot results, and save
=#

println("Building model and solving...  ")

model, prob, sol, plt = HelpFuncs.build_symbolic_model(Np=Np, Nc=Nc,
    my_kp=kp_init, my_kc=kc_init, my_u=u_init, my_v=v_init, my_α=my_α, my_kc_jam=my_kc_jam,
    γ_version=demand_function, my_period=my_period, t_end=my_t_end, my_v_f=my_v_f, my_a=my_a, my_L=my_L,
    my_r=my_r, my_x0=my_x0, my_shift=my_shift)

# Suffix to specify test
suffix = "_new_speed_density_flux" #"_original" #"_shut_off_entry_flux" #"_new_speed_density" #"_shut_off_entry_flux" #"_shut_off_exit_flux"

# Save plot
if demand_function == "periodic_logistic"
    savefig(plt, "traffic_pop_curve_diurnal_test_" * demand_function * "_Nc_$Nc" * "_kc_jam_base_$kc_jam_base" * "_t_end_$my_t_end" * "_L_$my_L" * "_r_$my_r" * "_x0_$my_x0" * "_shift_$my_shift" * "_timescale_$scale" * suffix * ".png")
else
    savefig(plt, "traffic_pop_curve_diurnal_test_" * demand_function * "_Nc_$Nc" * "_kc_jam_base_$kc_jam_base" * "_t_end_$my_t_end" * "_timescale_$scale" * suffix * ".png")
end

display(plt)