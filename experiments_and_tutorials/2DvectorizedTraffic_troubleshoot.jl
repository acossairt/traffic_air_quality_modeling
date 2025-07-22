# Fast call for helper functions to run model
aurora_path_to_help_funcs = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/help_funcs.jl"
marty_path_to_help_funcs = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/help_funcs.jl"
include(aurora_path_to_help_funcs) # CHANGE FOR YOURSELF
using .HelpFuncs
using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra

#=
Set initial conditions, parameters, and other arguments to pass to 
model solver
=#

Np = 2 # Number of patches
Nc = 1 # Number of corridors

# Set initial conditions
kp_init = [1.0; 0.0]
kc_init = zeros(Np, Np, Nc)
u_init = 1
v_init = 0

# Set timescale
my_period = 1440    # 24 hours, in minutes
my_t_end = 1440     # Change if you only want to examine part of a day (< my_period), or multiple days (> my_period)
scale = my_period / 1440 # if using a period other than 1440, this will rescale the other parameters

# Set parameters
# Tolerance for congestion (for now, assume same for patch 1 and patch 2)
my_α = [1]
# Inverse road capacity (for now, assume same for all corridors)
my_kc_jam = 10 * ones(Np, Np, Nc) * scale
# set diagonal values to almost Inf
my_kc_jam[[CartesianIndex(i, i, k) for i in 1:Np, k in 1:Nc]] .= 1e9
display(my_kc_jam)

# Arguments for demand function (which is f() = L / (1 + exp(-r * (x-x0))) )
demand_function = "periodic_logistic" # Other option is "static"
my_L = 0.6
my_r = 100
my_x0 = 0.99
my_shift = 0 * 60   # in units of minutes, set to 0 for no shift

#=
Build model, solve problem, plot results, and save
=#

model, prob, sol, plt = HelpFuncs.build_symbolic_model_diurnal(Np=Np, Nc=Nc,
    my_kp=kp_init, my_kc=kc_init, my_u=u_init, my_v=v_init, my_α=my_α, my_kc_jam=my_kc_jam,
    γ_version=demand_function, my_period=my_period, t_end=my_t_end, my_L=my_L,
    my_r=my_r, my_x0=my_x0, my_shift=my_shift)

# Save plot
if demand_function == "periodic_logistic"
    savefig(plt, "traffic_pop_curve_diurnal_test_" * demand_function * "_Nc_$Nc" * "_t_end_$my_t_end" * "_L_$my_L" * "_r_$my_r" * "_x0_$my_x0" * "_shift_$my_shift" * ".png")
else
    savefig(plt, "traffic_pop_curve_diurnal_test_" * demand_function * "_Nc_$Nc" * "_t_end_$my_t_end.png")
end

display(plt)