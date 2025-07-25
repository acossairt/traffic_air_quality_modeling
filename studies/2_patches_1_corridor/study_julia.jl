# Fast call for helper functions to run model
#aurora_path_to_model = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_model.jl"
#aurora_path_to_plotting_funcs = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_plots.jl"
marty_path_to_model = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_model.jl"
marty_path_to_plotting_funcs = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_plots.jl"
include(marty_path_to_model) # CHANGE FOR YOURSELF
include(marty_path_to_plotting_funcs) # CHANGE FOR YOURSELF

# For some reason doesn't work if I include and also call the module? 
# So I just commented out the include

using Printf, ModelingToolkit, DifferentialEquations, LinearAlgebra, Interpolations
using IfElse, Plots
using .traffic_model
using .traffic_plots

#=
Set initial conditions, parameters, and other arguments to pass to model solver
=#

println("Constructing ode...")
ode = traffic_model.construct_ode_system() # Why does this get called before everything else?
#ode = construct_ode_system() # Do I have to start with the module name?

# Why didn't I need this in the previous version?
#@parameters t
#D = Differential(t)

# ? Do I actually need to do this in hard code?
# Set number of patches and corridors
Np = 2
Nc = 1

@printf "Number of patches: %d\n" Np
@printf "Number of corridors: %d\n" Nc

# Set timescale
my_period = 1440         # 24 hours, in minutes
t_end = my_period        # Change if you only want to examine part of a day (< my_period), or multiple days (> my_period)
scale = 1440 / my_period # if using a period other than 1440, this will rescale the other parameters
tspan = (0.0, t_end);    #time span in minutes.

# Set initial conditions
kp_init = [1.0; 0.0]
kc_init = zeros(Np, Np, Nc)
u_init = 1
v_init = 0
γ_init = [1.0; 0.0]

println("Initial patch densities (kp_init):")
display(kp_init)
println("Initial corridor densities (kc_init):")
display(kc_init)

u0 = [kc => kc_init, kp => kp_init, γ => γ_init, u => u_init, v => v_init];

println("Setting up parameters...")
# Set parameters
my_α = [1]   # Tolerance for congestion (for now, assume same for patch 1 and patch 2)
kc_half_jam_base = 150
my_kc_half_jam = (1 / kc_half_jam_base) / 2 * ones(Np, Np, Nc)  # 1/2 jam density for each corridor
my_kc_half_jam[[CartesianIndex(i, i, k) for i in 1:Np, k in 1:Nc]] .= 1e9 # set diagonal values to almost Inf

# Parameters for calculating demand (diurnal function)
my_L = 0.6
my_r = 100
my_x0 = 0.99
my_shift = 0 * 60    # in units of minutes, set to 0 for no shift

# Parameters for calculating average space-mean speed, and therefore fluxes
my_λ = 1
my_a = 50
v_f_base = 30
my_v_f = v_f_base / 60 * ones(Np, Np, Nc)  # Divided by 60 so that this is km per min
my_v_f[[CartesianIndex(i, i, k) for i in 1:Np, k in 1:Nc]] .= 0 # no movement in loops
println("Free-flow velocities (my_v_f):")
display(my_v_f)

p = [v_f => my_v_f, a => my_a, λ => my_λ,
    r => my_r, x_0 => my_x0, L => my_L, shift => my_shift,
    α => my_α, kc_half_jam => my_kc_half_jam, period => my_period
];

# Simplify structure and set up model
println("Simplifying structure...") # because the system has vector unknowns, we need to simplify it
simple_sys = structural_simplify(ode)

println("Creating problem...")
prob = traffic_model.construct_ode_problem(simple_sys, tspan, u0, p)

println("Model setup done!")

# put it into a function, since according to documentation that's faster
function construct_solve(prob)
    return solve(prob, Tsit5())
end

sol = construct_solve(prob)
println("Runs done!")

#=
Calculate average speeds
    Currently have to create an object `long_kc_half_jam` so that the broadcasting works. 
    Also, even after all that, I still am computing the average speeds locally in my
    plotting functions...
=#

#=
println("Calculate average space-mean speeds...")
N = length(sol)
long_kc_half_jam = [my_kc_half_jam for i in 1:N]
long_v_f = [my_v_f for i in 1:N] # should I use sol.prob.ps[:v_f] instead of my_v_f?
avg_v = avg_speed.(my_λ, long_v_f, my_a, sol[kc], long_kc_half_jam)

println("Speeds")
display(avg_v)

=#

println("Plotting results...")
traffic_plots.plot_results(sol)