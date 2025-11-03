# Fast call for helper functions to run model
aurora_path_to_model = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_model.jl"
aurora_path_to_plotting_funcs = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_plots.jl"
marty_path_to_model = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_model.jl"
marty_path_to_plotting_funcs = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_plots.jl"
include(aurora_path_to_model) # CHANGE FOR YOURSELF
include(aurora_path_to_plotting_funcs) # CHANGE FOR YOURSELF

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
NumPatches = 2
NumCors = 1

@printf "Number of patches: %d\n" NumPatches
@printf "Number of corridors: %d\n" NumCors

# Set timescale
my_period = 1440         # 24 hours, in minutes
t_end = my_period        # Change if you only want to examine part of a day (< my_period), or multiple days (> my_period)
scale = 1440 / my_period # if using a period other than 1440, this will rescale the other parameters
tspan = (0.0, t_end);    #time span in minutes.

# Set initial conditions
np_init = [150000.0; 0.0]
nc_init = zeros(NumPatches, NumPatches, NumCors)
u_init = 1
v_init = 0
γ_init = [1.0; 0.0]

println("Initial patch densities (np_init):")
display(np_init)
println("Initial corridor densities (nc_init):")
display(nc_init)

u0 = [nc => nc_init, np => np_init, γ => γ_init, u => u_init, v => v_init];

println("Setting up parameters...")

# Parameters for calculating average space-mean speed, and therefore fluxes
my_λ = 2
my_ψ = 24.85
my_Le = 123.92
v_f_base = 117.4 # free-flow velocity in kmh
my_v_f = v_f_base / 60 * ones(NumPatches, NumPatches, NumCors)  # Divided by 60 so that this is km per min
my_v_f[[CartesianIndex(i, i, k) for i in 1:NumPatches, k in 1:NumCors]] .= 0 # no movement in loops
println("Free-flow velocities (my_v_f):")
display(my_v_f)

# Parameters for calculating demand (diurnal function)
my_L = 0.6
my_r = 100
my_x0 = 0.99
my_shift = 0 * 60    # in units of minutes, set to 0 for no shift

# Parameters for population model
nc_half_jam_base = 69715
my_nc_half_jam = nc_half_jam_base * ones(NumPatches, NumPatches, NumCors)
my_nc_half_jam[[CartesianIndex(i, i, k) for i in 1:NumPatches, k in 1:NumCors]] .= 1e9 # set diagonal values to almost Inf

# Parameters for experimenting with flux functions
my_entry_on = 1
my_exit_on = 1

# Make parameters list
p = [v_f => my_v_f, λ => my_λ,
    r => my_r, x_0 => my_x0, L => my_L, shift => my_shift,
    ψ => my_ψ, Le => my_Le,
    nc_half_jam => my_nc_half_jam, period => my_period,
    exit_on => my_exit_on, entry_on => my_entry_on,
];

# Check params in original system
println("Check parameters in original system...")
og_params = parameters(ode)
println(og_params)

# Test
#play1 = traffic_model.EntryFlux(0, 150000, 0.1, 500, 90 / 60, 3, 100, 200, 1)
#println("Practice 1: ", play1)
#play2 = traffic_model.EntryFlux(nc_init, np_init, γ_init, my_nc_half_jam, my_v_f, my_λ, my_Le, my_ψ, my_exit_on)
#println("Practice 2: ", play2)
#println(play2[1, 2, 1])

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

#println(sol[ϕ_in[1, 2, 1]])
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

# Test some functions
#println("Testing calc_avg_speed()")
#example_nc_half_jam = sol.prob.ps[:nc_half_jam][1, 2, 1]
#avg = traffic_plots.avg_speed(example_nc_half_jam)
#println("avg for example_nc_half_jam:", avg)


println("Plotting results...")
#plt = traffic_plots.plot_entry_exit(sol)

plt = traffic_plots.plot_results(sol)

# Save plot
println("Saving plot...")
demand_function = "periodic_logistic"
prefix = "study_basic"
suffix = "_new_avg_speed"

display(plt)
savefig(plt, "studies/graphics/$prefix" * "_traffic_pop_curve_diurnal_test_" * demand_function * "_NumCors_$NumCors" * "_nc_half_jam_denom_$nc_half_jam_base" * "_L_$my_L" * "_r_$my_r" * "_x0_$my_x0" * "_shift_$my_shift" * "_sharpness_$my_λ" * suffix * ".png")

# Display plot
println("Displaying plot...")
display(plt)


# TODO be awesome