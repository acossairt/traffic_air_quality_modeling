## Functions to set up model ##

# help_funcs.jl
module HelpFuncs

export build_symbolic_model, plot_populations  # list the functions you want to make accessible

# Actually build those functions
using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

function build_symbolic_model(; Np=2, Nc=2, make_prob=true, make_plot=true, pm, cm, my_α, my_β, my_d)
    #=
    Arguments
        - Np (int): number of patches
        - Nc (int): maximum number of corridors between any two patches
    =#

    # Create variables
    @variables p(t)[1:Np]           # patches
    @variables c(t)[1:Np, 1:Np, 1:Nc] # corridors

    # Create parameters
    @parameters α[1:Np]   # tolerance for congestion
    @parameters β[1:Nc]   # inverse road capacity
    @parameters d[1:Np]   # desire to leave patch i

    # Define the corridor flux matrices
    #notice the use of `.` notation before each arithmetic operation, but NOT before `=`
    EnFlx(p, c, d, α, β) = exp.(-α .* β .* c) .* p .* d
    ExFlx(c, β) = exp.(-β .* c) .* c

    # Define equations
    # notice there are only two equations, regardless of the size of Np or Nc
    eqs = [
        D.(p) ~ [sum(ExFlx(c, β[1])[:, i, :]) for i in 1:Np] .- [sum(EnFlx(p, c, d, α[1], β[1])[i, :, :]) for i in 1:Np],
        D.(c) ~ collect(-ExFlx(c, β[1]) + EnFlx(p, c, d, α[1], β[1])),
    ]
    #=
    Notice we only used "collect()" in the equation for D.(c). 
    For D.(p), we summed over the columns of the EnFlx matrix and over the columns of the ExFlx matrix, 
    then took their difference. If you don't do this (summing), you will get an error because you will be 
    attempting to save an object of shape (2,2,1) into an object of shape (2,).
    =#

    # Build model symbolically
    @mtkbuild model = ODESystem(eqs, t)

    if make_prob
        prob = ODEProblem(model, [p => pm, c => cm], (0.0, 30), [α => my_α, β => my_β, d => my_d])

        if make_plot
            sol = solve(prob, Tsit5())
            plt = plot()
            for i in 1:Np
                plot!(sol, idxs=(p[i]), label="p[$i]")
            end
            for i in 1:Np
                for j in 1:Np
                    if j != i
                        for k in 1:Nc
                            plot!(sol, idxs=(c[i, j, k]), label="c[$i,$j,$k]")
                        end
                    end
                end
            end
            return model, prob, plt

        else
            return model, prob
        end
    else
        return model, p, c
    end
end

function plot_populations(prob, p, c, Np, Nc)
    sol = solve(prob, Tsit5())
    for i in 1:Np
        plot(sol, idxs=(p[i]), label="p[$i]")
    end
    for i in 1:Np
        for j in 1:Np
            if j != i
                for k in 1:Nc
                    plot!(sol, idxs=(c[i, j, k]), label="c[$i,$j,$k]")
                end
            end
        end
    end
end

end  # module

#=
##############################################
## Convert average speeds to emission rates ##
##############################################

# CHECK IF YOU REALLY NEED ALL THESE PACKAGES
using Plots
using LinearAlgebra
using LaTeXStrings
using Roots
using ModelingToolkit

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

# Calculate emissions from a given array of speeds
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


###############################################
## Convert Emission Rates to Total Emissions ##
###############################################

function integrate_emissions(x, y, a, b)
    # Create an interpolation function over data points
    interp_func = LinearInterpolation(x, y, extrapolation_bc=Line())

    # Integrate the interpolated function from a to b
    result, error = quadgk(interp_func, a, b)

    return result, error
end

=#

