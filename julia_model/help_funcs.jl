## Functions to set up model ##

# help_funcs.jl
module HelpFuncs

export build_symbolic_model, plot_populations  # list the functions you want to make accessible

# Actually build those functions
using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using Interpolations
using IfElse

function build_symbolic_model_diurnal(; Np=2, Nc=2, make_prob=true, make_plot=true, pm, cm, my_α, my_β, my_u=1, my_v=0, my_period, d_version="default", t_end=120, v_f=90)
    #=
    Arguments
        - Np (int): number of patches
        - Nc (int): maximum number of corridors between any```````` two patches
        - pm (array of floats): initial conditions for variable p
        - cm (array of floats): initial conditions for variable c
        - my_α (array of floats): values for parameter α
        - my_β (array of floats): values for parameter β
        - my_u (int): initial conditions for variable u (the cos var)
        - my_v (int): initial conditions for variable v (the sin var)
        - my_period (int): value for parameter period
        - d_version (string): specifies version of demand function to use
            - "default": d1 = 1 always, d2 = 0 always
            - "periodic": d1 peaks at theta = 3pi/4, d2 peaks at theta = 5pi/4
            - "step": d1 = 1 for pi/2 \leq theta \leq pi, 0 everywhere else
                      d2 = 1 for -pi \leq theta \leq -pi/2, 0 everywhere else
        - t_end (int): number of time steps for dynamical solver to take
        - v_f (float): free-flow-velocity (argument to be passed to calc_space_mean_speed_alternative)
    =#

    # Create variables for population model and internal "clock" model
    @variables p(t)[1:Np]             # patches
    @variables c(t)[1:Np, 1:Np, 1:Nc] # corridors
    @variables d(t)[1:Np]             # demand to leave each patch
    @variables u(t)                   # cosine function (for internal clock)
    @variables v(t)                   # sine function (for internal clock)

    # Create parameters
    @parameters α[1:Np]               # tolerance for congestion
    @parameters β[1:Np, 1:Np, 1:Nc]   # inverse road capacity
    @parameters period                # period for internal clock

    # Define demand function (to create "d" for Flux functions)
    theta = atan(v, u)
    if d_version == "periodic_0"
        d_eqs = [
            d[1] ~ max(0, cos(theta - (3pi / 4))),
            d[2] ~ max(0, cos(theta - (5pi / 4))) # seems basically the same as -v?
        ]
    elseif d_version == "periodic_small"
        d_eqs = [
            d[1] ~ max(0.00001, cos(theta - (3pi / 4))),
            d[2] ~ max(0.00001, cos(theta - (5pi / 4))) # seems basically the same as -v?
        ]
    elseif d_version == "step"
        d_eqs = [
            d[1] ~ IfElse.ifelse(theta ≥ pi / 2, IfElse.ifelse(theta ≤ pi, 1.0, 0.0), 0.0),
            d[2] ~ IfElse.ifelse(theta ≥ -pi, IfElse.ifelse(theta ≤ -pi / 2, 1.0, 0.0), 0.0)
        ]
    elseif d_version == "static"
        d_eqs = [
            d[1] ~ 1,
            d[2] ~ 0,
        ]
    else
        println("Please specify a valid version code for the demand function. Defaulting to \"static\"")
        d_eqs = [
            d[1] ~ 1,
            d[2] ~ 0,
        ]
    end

    # Define the corridor flux matrices
    #notice the use of `.` notation before each arithmetic operation, but NOT before `=`
    EnFlx(p, c, d, α, β) = exp.(-α .* β .* c) .* p .* d
    ExFlx(c, β) = exp.(-β .* c) .* c

    # Define equations (for internal clock and population model)
    # notice there are only two equations in the population model, regardless of the size of Np or Nc
    eqs = [
        D(u) ~ u * (1 - u^2 - v^2) - (2 * pi / period) * (v),
        D(v) ~ v * (1 - u^2 - v^2) + (2 * pi / period) * (u),
        D.(p) ~ [sum(ExFlx(c, β)[:, i, :]) for i in 1:Np] .- [sum(EnFlx(p, c, d, α[1], β)[i, :, :]) for i in 1:Np],
        D.(c) ~ collect(-ExFlx(c, β) + EnFlx(p, c, d, α[1], β)),
        d_eqs...   # <-- include the triple dots to splice the equations into the list
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
        prob = ODEProblem(model, [p => pm, c => cm, u => my_u, v => my_v], (0.0, t_end), [α => my_α, β => my_β, period => my_period])
        # Notice we don't pass `d` as a parameter or variable, because it is defined within our equations as an `observed` variable
        # Also, d does not require initial conditions because it is a function of u and v, which do have initial conditions assigned

        if make_plot
            sol = solve(prob, Tsit5())

            # Define accessories
            my_colors = [:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2]
            my_linestyles = [:dash, :dashdotdot]
            legend_loc = Nc > 1 ? :right : :topleft

            # Set up subplots
            plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=18, xguidefontsize=14, yguidefontsize=14, ylabel="population", legend=legend_loc) # palette=palette(my_colors)
            plt2 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="v (km/h)", legend=:right) #, ylims=(35, 50)) # palette=palette(my_colors)
            plt3 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", legend=:right)

            # Plot patch populations (first subplot)
            for i in 1:Np
                plt1 = plot!(plt1, sol, idxs=(p[i]), label="p$i", linewidth=3, color=my_colors[i])
                plt3 = plot!(plt3, sol, idxs=(d[i]), label="d$i", linewidth=3, linestyle=:dash, color=my_colors[i])
            end

            # Calculate average speeds (for second subplot)
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)

            # Track total populations to ensure population is conserved
            total_pop = zeros(length(sol.t))
            total_pop_patches = zeros(length(sol.t))

            # Plot corridor populations (first subplot) and average speeds (second subplot)
            for j in 1:Np # to patch j
                my_colors1 = palette([my_colors[j], :snow2], Np * Nc + 2) # padding because I won't use the first or last colors in this palette
                display(my_colors1)

                total_pop_patches .+= sol[p[j]]
                total_pop .+= sol[p[j]]

                for i in 1:Np # from patch i
                    if j != i
                        for k in 1:Nc # via corridor k
                            this_β = my_β[i, j, k]
                            this_color = my_colors1[i+k] # always skips the first (brightest) shade
                            plt1 = plot!(plt1, sol, idxs=(c[i, j, k]), label="c$k, p$i→p$j (C_jam=1/$this_β)", linewidth=3, linestyle=my_linestyles[k], color=this_color)
                            plt2 = plot!(plt2, sol.t, C_speeds[:, i, j, k], label="u: c$k, p$i → p$j", title="average speeds (assume v_f=90 kmh)", linewidth=3, linestyle=my_linestyles[k], color=this_color)
                            total_pop .+= sol[c[i, j, k]]
                        end
                        #else
                        #    for k in 1:Nc
                        #        total_pop .+= sol[c[i, j, k]]
                        #    end
                        #end
                    end
                end
            end

            # Plot total populations (should always be 1)
            plt1 = plot!(plt1, sol.t, total_pop, label="total population", color="black")
            #plt1 = plot!(plt1, sol.t, total_pop_patches, label="total population (patches only)", color="red")

            if Nc > 1
                title!(plt1, "Daily commute: $Np patches, $Nc corridors")
            else
                title!(plt1, "Daily commute: $Np patches, $Nc corridor")
            end

            plt1 = title!(plt1, "Daily commute: $Np patches, $Nc $(Nc > 1 ? "corridors" : "corridor")")

            # Last subplot, diurnal pattern
            plt3 = plot!(plt3, sol, idxs=(u), label="u", linewidth=3)
            plt3 = plot!(plt3, sol, idxs=(v), label="v", linewidth=3)
            plt3 = plot!(plt3, sol.t, atan.(sol[v], sol[u]), label="theta")
            plt3 = title!(plt3, "Demand function: " * d_version)

            #=
            # Third subplot, emission rates
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)
            C1_emissions = calc_emissions_from_speed(sol[c], C_speeds, interp_fn)
            plt3 = plot(sol.t, C1_emissions, xlabel="t", ylabel="CO2 (G / hour)", label="CO2", title="average emissions rate")
            =#

            # --- Combine into final layout ---
            plt = plot(plt1, plt2, plt3, layout=@layout([a; b; c]), link=:x, size=(800, 900))

            # I want to also plot the sum of all the populations, but not sure how to do that with these idx's?

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

##################################
# Extract variable vals from sol #
##################################

function extract_pvals_from_sol(p_vals, i)
    my_vals = [p_vals[t][i] for t in 1:length(p_vals)]
    return my_vals
end

function extract_cvals_from_sol(c_vals, i, j, k)
    my_vals = [c_vals[t][i, j, k] for t in 1:length(c_vals)]
    return my_vals
end

#####################################
## Calculate speeds from densities ##
#####################################

function calc_space_mean_speed_alternative(v_f, C, my_β; a=1, Np=2, Nc=2)
    #=
        Returns:
            - u_s: float
                average speed for vehicles in a given traffic flow. If negative, return 0.
        Arguments:
            - v_f: float
                free-flow velocity
            - C: float
                vehicle density (in corridor C)
            - C_half:
                threshold value of vehicle density, where u_s = 1/2 * v_f
    =#
    C_jam = 1 ./ my_β
    C_half = C_jam ./ 2
    t_steps = length(C)
    u_s = Array{Float64}(undef, t_steps, Np, Np, Nc) # is this the right order? Or should Nc go first?
    for i in 1:Np
        for j in 1:Np
            for k in 1:Nc
                u_s[:, i, j, k] .= [-(v_f / pi) * atan.(a * (C[t][i, j, k] - C_half[i, j, k])) + (v_f / 2) for t in 1:length(C)]
            end
        end
    end
    #u_s .= [-(v_f ./ pi) .* atan.(a .* (C[t] .- C_half)) .+ (v_f ./ 2) for t in 1:length(C)]
    return max.(u_s, 0.0)
end

######################
# Get emission rates #
######################

function get_emission_rates(x_1) # how to specify that this function should work with vector arguments?
    x_2 = ifelse.(i * f * x_1 .> 60, x_1 - 60 * i * f * x_1, 0.0) # what are i and f?
    x_3 = ifelse.(i * f * x_1 .> 80, x_1 - 80 * i * f * x_1, 0.0)
    Y = 0.027 + 6.8e-05 * x_1 + 0.00032 * x_2 + 0.00050 * x_3
    return Y
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

