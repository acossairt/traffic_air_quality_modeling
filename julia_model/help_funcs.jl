## Functions to set up model ##

# help_funcs.jl
module HelpFuncs

export build_symbolic_model, plot_populations  # list the functions you want to make accessible

# Actually build those functions
using ModelingToolkit, DifferentialEquations, Plots, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D
using Interpolations
using IfElse

function demand(u, v; u_low=-0.9, u_high=-0.1, v_low=0.1, v_high=0.9)
    if u > u_low && u < u_high
        du = 1
        if v > v_low && v < v_high
            dv = 1
            d1 = 1
        else
            dv = 0.1
            d1 = 0.1
        end
    else
        du = 0.1
        if v > v_low && v < v_high
            dv = 1
            d1 = 0.1
        else
            dv = 0.1
            d1 = 0.1
        end
    end
    return d1, du, dv
end

function build_diurnal_model(; my_u, my_v, my_period, t_end=120)

    # DEFINE VARS AND PARS FOR DIURNAL MODEL
    @variables u(t)
    @variables v(t)
    @variables d(t)
    @parameters period

    # Define equations (for both population and diurnal model)
    # notice there are only two equations, regardless of the size of Np or Nc
    eqs = [
        D(u) ~ u * (1 - u^2 - v^2) - (2 * pi / period) * (v),
        D(v) ~ v * (1 - u^2 - v^2) + (2 * pi / period) * (u),
        #D(d) ~ d * (1 - u^2 - d^2) - (2 * pi / period) * (d), # Results in interesting shapes which I don't understand
    ]

    # Build model symbolically and solve it
    @mtkbuild model = ODESystem(eqs, t)
    prob = ODEProblem(model, [u => my_u, v => my_v, d => my_v], (0.0, t_end), [period => my_period])
    sol = solve(prob, Tsit5())

    # Define demand as a function of time
    # Broadcast the function to get a vector of 3-tuples
    results = demand.(sol[u], sol[v])
    # Unpack into separate vectors
    out_d1 = getindex.(results, 1)
    out_du = getindex.(results, 2)
    out_dv = getindex.(results, 3)
    println(out_d1)

    # Last subplot, diurnal pattern
    plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=18, xguidefontsize=14, yguidefontsize=14, xlabel="t")
    plt1 = plot!(plt1, sol, idxs=(u), label="u", linewidth=3)
    plt1 = plot!(plt1, sol, idxs=(v), label="v", linewidth=3)
    #plt1 = plot!(plt1, sol, idxs=(d), label="d", linewidth=3)
    plt1 = plot!(plt1, sol.t, out_du, label="du", linewidth=3)
    plt1 = plot!(plt1, sol.t, out_dv, label="dv", linewidth=3)
    plt1 = plot!(plt1, sol.t, out_d1, label="d1", linewidth=3, linestyle=:dash)

    # Check sum
    diurnal_sum = sol[u] .^ 2 + sol[v] .^ 2
    #plt1 = plot!(plt1, sol.t, sol[u] .* (1 .- diurnal_sum), label="u ( 1 - u^2 - v^2)", linewidth=3)
    #plt1 = plot!(plt1, sol.t, sol[v] .* (1 .- diurnal_sum), label="v ( 1 - u^2 - v^2)", linewidth=3)
    #println("\ndiurnal sum:\n", diurnal_sum)
    #println("\nu*(1 - u^2 - v^2):\n", sol[u] .* (1 .- diurnal_sum))
    #println("\nv*(1 - u^2 - v^2):\n", sol[v] .* (1 .- diurnal_sum))

    #osc_term_u = -(2 * pi / my_period) .* sol[v]
    #osc_term_v = (2 * pi / my_period) .* sol[u]
    #plt1 = plot!(plt1, sol.t, osc_term_u, label="osc_term_u (-2pi/period * v)", linewidth=3, linestyle=:dash)
    #plt1 = plot!(plt1, sol.t, osc_term_v, label="osc_term_v (2pi/period * u)", linewidth=3, linestyle=:dash)
    return model, prob, plt1

end

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
    if d_version == "periodic"
        d_eqs = [
            d[1] ~ max(0, cos(theta - (3pi / 4))),
            d[2] ~ max(0, cos(theta - (5pi / 4))) # seems basically the same as -v?
        ]
    elseif d_version == "step"
        d_eqs = [
            d[1] ~ IfElse.ifelse(theta ≥ pi / 2, IfElse.ifelse(theta ≤ pi, 1.0, 0.0), 0.0),
            d[2] ~ IfElse.ifelse(theta ≥ -pi, IfElse.ifelse(theta ≤ -pi / 2, 1.0, 0.0), 0.0)
        ]
    else
        println("Please specify a valid version code for the demand function (options: \"periodic\", \"step\")")
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
            display(palette(my_colors))

            # Set up subplots
            plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=18, xguidefontsize=14, yguidefontsize=14, ylabel="population", legend=legend_loc) # palette=palette(my_colors)
            plt2 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="v (km/h)", legend=:right, ylims=(35, 50)) # palette=palette(my_colors)
            plt3 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", legend=:right)

            # Plot patch populations (first subplot)
            for i in 1:Np
                plt1 = plot!(plt1, sol, idxs=(p[i]), label="p$i", linewidth=3, color=my_colors[i])
                plt3 = plot!(plt3, sol, idxs=(d[i]), label="d$i", linewidth=3, linestyle=:dash, color=my_colors[i])
            end

            # Calculate average speeds (for second subplot)
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)

            # Plot corridor populations (first subplot) and average speeds (second subplot)
            for j in 1:Np # to patch j
                my_colors1 = palette([my_colors[j], :snow2], Np * Nc + 2) # padding because I won't use the first or last colors in this palette
                display(my_colors1)

                # Old palette approach
                #my_colors1 = palette([my_colors[j], :snow2], Np * Nc + 1)
                #my_colors2 = palette([my_colors1[2], :snow2], Np * Nc + 1)
                #plt1 = plot!(plt1, palette=my_colors2)
                #plt2 = plot!(plt2, palette=my_colors2)
                # color=my_colors2[i+k-1]

                for i in 1:Np # from patch i
                    if j != i
                        for k in 1:Nc # via corridor k
                            this_β = my_β[i, j, k]
                            this_color = my_colors1[i+k] # always skips the first (brightest) shade
                            plt1 = plot!(plt1, sol, idxs=(c[i, j, k]), label="c$k, p$i→p$j (C_jam=1/$this_β)", linewidth=3, linestyle=my_linestyles[k], color=this_color)
                            plt2 = plot!(plt2, sol.t, C_speeds[:, i, j, k], label="u: c$k, p$i → p$j", title="average speeds (assume v_f=90 kmh)", linewidth=3, linestyle=my_linestyles[k], color=this_color)
                        end
                    end
                end
            end

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

            # Get demand function
            #results = demand.(sol[u], sol[v])
            # Unpack into separate vectors
            #out_d1 = getindex.(results, 1)
            #out_du = getindex.(results, 2)
            #out_dv = getindex.(results, 3)

            # Check sum
            #diurnal_sum = sol[u] .^ 2 + sol[v] .^ 2
            #plt3 = plot!(plt3, sol.t, diurnal_sum, label="u^2 + v^2", linewidth=3)

            #osc_term_u = -(2 * pi / my_period) .* sol[v]
            #osc_term_v = (2 * pi / my_period) .* sol[u]
            #plt3 = plot!(plt3, sol.t, osc_term_u, label="osc_term_u (-2pi/period * v)", linewidth=3, linestyle=:dash)
            #plt3 = plot!(plt3, sol.t, osc_term_v, label="osc_term_v (2pi/period * u)", linewidth=3, linestyle=:dash)

            #=
            # Third subplot, emission rates
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)
            C1_emissions = calc_emissions_from_speed(sol[c], C_speeds, interp_fn)
            plt3 = plot(sol.t, C1_emissions, xlabel="t", ylabel="CO2 (G / hour)", label="CO2", title="average emissions rate")

            # --- Combine into final layout ---
            plt = plot(plt1, plt2, plt3, layout=@layout([a; b; c]), link=:x, size=(800, 900))

            plt = plot(plt1, plt2, layout=@layout([a; b]), link=:x, size=(800, 600))
            =#

            # --- Combine into final layout ---
            plt = plot(plt1, plt2, plt3, layout=@layout([a; b; c]), link=:x, size=(800, 900))

            # Display or save the plot
            # savefig(final_plot, "my_subplots.png")  # if needed

            # I want to also plot the sum of all the populations, but not sure how to do that with these idx's?

            return model, prob, plt


        else
            return model, prob
        end
    else
        return model, p, c
    end
end

function build_symbolic_model(; Np=2, Nc=2, make_prob=true, make_plot=true, pm, cm, my_α, my_β, my_d, t_end=120, v_f=90)
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
    @parameters β[1:Np, 1:Np, 1:Nc]   # inverse road capacity
    @parameters d[1:Np]   # desire to leave patch i

    # Define the corridor flux matrices
    #notice the use of `.` notation before each arithmetic operation, but NOT before `=`
    EnFlx(p, c, d, α, β) = exp.(-α .* β .* c) .* p .* d
    ExFlx(c, β) = exp.(-β .* c) .* c

    # Define equations
    # notice there are only two equations, regardless of the size of Np or Nc
    eqs = [
        D.(p) ~ [sum(ExFlx(c, β)[:, i, :]) for i in 1:Np] .- [sum(EnFlx(p, c, d, α[1], β)[i, :, :]) for i in 1:Np],
        D.(c) ~ collect(-ExFlx(c, β) + EnFlx(p, c, d, α[1], β)),
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
        prob = ODEProblem(model, [p => pm, c => cm], (0.0, t_end), [α => my_α, β => my_β, d => my_d])

        if make_plot
            sol = solve(prob, Tsit5())

            my_colors = [:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2] #palette([:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2])
            my_linestyles = [:dash, :dashdotdot]
            legend_loc = Nc > 1 ? :right : :topleft
            plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=18, xguidefontsize=14, yguidefontsize=14, ylabel="population", palette=palette(my_colors), legend=legend_loc)
            plt2 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="v (km/h)", palette=palette(my_colors), legend=:right, ylims=(35, 50))

            # First subplot, patch populations
            for i in 1:Np
                plt1 = plot!(plt1, sol, idxs=(p[i]), label="p$i", linewidth=3)
            end

            # Calculate for second subplot, average speeds
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)

            # First and second subplots, corridor patches and average speeds
            for j in 1:Np # to patch j
                my_colors1 = palette([my_colors[j], :snow2], Np * Nc + 1)
                display(my_colors1)
                my_colors2 = palette([my_colors1[2], :snow2], Np * Nc + 1)
                display(my_colors2)
                plt1 = plot!(plt1, palette=my_colors2)
                plt2 = plot!(plt2, palette=my_colors2)
                for i in 1:Np # from patch i
                    if j != i
                        for k in 1:Nc # via corridor k
                            this_β = my_β[i, j, k]
                            plt1 = plot!(plt1, sol, idxs=(c[i, j, k]), label="c$k, p$i→p$j (C_jam=1/$this_β)", linewidth=3, linestyle=my_linestyles[k], color=my_colors2[i+k-1]) # label="c[$i,$j,$k], k_jam=1/$this_β"
                            if j == 2 # TEMP: only plotting speeds for corridors to patch 2
                                plt2 = plot!(plt2, sol.t, C_speeds[:, i, j, k], label="u: c$k, p$i → p$j", title="average speeds (assume v_f=90 kmh)", linewidth=3, linestyle=my_linestyles[k], color=my_colors2[i+k-1])
                            end
                        end
                    end
                end
            end

            if Nc > 1
                title!(plt1, "Daily commute: $Np patches, $Nc corridors")
            else
                title!(plt1, "Daily commute: $Np patches, $Nc corridor")
            end

            plt1 = title!(plt1, "Daily commute: $Np patches, $Nc $(Nc > 1 ? "corridors" : "corridor")")

            #=
            # Third subplot, emission rates
            C_speeds = calc_space_mean_speed_alternative(v_f, sol[c], my_β, a=1, Np=Np, Nc=Nc)
            C1_emissions = calc_emissions_from_speed(sol[c], C_speeds, interp_fn)
            plt3 = plot(sol.t, C1_emissions, xlabel="t", ylabel="CO2 (G / hour)", label="CO2", title="average emissions rate")

            # --- Combine into final layout ---
            plt = plot(plt1, plt2, plt3, layout=@layout([a; b; c]), link=:x, size=(800, 900))

            =#
            plt = plot(plt1, plt2, layout=@layout([a; b]), link=:x, size=(800, 600))

            # Display or save the plot
            # savefig(final_plot, "my_subplots.png")  # if needed

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

#v_f = 90             # free-flow velocity, 90 km/hr, same for C1 and C2
#C1_jam = 1 / β₁      # jam density for C1 (causes avg speed = 0)
#C2_jam = 1 / β₂      # jam density for C2 (causes avg speed = 0)
#C1_half = C1_jam / 2 # threshold density for C1 (causes avg speed = 1/2 free-flow speed)
#C2_half = C2_jam / 2 # threshold density for C2 (causes avg speed = 1/2 free-flow speed)

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

###################################
# OLD VERSION, KEEPING FOR RECORD #
###################################

function calc_space_mean_speed_alternative_old(v_f, C, C_half; a=1)
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
    u_s = -(v_f ./ pi) .* atan.(a .* (C .- C_half)) .+ (v_f ./ 2)
    return max.(u_s, 0.0)
    #@. return u_s > 0 ? u_s : 0
    #return ifelse.(u_s .> 0, u_s, 0.0)
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

################
# Old approach #
################

# Make a U-shaped curve using data from the California paper
start = 5
my_step = 5
stop = 100
mph_to_kmh = 1.60934
speed_arr = collect(start:my_step:stop) * mph_to_kmh # convert mph to kmh
emissions_arr = [1200, 950, 700, 500, 425, 350, 325, 310, 309, 308, 308, 308, 309, 320, 330, 350, 375, 400, 450, 550] * mph_to_kmh

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

