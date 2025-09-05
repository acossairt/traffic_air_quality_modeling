module traffic_plots

using Plots

aurora_path_to_model = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_model.jl"
marty_path_to_model = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_model.jl"
include(aurora_path_to_model)
using .traffic_model

export plot_results

function plot_results(sol)
    #=
    Plot results:
        Subplot 1: population (densities) of patches and corridors vs. time (in hours)
        Subplot 2: average vehicle speeds vs. time (in hours)
        Subplot 3: diurnal clock (`u` and `v`) and demand function `γ` vs. time (in hours)
        Subplot 4: speed-density curve
        Subplot 5: zoom in on rush hour 1
    =#

    # Extract some parameter values which will be displayed or plotted explicitly
    #println(sol.prob.ps)
    my_L = sol.prob.ps[:L]                                    # subplot 3
    my_r = sol.prob.ps[:r]                                    # subplot 3
    my_x_0 = sol.prob.ps[:x_0]                                # subplot 3
    my_shift = sol.prob.ps[:shift]                            # subplot 3
    example_kc_half_jam = sol.prob.ps[:kc_half_jam][1, 2, 1]  # subplot 4
    example_v_f = sol.prob.ps[:v_f][1, 2, 1]                  # subplot 4
    my_λ = sol.prob.ps[:λ]                                    # subplot 4
    #my_version = sol.prob.ps[:version]                        # subplot 2
    my_shift = sol.prob.ps[:shift]                            # subplot 5
    #my_kc_crit = sol.prob.ps[:kc_crit]
    #println("my_version", my_version)

    # Define accessories for plots (colors, linestyles, and legend specs)
    my_colors = [:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2]
    my_linestyles = [:dash, :dashdotdot, :dot]
    legend_loc = Nc > 1 || Np > 2 ? :right : :topleft
    rush_hour_1_peak = my_shift + 360 # to define xlims for subplot 5

    # Set up subplots
    plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=14, xguidefontsize=14, yguidefontsize=14, ylabel="population density (k)", legend=legend_loc, legend_background_color=RGBA(1, 1, 1, 0.5))
    plt2 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="v (km/h)", legend=:right, legend_background_color=RGBA(1, 1, 1, 0.5))
    plt3 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", legend=:right, legend_background_color=RGBA(1, 1, 1, 0.5), palette=:lightrainbow)
    plt4 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="population density (k)", ylabel="avg speeds (v)", xlims=(0, 4 * sol.prob.ps[kc_half_jam][1, 2, 1]))
    plt5 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="population density (k)", xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))
    plt6 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="population density (k)", xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))

    #####################################
    # Subplot 3, diurnal demand pattern #
    #####################################
    plt3 = plot!(plt3, sol.t ./ 60, sol[u], label="u", linewidth=2)
    plt3 = plot!(plt3, sol.t ./ 60, sol[v], label="v", linewidth=2)

    # Annotate subplot 3 (specify demand function)
    my_shift_normalized = my_shift / 60
    # Create parameter text
    param_text = """
        Demand function:
        f = L / (1 + exp^(-r*(x-x0)))
        Parameters:
        L: $my_L | r: $my_r | x0: $my_x_0
        x shifted by: $my_shift_normalized hours
        """

    # Use relative positioning (0-1 scale)
    x_pos = 0.05
    y_pos = 0.15
    xlims = Plots.xlims(plt3)
    ylims = Plots.ylims(plt3)
    x_actual = xlims[1] + (xlims[2] - xlims[1]) * x_pos
    y_actual = ylims[1] + (ylims[2] - ylims[1]) * y_pos

    # Attach annotations to subplot 3
    plt3 = annotate!(plt3, (x_actual, y_actual, text(param_text, :left, 8, :black, RGBA(0, 0, 0, 1))))

    # Set subplot 1 title
    plt3 = title!(plt3, "demand function")

    ###########################################
    # Patch-specific plots (subplots 1 and 3) #
    ###########################################
    for i in 1:Np
        plt1 = plot!(plt1, sol.t ./ 60, sol[kp[i]], label="kp$i", linewidth=3, color=my_colors[i])
        plt3 = plot!(plt3, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
        plt5 = plot!(plt5, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
        plt6 = plot!(plt6, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
    end

    ##############################################
    # Corridor specific plots (subplots 1, 2, 5) #
    ##############################################

    # Track total populations to ensure population conserved within REAL corridors
    total_pop = zeros(length(sol.t))
    total_pop_patches = zeros(length(sol.t))

    # Plot corridor-specific values (subplots 1 and 2)
    for j in 1:Np # to patch j
        my_colors1 = palette([my_colors[j], :snow2], Np * Nc + 2) # padding because I won't use the first or last colors in this palette

        # Update total population counts
        total_pop_patches .+= sol[kp[j]]
        total_pop .+= sol[kp[j]]

        for i in 1:Np # from patch i
            if j != i # skip diagonals (non-existent "loop" corridors)
                for k in 1:Nc # via corridor k
                    this_color = my_colors1[i+k] # always skips the first (brightest) shade

                    # Extract parameter values assigned in this problem
                    this_kc_half_jam = sol.prob.ps[:kc_half_jam][i, j, k] #kc_half_jam[i, j, k]
                    this_v_f = sol.prob.ps[:v_f][i, j, k]
                    this_v_f_kmh = this_v_f * 60 # convert to kmh
                    this_λ = sol.prob.ps[:λ]

                    # Plot corridor population densities (subplot 1)
                    plt1 = plot!(plt1, sol.t ./ 60, sol[kc[i, j, k]],
                        label="kc$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color)

                    # Same for subplot 5 (zoomed in)
                    plt5 = plot!(plt5, sol.t ./ 60, sol[kc[i, j, k]],
                        label="kc$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) ./ 60, (rush_hour_1_peak + 180) ./ 60)) # Need to change the ylims later, in case of different kc_half_jams

                    # Same for subplot 6 (zoomed WAY in)
                    plt6 = plot!(plt6, sol.t ./ 60, sol[kc[i, j, k]],
                        label="kc$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) ./ 60, (rush_hour_1_peak + 180) ./ 60),
                        ylims=(0, 5 * this_kc_half_jam))

                    # Plot jam densities for each corridor (subplot 1)
                    plt1 = hline!(plt1, [this_kc_half_jam],
                        label="kc_half_jam[$i,$j,$k]=$(round(this_kc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color)

                    # Same for subplot 5 (zoomed in)
                    plt5 = hline!(plt5, [this_kc_half_jam],
                        label="kc_half_jam[$i,$j,$k]=$(round(this_kc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))

                    # Same for subplot 6 (zoomed WAY in)
                    plt6 = hline!(plt6, [this_kc_half_jam],
                        label="kc_half_jam[$i,$j,$k]=$(round(this_kc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))

                    # Calculate average speeds... would prefer to do this outside the loop...
                    # NOTE: multiplied by 60 for plotting purposes??
                    my_C_speeds = avg_speed(sol[kc[i, j, k]], this_v_f, this_λ, this_kc_half_jam) .* 60

                    # Plot average vehicle speeds in each corridor (subplot 2)
                    plt2 = plot!(plt2, sol.t ./ 60, my_C_speeds,
                        label="u: c$k, p$i → p$j",
                        linewidth=3, linestyle=my_linestyles[k], color=this_color)

                    # Update total population count
                    total_pop .+= sol[kc[i, j, k]] # If this dips below 1, you are probably losing population to the loops
                end
            end
        end
    end

    # Set subplot titles
    plt1 = title!(plt1, "Daily commute: $Np patches, $Nc $(Nc > 1 ? "corridors" : "corridor")")
    plt2 = title!(plt2, "average speeds (assume v_f=$example_v_f kmh)")
    plt5 = title!(plt5, "First rush hour (zoomed in) with demand")
    plt6 = title!(plt6, "First rush hour (zoomed WAY in)")

    # Plot total populations (within REAL corridors) - should always be 1! (subplot 1)
    plt1 = plot!(plt1, sol.t ./ 60, total_pop, label="total population", color="black")

    ###################################################################
    # Create subplot 4, speed-density curve (example for c1: p1 → p2) #
    ###################################################################
    my_C_range = [0.0:0.001:1.0;]
    example_half_v = avg_speed(example_kc_half_jam, example_v_f, my_λ, example_kc_half_jam) .* 60

    # Test on corridor 1, p1 → p2
    # NOTE: multiplied everything by 60 for plotting??
    my_speeds_test = avg_speed(my_C_range, example_v_f, my_λ, example_kc_half_jam) .* 60
    #avg_speed.(example_v_f, my_λ, my_C_range, example_kc_half_jam) .* 60
    my_speeds_direct = avg_speed(sol[kc[1, 2, 1]], example_v_f, my_λ, example_kc_half_jam) .* 60
    #avg_speed.(example_v_f, my_λ, sol[kc[1, 2, 1]], example_kc_half_jam) .* 60
    #my_speeds_old_test = old_avg_speed.(my_C_range, example_kc_half_jam) .* 60
    #my_speeds_old_direct = old_avg_speed.(sol[kc[1, 2, 1]], example_kc_half_jam) .* 60
    plt4 = plot(my_C_range, my_speeds_test, label="speed-density relation", title="speed-density curve: c1, p1 → p2", ylabel="avg speed", xlabel="density") # add back version
    #plt4 = plot!(plt4, my_C_range, my_speeds_old_test, label="old speed-density relation (λ=0)")
    plt4 = vline!(plt4, [example_kc_half_jam], label="my C_half")
    plt4 = scatter!(plt4, sol[kc[1, 2, 1]], my_speeds_direct, label="my speeds (directly calculated)", linewidth=3)
    plt4 = hline!(plt4, [example_half_v], label="v at kc_half_jam: $(round(example_half_v; digits = 3))")
    #plt4 = plot!(plt4, sol[kc[1, 2, 1]], my_speeds_old_direct, label="my speeds OLD (directly calculated)", linewidth=3, xlims=(0, 10 * sol.prob.ps[kc_half_jam][1, 2, 1])) # Why do I have to put this all the way at the end?
    plt4 = xlims!(plt4, (0, 10 * sol.prob.ps[kc_half_jam][1, 2, 1]))

    #println("is speed ever 0?\n", my_speeds_test)

    #################
    # Finalize plot #
    #################
    layout = @layout [
        a d
        b e
        c f
    ]

    plt = plot(plt1, plt4, plt2, plt5, plt3, plt6, layout=layout, size=(1200, 1200))

    return plt

end # function

end # module