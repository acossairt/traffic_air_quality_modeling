module traffic_plots

using Plots
using Measures

aurora_path_to_model = "/Users/auroracossairt/opt/traffic_air_quality_modeling/julia_model/traffic_model.jl"
marty_path_to_model = "/Users/janderie/work/cmpting/git/traffic_air_quality_modeling/julia_model/traffic_model.jl"
include(aurora_path_to_model)
using .traffic_model

export plot_results

function plot_results(sol)
    #println(sol[ϕ_in])
    #=
    Plot results:
        Subplot 1: population (densities) of patches and corridors vs. time (in hours)
        Subplot 2: average vehicle speeds vs. time (in hours)
        subplot 3: entry and exit fluxes (per corridor)
        subplot 4: diurnal clock (`u` and `v`) and demand function `γ` vs. time (in hours)
        subplot 5: speed-density curve
        subplot 6: zoom in on rush hour 1
    =#

    # Extract some parameter values which will be displayed or plotted explicitly
    #println(sol.prob.ps)
    my_L = sol.prob.ps[:L]                                    # subplot 4
    my_r = sol.prob.ps[:r]                                    # subplot 4
    my_x_0 = sol.prob.ps[:x_0]                                # subplot 4
    my_shift = sol.prob.ps[:shift]                            # subplot 4
    example_nc_half_jam = sol.prob.ps[:nc_half_jam][1, 2, 1]  # subplot 5
    example_v_f = sol.prob.ps[:v_f][1, 2, 1]                  # subplot 5
    my_λ = sol.prob.ps[:λ]                                    # subplot 5
    my_ψ = sol.prob.ps[:ψ]
    my_Le = sol.prob.ps[:Le]
    #my_version = sol.prob.ps[:version]                       # subplot 2
    my_shift = sol.prob.ps[:shift]                            # subplot 6
    #my_nc_crit = sol.prob.ps[:nc_crit]
    #println("my_version", my_version)

    println("my_ψ: ", my_ψ, "\nmy_Le: ", my_Le)

    # Define accessories for plots (colors, linestyles, and legend specs)
    my_colors = [:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2]
    my_linestyles = [:dash, :dashdotdot, :dot]
    legend_loc = NumCors > 1 || NumPatches > 2 ? :outerright : :outertopleft
    rush_hour_1_peak = my_shift + 360 # to define xlims for subplot 6

    # Set up subplots
    plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=14, xguidefontsize=14, yguidefontsize=14, ylabel="population (# cars)", legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5))
    plt2 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="flux (# cars)", legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5)) # right
    plt3 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="v (km/h)", legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5)) # right
    plt4 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5), palette=:lightrainbow) # right
    plt5 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, ylabel="avg speed (v)", xlabel="population density (k)", xlims=(0, 4 * sol.prob.ps[nc_half_jam][1, 2, 1]), legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5))
    plt6 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="population density (k)", xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60), legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5))
    plt7 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, xguidefontsize=14, yguidefontsize=14, xlabel="t", ylabel="population density (k)", xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60), legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5))

    #####################################
    # subplot 4, diurnal demand pattern #
    #####################################
    plt4 = plot!(plt4, sol.t ./ 60, sol[u], label="u", linewidth=2)
    plt4 = plot!(plt4, sol.t ./ 60, sol[v], label="v", linewidth=2)

    # Annotate subplot 4 (specify demand function)
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
    xlims = Plots.xlims(plt4)
    ylims = Plots.ylims(plt4)
    x_actual = xlims[1] + (xlims[2] - xlims[1]) * x_pos
    y_actual = ylims[1] + (ylims[2] - ylims[1]) * y_pos

    # Attach annotations to subplot 4
    plt4 = annotate!(plt4, (x_actual, y_actual, text(param_text, :left, 8, :black, RGBA(0, 0, 0, 1))))

    # Set subplot 1 title
    plt4 = title!(plt4, "demand function")

    ###########################################
    # Patch-specific plots (subplots 1 and 3) #
    ###########################################
    for i in 1:NumPatches
        plt1 = plot!(plt1, sol.t ./ 60, sol[np[i]], label="np$i", linewidth=3, color=my_colors[i])
        plt4 = plot!(plt4, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
        plt6 = plot!(plt6, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
        plt7 = plot!(plt7, sol.t ./ 60, sol[γ[i]], label="γ$i", linewidth=3, linestyle=:dot, color=my_colors[i])
    end

    #################################################
    # Corridor specific plots (subplots 1, 2, 3, 5) #
    #################################################

    # Track total populations to ensure population conserved within REAL corridors
    total_pop = zeros(length(sol.t))
    total_pop_patches = zeros(length(sol.t))

    # Plot corridor-specific values (subplots 1 and 2)
    for j in 1:NumPatches # to patch j
        my_colors1 = palette([my_colors[j], :snow2], NumPatches * NumCors + 2) # padding because I won't use the first or last colors in this palette

        # Update total population counts
        total_pop_patches .+= sol[np[j]]
        total_pop .+= sol[np[j]]

        for i in 1:NumPatches # from patch i
            if j != i # skip diagonals (non-existent "loop" corridors)
                for k in 1:NumCors # via corridor k
                    this_color = my_colors1[i+k] # always skips the first (brightest) shade

                    # Extract parameter values assigned in this problem
                    this_nc_half_jam = sol.prob.ps[:nc_half_jam][i, j, k] #nc_half_jam[i, j, k]
                    this_v_f = sol.prob.ps[:v_f][i, j, k]
                    this_v_f_kmh = this_v_f * 60 # convert to kmh
                    this_λ = sol.prob.ps[:λ]

                    # Plot corridor population densities (subplot 1)
                    plt1 = plot!(plt1, sol.t ./ 60, sol[nc[i, j, k]],
                        label="c$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color)

                    # Same for subplot 6 (zoomed in)
                    plt6 = plot!(plt6, sol.t ./ 60, sol[nc[i, j, k]],
                        label="c$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) ./ 60, (rush_hour_1_peak + 180) ./ 60)) # Need to change the ylims later, in case of different nc_half_jams

                    # Same for subplot 7 (zoomed WAY in)
                    plt7 = plot!(plt7, sol.t ./ 60, sol[nc[i, j, k]],
                        label="c$k, p$i→p$j", linewidth=3,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) ./ 60, (rush_hour_1_peak + 180) ./ 60),
                        ylims=(0, 5 * this_nc_half_jam))

                    # Plot jam densities for each corridor (subplot 1)
                    plt1 = hline!(plt1, [this_nc_half_jam],
                        label="nc_half_jam[$i,$j,$k]=$(round(this_nc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color)

                    # Same for subplot 6 (zoomed in)
                    plt6 = hline!(plt6, [this_nc_half_jam],
                        label="nc_half_jam[$i,$j,$k]=$(round(this_nc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))

                    # Same for subplot 7 (zoomed WAY in)
                    plt7 = hline!(plt7, [this_nc_half_jam],
                        label="nc_half_jam[$i,$j,$k]=$(round(this_nc_half_jam; digits = 3))", linewidth=1,
                        linestyle=my_linestyles[k], color=this_color,
                        xlims=((rush_hour_1_peak - 180) / 60, (rush_hour_1_peak + 180) / 60))

                    # Calculate average speeds... would prefer to do this outside the loop...
                    # NOTE: multiplied by 60 for plotting purposes??
                    #practice_calc_b(nc_half_jam, λ) = (nc_half_jam) .* (1 ./ log(2)) .^ (1 ./ λ)
                    #practice_avg_speed(nc, ψ, Le, v_f, λ, nc_half_jam) = v_f .* exp.(-1 .* (nc ./ practice_calc_b(nc_half_jam, λ)) .^ λ)

                    my_C_speeds = avg_speed(sol[nc[i, j, k]], my_ψ, my_Le, this_v_f, this_λ, this_nc_half_jam) .* 60
                    #my_C_speeds = avg_speed(sol[nc[i, j, k]]) .* 60

                    #println("my_C_speeds:\n", my_C_speeds)

                    # Plot average vehicle speeds in each corridor (subplot 2)
                    plt3 = plot!(plt3, sol.t ./ 60, my_C_speeds,
                        label="u: c$k, p$i → p$j",
                        linewidth=3, linestyle=my_linestyles[k], color=this_color)

                    # Plot entry and exit fluxes for each corridor and direction
                    plt2 = plot!(plt2, sol.t ./ 60, sol[ϕ_in[i, j, k]], linewidth=3,
                        linestyle=:dot, color=this_color, label="Entry Flux: c$k, p$i→p$j")
                    plt2 = plot!(plt2, sol.t ./ 60, sol[ϕ_out[i, j, k]], linewidth=3,
                        linestyle=:dash, color=this_color, label="Exit Flux: c$k, p$i→p$j")

                    # Update total population count
                    total_pop .+= sol[nc[i, j, k]] # If this dips below 1, you are probably losing population to the loops
                end
            end
        end
    end

    # Set subplot titles
    plt1 = title!(plt1, "Daily commute: $NumPatches patches, $NumCors $(NumCors > 1 ? "corridors" : "corridor")")
    plt2 = title!(plt2, "Entry and Exit fluxes")
    plt3 = title!(plt3, "average speeds (assume v_f=$(round(example_v_f; digits = 2)) km/min)")
    plt6 = title!(plt6, "First rush hour (zoomed in) with demand")
    plt7 = title!(plt7, "First rush hour (zoomed WAY in)")

    # Plot total populations (within REAL corridors) - should always be 1! (subplot 1)
    plt1 = plot!(plt1, sol.t ./ 60, total_pop, label="total population", color="black")

    ###################################################################
    # Create subplot 5, speed-density curve (example for c1: p1 → p2) #
    ###################################################################
    my_C_range = [0.0:0.001:1.0;]
    example_half_v = avg_speed(example_nc_half_jam, my_ψ, my_Le, example_v_f, my_λ, example_nc_half_jam) .* 60
    #example_half_v = avg_speed(example_nc_half_jam) .* 60

    # Test on corridor 1, p1 → p2
    # NOTE: multiplied everything by 60 for plotting??
    my_speeds_test = avg_speed(my_C_range, my_ψ, my_Le, example_v_f, my_λ, example_nc_half_jam) .* 60
    my_speeds_direct = avg_speed(sol[nc[1, 2, 1]], my_ψ, my_Le, example_v_f, my_λ, example_nc_half_jam) .* 60
    #my_speeds_test = avg_speed(my_C_range) .* 60
    #my_speeds_direct = avg_speed(sol[nc[1, 2, 1]]) .* 60
    plt5 = plot(my_C_range, my_speeds_test, label="speed-density relation", title="speed-density curve: c1, p1 → p2") # add back version
    plt5 = vline!(plt5, [example_nc_half_jam], label="my C_half")
    plt5 = scatter!(plt5, sol[nc[1, 2, 1]], my_speeds_direct, label="my speeds (directly calculated)", linewidth=3)
    plt5 = hline!(plt5, [example_half_v], label="v at nc_half_jam: $(round(example_half_v; digits = 3))")
    plt5 = xlims!(plt5, (0, 10 * sol.prob.ps[nc_half_jam][1, 2, 1]))

    #################
    # Finalize plot #
    #################

    #layout = @layout [
    #    a d
    #    b e
    #    c f
    #    g h
    #]

    layout = grid(4, 2, heights=[0.25, 0.25, 0.25, 0.25], widths=[0.5, 0.5])

    plt = plot(plt1, plt5, plt2, plt6, plt3, plt7, plt4, layout=layout, size=(1800, 2400))

    return plt

end # function

function plot_entry_exit(sol)
    # Define accessories for plots (colors, linestyles, and legend specs)
    my_colors = [:darkcyan, :chocolate2, :purple, :dodgerblue4, :orangered2]
    my_linestyles = [:dash, :dashdotdot, :dot]
    legend_loc = NumCors > 1 || NumPatches > 2 ? :right : :topleft

    plt1 = plot(legendfontsize=10, xtickfontsize=12, ytickfontsize=12, titlefontsize=14, xguidefontsize=14, yguidefontsize=14, ylabel="flux", legend=:outerbottom, legend_background_color=RGBA(1, 1, 1, 0.5))

    for k in 1:NumCors # via corridor k
        plt1 = plot!(plt1, sol.t ./ 60, sol[ϕ_in[k]], label="Entry Flux")
        plt1 = plot!(plt1, sol.t ./ 60, sol[ϕ_out[k]], label="Exit Flux")
    end # loop

    return plt1
end # function

end # module