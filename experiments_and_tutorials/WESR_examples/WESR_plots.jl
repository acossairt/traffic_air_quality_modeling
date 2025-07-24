module WESR_plots_mod


using Plots


include("WESR_std_model.jl")
using .WESR_std_model_mod

function plot_sol(sol)
    plot(sol)
end

function plot_all_trajs(sol, N_samp, zoom::Bool)
    fig = plot()
    for i in 1:N_samp
        plot!(fig, sol[i][G*100], sol[i][Y₁+Y₂], xlabel="Atmoshperic carbon concentration",
        ylabel="Total GNI, World", alpha=0.5, legend=false)
    end
    if zoom
        xlims!(275,700)
        ylims!(50,155)
    end
    return fig
end


function plot_with_gradient(sol, N_samp, combs, param, name, c)
    """
    sol: ode/sde ensemble solution
    N_samp: number of ensemble runs
    combs: combinations of parameters
    param: parameter which should give the gradient
    c: colormap name
    """
    cmap = cgrad(c, length(param), categorical = true)
    fig=plot()
    for i in 1:N_samp
        index = 0
        for (j, v) in enumerate(param)
            if v == combs[i][2]
            index = j
            end
        end
        if i == 1 || i == N_samp
            G0=combs[i][2]
            plot!(fig,sol[i][G*100],sol[i][Y₁+Y₂], xlabel="Atmoshperic carbon concentration",
            ylabel="Total GNI, World",xlims=[275,700],ylims=[50,155], alpha=0.5, c = cmap[index], label="$name $G0")
        end
        plot!(fig,sol[i][G*100],sol[i][Y₁+Y₂], xlabel="Atmoshperic carbon concentration",
        ylabel="Total GNI, World",xlims=[275,700],ylims=[50,155], alpha=0.5, c = cmap[index], label="")
    end
    return fig
end


function plot_resilience_index_matrix(m, p1, p2, name1, name2)
    """
    m = data from resilience_index_matrix function
    p1 = first parameter that is varied
    p2 = second
    """
    fig = heatmap(p1*100, p2*100, m .* 100, c = cgrad(:vik, rev=true), legend=true, colorbar_title="Resilience Index [%]", clim=(0,100))
    xlabel!(name1)
    ylabel!(name2)
    return fig
end

function plot_gridded_resilience_index_matrix(m, p1, p2, name1, name2)
    """
    m = data from resilience_index_matrix function
    p1 = first parameter that is varied
    p2 = second
    """
    fig = heatmap(p1*100, p2*100, m .* 100, c = cgrad(:vik, rev=true), legend=true, colorbar_title="Resilience Index [%]", clim=(0,100), grid=true)
    xlabel!(name1)
    ylabel!(name2)
    return fig
end

function overall_resilience(compressed_data)
    return sum(compressed_data)/length(compressed_data)
end

function resilience_scatter(sol, p1, p2, name1, name2)
    R_i=[]
    for i in eachindex(sol)
        if sol[i](600)[5] < 5
            push!(R_i, 1)
        else
            push!(R_i, 0)
        end
    end
    markercolor = []
    for i in eachindex(sol)
        if R_i[i] == 1
            push!(markercolor, :blue)
        else
            push!(markercolor, :red)
        end
    end
    fig=plot()
    scatter!(fig, p1, p2, markercolor=markercolor)
    xlabel!(name1)
    ylabel!(name2)
    return fig
end # end of function


function plot_3d_scatter(sol, parameters, names)
    """
    function that plots 4 parameters in a 3d scatter plot, where the 4th determines marker size
    parameters = list of lists of parameter ranges
    names = list of names of the parameters
    """
    x_res = []
    y_res = []
    z_res = []
    x_non = []
    y_non = []
    z_non = []
    four_res = []
    four_non = []
    R_i=[]
    for i in eachindex(sol)
        if sol[i](600)[5] < 5
            push!(R_i, 1)
        else
            push!(R_i, 0)
        end
    end
    for i in eachindex(sol)
        if R_i[i] == 1
            push!(x_res, parameters[1][i])
            push!(y_res, parameters[2][i])
            push!(z_res, parameters[3][i])
            push!(four_res, parameters[4][i])
        else
            push!(x_non, parameters[1][i])
            push!(y_non, parameters[2][i])
            push!(z_non, parameters[3][i])
            push!(four_non, parameters[4][i])
        end
    end
    fig=plot()
    scatter!(fig, x_res, y_res, z_res, markercolor=:blue, markersize=four_res*50, label="Resilient Trajectories", legend=:outerbottom)
    scatter!(fig, x_non, y_non, z_non, markercolor=:red, markersize=four_non*50,label="Non-Resilient Trajectories", legend=:outerbottom)
    xlabel!(names[1])
    ylabel!(names[2])
    zlabel!(names[3])
    return fig
end


function plot_ensembles_in_colors_1(sol, combs, N_traj, N_ens, zoom::Bool)
    # plot all trajs
    cmap=cgrad(:rainbow, Integer(N_traj/N_ens), categorical=true)
    fig=plot()
    k = 1
    for i in 1:N_traj
        if i % (N_ens+1) == 0 #the have every N_ens steps (array starts at 1)
            k += 1
        end
        plot!(fig,
        sol[i][G*100],sol[i][Y₁+Y₂], label=false,
        alpha=0.6, c=cmap[k])
        if i % N_ens == 0
            a = combs[k][1]
            b = combs[k][2]
            lbl = (raw"Tᵧ = " * "$a" * ", G₀ = " * "$b")
            plot!(fig, sol[i][G*100],sol[i][Y₁+Y₂], label=lbl,
            legend=true, c=cmap[k], alpha=0.8)
        end
    end
    plot!(fig, legend=:outerbottom, legendcolumns=3)
    xlabel!("Atmoshperic carbon concentration")
    ylabel!("Total GNI, World")
    if zoom
        xlims!(275,800)
        ylims!(0,170)
    end
    return fig
end


end # end of module