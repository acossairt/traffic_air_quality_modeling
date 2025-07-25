using ModelingToolkit, DifferentialEquations, Plots, JLD2, DataFrames, Arrow, Random

"""Plots mutiple trajectory files with same values for sigma and re for same regions."""

include("../../WESR/WESR_std_model.jl")
include("../../WESR/WESR_analyses.jl")
include("../../WESR/WESR_plots.jl")

using .WESR_std_model_mod
using .WESR_analyses_mod
using .WESR_plots_mod

println("Setting up Model...")
ode = WESR_std_model_mod.construct_ode_system()

# lets define the parameters that we want in this study script

tspan = (0.0, 1000.0); #time span in years.

u0 = [P₁ => 0.24, P₂ => 0.24, K₁ => 0, K₂ => 0,
        G => 2.8, e₁ => 0.0004, e₂ => 0.0004, z => 0]; #initial conditions

p = [r₁ =>
                0.038, Kp₁ => 1.5, #population growth, carrying capacity
        r₂ => 0.042, Kp₂ => 9.7, #population growth, carrying capacity
        ml₁₂ => 10, ml₂₁ => 10, mp₁₂ => 10, mp₂₁ => 10, #migration preference parameters
        r₁₂ => 0, r₂₁ => 0, #migration rates
        α₁ => 0.5, α₂ => 0.5, #capital factor productivity
        a₁ => 2.7, a₂ => 1.7, #total factor productivity
        δ₁ => 0.05, δ₂ => 0.05, #capital entropic decay rates  
        s₁ => 0.25, s₂ => 0.21, #savings rates (of disposable income)
        Cₘ₁ => 0.7, Cₘ₂ => 0.7, #minimum subsistence consumption
        σ₁ => 0.03, σ₂ => 0.03, dmₓ => 100, #climate damages on infrastructure
        re₁ => 0.1, re₂ => 0.1, eb => 0.00004, Tg => 4.2,  #decarb, #carbon intensities
        η => 1, #rate at which decarbonization is initiated
        u => 0.0025, α => 0.1, #earth system carbon uptake and release
        G₁ => 5, G₀ => 4.6, Gₘ => 20 #G1=damage threshold, G₀ = climate threshold, Gₘ = max atm carbon 
];

prob = WESR_std_model_mod.construct_ode_problem(ode, tspan, u0, p)

# construct arrays for changing params
Tg_v = [3.0, 3.4, 3.8, 4.2, 4.6, 5, 5.4, 5.8, 6.2, 6.6, 7.0]
G0_v = [7.0]
sigma1_v = [0.03]
re1_v = [0.1]
combs1 = unique(collect(Iterators.product(Tg_v, G0_v, sigma1_v, re1_v)))

println("Model setup done!")


println("Start runs...")
function prob_func1(prob, i, repeat)
        #rand_Tg = rand(pa)
        #rand_G_0 = rand(Gms)
        #rand_Tg = pa[i]
        #rand_G_0 = Gms[i]
        pnew = [r₁ => 0.038, Kp₁ => 1.5, r₂ => 0.042, Kp₂ => 9.7, ml₁₂ => 10,
                ml₂₁ => 10, mp₁₂ => 10, mp₂₁ => 10, r₁₂ => 0, r₂₁ => 0, α₁ => 0.5, α₂ => 0.5, a₁ => 2.7,
                a₂ => 1.7, δ₁ => 0.05, δ₂ => 0.05, s₁ => 0.25, s₂ => 0.21, Cₘ₁ => 0.7, Cₘ₂ => 0.7,
                σ₁ => combs1[i][3], σ₂ => combs1[i][3], dmₓ => 100, re₁ => combs1[i][4],
                re₂ => combs1[i][4], eb => 0.00004, Tg => combs1[i][1],
                η => 1, u => 0.0025, α => 0.1, G₁ => 5, G₀ => combs1[i][2], Gₘ => 20]
        # push!(rand_values, [rand_Tg, rand_G_0])
        remake(prob, p=pnew) #change Tg and G_0
end

N_samp = length(combs1)


# put it into a function, since according to documentation thats faster
function construct_solve(prob, prob_func, N_samp)
        ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
        return solve(ensemble_prob, trajectories=N_samp, saveat=1)
end

sol1 = construct_solve(prob, prob_func1, N_samp)
println("Runs done!")

# set the backend to plotly to investigate well
plotly()

# create a randomly mixed colormap but with batlow to be consistent
cmap = cgrad(:roma, N_samp, categorical=true)
array = shuffle(1:N_samp)
pp1 = plot()
for i in 1:N_samp
        plot!(pp1, sol1[i][G*100], sol1[i][Y₁+Y₂], xlabel="Climate Threshold G₀ [ppm]",
                ylabel="Total Gross National Income [Trillion \$]", alpha=1.0, legend=false, xlim=[280, 700], ylim=[0, 130], c=cmap[rand(array)])
end
display(pp1)

cmap = cgrad(:roma, N_samp, categorical=true)
array = shuffle(1:N_samp)
pp1 = plot()
for i in 1:N_samp
        plot!(pp1, sol1[i][G*100], sol1[i][Y₁+Y₂], xlabel="Climate Threshold G₀ [ppm]",
                ylabel="Total Gross National Income [Trillion \$]", alpha=1.0, legend=false, xlim=[200, 700], ylim=[0, 200], c=cmap[rand(array)])
end
display(pp1)

# lets investigate the whole solution
plot(sol1[1], legend=true)

# National Incomes
plot(sol1[1][Y₁], label="Y₁")
plot!(sol1[1][Y₂], label="Y₂")
plot!(sol1[1][G*10], label="G*10")

# plot phase space
plot(sol1[1][G*100], sol1[1][Y₁], label="region 1", ylabel="Gross National Income [Trillion \$]")
plot!(sol1[1][G*100], sol1[1][Y₂], label="region 2", xlabel="G₀")
vline!(G0_v * 100, label="Threshold")