# polynomial_fit_emissions.jl
using Pkg
# uncomment if you need to install packages
 Pkg.add("Polynomials")
# Pkg.add("Plots")
 Pkg.add("Statistics")

using Polynomials, Plots, Statistics

# Data you provided
mph_to_kmh = 1.60934
speed_range = range(5*mph_to_kmh, step=5*mph_to_kmh, stop=100*mph_to_kmh)
x = collect(speed_range)    # numeric x vector (length 20)

ygpm = [1200, 950, 700, 500, 425, 350, 325, 310, 309, 308,
     308, 308, 309, 320, 330, 350, 375, 400, 450, 550]

# convert from g/mi to kg/km
ykgpkm = ygpm ./ (1000.0 * mph_to_kmh)

#convert from kg/km to kg/hour
y = collect(ykgpkm) .* x

# LOOCV score for a given polynomial degree (mean squared error)
function loocv_mse(x, y, deg)
    n = length(x)
    tot = 0.0
    for i in 1:n
        xi = copy(x); yi = copy(y)
        deleteat!(xi, i); deleteat!(yi, i)   # leave one out
        p = fit(xi, yi, deg)                 # Polynomials.fit (least squares)
        pred = p(x[i])
        tot += (y[i] - pred)^2
    end
    return tot / n
end

# Try a range of degrees (0..8 by default). Don't exceed n-1.
maxdeg = min(8, length(x)-1)
degrees = 0:maxdeg

scores = [loocv_mse(x, y, d) for d in degrees]

# choose best degree (lowest LOOCV MSE)
bestdeg = degrees[argmin(scores)]
bestdeg_score = minimum(scores)
@show bestdeg, bestdeg_score

# Fit final polynomial on all data
p_best = fit(x, y, bestdeg)

# Diagnostics: training RMS error and LOOCV RMSE
train_rmse = sqrt(mean((y .- p_best.(x)).^2))
loocv_rmse = sqrt(bestdeg_score)
println("Selected degree = $bestdeg")
println("Training RMSE = $(round(train_rmse, sigdigits=6))")
println("LOOCV RMSE     = $(round(loocv_rmse, sigdigits=6))")

# Print polynomial coefficients (lowest → highest)
println("\nPolynomial coefficients (constant term first):")
println(coeffs(p_best))

# Plot data and fitted curve
xq = range(first(x), last(x), length=400)
plot1 = plot(x, y, seriestype=:scatter, ylims=(0, maximum(y)*1.1), 
    xlabel="Speed (km/h)", ylabel="Emissions (kg/hr)")
plot1 = plot!(xq, p_best.(xq), label="poly degree $bestdeg", linewidth=2)
plot1 = title!("Emissions vs speed — polynomial fit (deg=$bestdeg)")

display(plot1)

#plot original data in g/mi vs speed in mph
xq_mph = [5:5:100]
plot2 = scatter(xq_mph, ygpm, 
    label="", title="Emissions vs speed — original data",
    xlabel="Speed (mph)", ylabel ="Emissions (g/mi)")

display(plot2)