# Configurable variant of traffic_v2_compact_emissions.jl.
# Usage: julia traffic_configurable.jl <config.json> <output.json>
# Requires: ModelingToolkit, DifferentialEquations, JSON3

using ModelingToolkit, DifferentialEquations, JSON3, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

if length(ARGS) < 2
    error("Usage: julia traffic_configurable.jl <config.json> <output.json>")
end
config_path = ARGS[1]
output_path = ARGS[2]

config = JSON3.read(read(config_path, String))

# patch id (user-facing) -> 1-based array index
patch_ids    = [Int(p.id) for p in config.patches]
patch_index  = Dict(id => i for (i, id) in enumerate(patch_ids))
patch_labels = Dict(Int(p.id) => String(p.label) for p in config.patches)
NumPatches   = length(patch_ids)
NumCors      = 1   # one corridor per directed pair for now

# parameter value matrices, with safe defaults where edges are absent.
# nc_half_ff/vsharp/etc. are nonzero to avoid div-by-zero in g(); γ is gated to 0
# for missing edges so demand-side defaults don't matter.
vff_v        = zeros(NumPatches, NumPatches, NumCors)
nc_half_ff_v = ones(NumPatches, NumPatches, NumCors)
vsharp_v     = ones(NumPatches, NumPatches, NumCors)
onff_v       = zeros(NumPatches, NumPatches, NumCors)
on_half_v    = ones(NumPatches, NumPatches, NumCors)
onsharp_v    = ones(NumPatches, NumPatches, NumCors)
demhf_v      = ones(NumPatches, NumPatches)
shift_v      = zeros(NumPatches, NumPatches)
dur_v        = ones(NumPatches, NumPatches)
dsharp_v     = ones(NumPatches, NumPatches)
bgd_v        = zeros(NumPatches, NumPatches)

edge_pairs = Set{Tuple{Int,Int}}()
for e in config.edges
    i = patch_index[Int(e.from)]
    j = patch_index[Int(e.to)]
    push!(edge_pairs, (i, j))
    vff_v[i,j,1]        = e.vff
    nc_half_ff_v[i,j,1] = e.nc_half_ff
    vsharp_v[i,j,1]     = e.vsharp
    onff_v[i,j,1]       = e.onff
    on_half_v[i,j,1]    = e.on_half
    onsharp_v[i,j,1]    = e.onsharp
    demhf_v[i,j]        = e.demhf
    shift_v[i,j]        = e.shift
    dur_v[i,j]          = e.dur
    dsharp_v[i,j]       = e.dsharp
    bgd_v[i,j]          = e.bgd
end

period_v = Float64(config.sim.period)
psi_v    = Float64(config.sim.psi)
L_v      = Float64(config.sim.L)
tspan    = (Float64(config.sim.tspan[1]), Float64(config.sim.tspan[2]))
saveat_v = Float64(config.sim.saveat)

# state and auxiliary variables
@variables np(t)[1:NumPatches]
@variables nc(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
@variables γ(t)[1:NumPatches, 1:NumPatches]
@variables v̄(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
@variables ec(t)[1:NumPatches, 1:NumPatches, 1:NumCors]
@variables et(t)
@variables u(t) v(t)

# parameters (assigned numeric values further down via p_vals dict)
@parameters period psi L
@parameters vff[1:NumPatches,1:NumPatches,1:NumCors]
@parameters nc_half_ff[1:NumPatches,1:NumPatches,1:NumCors]
@parameters vsharp[1:NumPatches,1:NumPatches,1:NumCors]
@parameters onff[1:NumPatches,1:NumPatches,1:NumCors]
@parameters on_half[1:NumPatches,1:NumPatches,1:NumCors]
@parameters onsharp[1:NumPatches,1:NumPatches,1:NumCors]
@parameters demhf[1:NumPatches,1:NumPatches]
@parameters shift[1:NumPatches,1:NumPatches]
@parameters dur[1:NumPatches,1:NumPatches]
@parameters dsharp[1:NumPatches,1:NumPatches]
@parameters bgd[1:NumPatches,1:NumPatches]

f(x, r, x0, m) = m / (1 + exp(-r * (x - x0)))
# Clamp x at 0 so small negative numerical noise on non-negative quantities
# (populations, densities) doesn't blow up under non-integer exponents `c`.
g(x, a, b, c) = a * exp(-log(2) * (max(x, 0) / b)^c)

v_shifted(s) = v * cos(2π * s / period) - u * sin(2π * s / period)

clock_eqs = [
    D(u) ~ u * (1 - u^2 - v^2) - (2π / period) * v,
    D(v) ~ v * (1 - u^2 - v^2) + (2π / period) * u,
]

demand_eqs = Equation[]
for i in 1:NumPatches, j in 1:NumPatches
    if (i, j) in edge_pairs
        push!(demand_eqs,
            γ[i,j] ~ (1 - g(np[i], 1.0, demhf[i,j], 4.0)) *
                    f(v_shifted(shift[i,j]), dsharp[i,j], dur[i,j], 1.0) +
                    bgd[i,j])
    else
        push!(demand_eqs, γ[i,j] ~ 0)
    end
end

vel_cor   = g.(nc ./ (psi * L), vff, nc_half_ff, vsharp)
ramp_flux = g.(nc ./ (psi * L), onff, on_half, onsharp)

CorEntryFlux = ramp_flux .* γ
CorExitFlux  = (nc ./ L) .* vel_cor

PatchExitFlux  = [sum(collect(CorEntryFlux[i, :, :])) for i in 1:NumPatches]
PatchEntryFlux = [sum(collect(CorExitFlux[:, i, :]))  for i in 1:NumPatches]

ec_coef = [0.9449071207432054, 0.8761499774070579,
           -0.029999531967464826, 0.0004680693625443692,
           -3.1760158923592917e-6, 8.050724090389599e-9]
em_fit(x) = ec_coef[1] + ec_coef[2]*x + ec_coef[3]*x^2 + ec_coef[4]*x^3 +
            ec_coef[5]*x^4 + ec_coef[6]*x^5

CorEmissionFlux = em_fit.(vel_cor) .* nc
TotEmissionFlux = sum(collect(CorEmissionFlux[:,:,:]))

eqs = [
    D.(np) ~ PatchEntryFlux .- PatchExitFlux,
    D.(nc) ~ CorEntryFlux .- CorExitFlux,
    demand_eqs...,
    clock_eqs...,
    v̄ ~ vel_cor,
    ec ~ CorEmissionFlux,
    D(et) ~ TotEmissionFlux
]

@named ode = ODESystem(eqs, t)

# parameter values
p_vals = Dict{Any,Any}(period => period_v, psi => psi_v, L => L_v)
for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors
    p_vals[vff[i,j,k]]        = vff_v[i,j,k]
    p_vals[nc_half_ff[i,j,k]] = nc_half_ff_v[i,j,k]
    p_vals[vsharp[i,j,k]]     = vsharp_v[i,j,k]
    p_vals[onff[i,j,k]]       = onff_v[i,j,k]
    p_vals[on_half[i,j,k]]    = on_half_v[i,j,k]
    p_vals[onsharp[i,j,k]]    = onsharp_v[i,j,k]
end
for i in 1:NumPatches, j in 1:NumPatches
    p_vals[demhf[i,j]]  = demhf_v[i,j]
    p_vals[shift[i,j]]  = shift_v[i,j]
    p_vals[dur[i,j]]    = dur_v[i,j]
    p_vals[dsharp[i,j]] = dsharp_v[i,j]
    p_vals[bgd[i,j]]    = bgd_v[i,j]
end

# initial conditions
u0 = Dict{Any,Any}()
np_init = zeros(NumPatches)
for p in config.patches
    np_init[patch_index[Int(p.id)]] = Float64(p.initial_pop)
end
for i in 1:NumPatches
    u0[np[i]] = np_init[i]
end
for i in 1:NumPatches, j in 1:NumPatches, k in 1:NumCors
    u0[nc[i,j,k]] = 0.0
end
u0[u]  = 1.0
u0[v]  = 0.0
u0[et] = 0.0

prob = ODEProblem(mtkcompile(ode), merge(u0, p_vals), tspan)
sol  = solve(prob, Tsit5(), saveat=saveat_v)

# write results
n = length(sol.t)
patches_out = Dict{String,Any}()
for i in 1:NumPatches
    pid = patch_ids[i]
    patches_out[string(pid)] = Dict(
        "label" => patch_labels[pid],
        "np"    => [sol[np[i]][k] for k in 1:n]
    )
end

edges_out = Dict{String,Any}()
for (i, j) in edge_pairs
    key = "$(patch_ids[i])->$(patch_ids[j])"
    edges_out[key] = Dict(
        "from" => patch_ids[i],
        "to"   => patch_ids[j],
        "nc"   => [sol[nc[i,j,1]][k] for k in 1:n],
        "v"    => [sol[v̄[i,j,1]][k] for k in 1:n],
        "ec"   => [sol[ec[i,j,1]][k] for k in 1:n],
    )
end

result = Dict(
    "time"    => sol.t,
    "patches" => patches_out,
    "edges"   => edges_out,
    "et"      => [sol[et][k] for k in 1:n]
)

open(output_path, "w") do io
    JSON3.write(io, result)
end

println("OK: wrote $output_path ($(n) timepoints)")
