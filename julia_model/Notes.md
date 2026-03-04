Marty's version (MV) -- traffic_v2_compact_emissions.jl
Aurora's original version (AV) -- traffic_model_v2.jl

Differences
- MV defines all variables and parameters at the top, AV defines them just before the functions where they are first used
- MV does not export variables and parameters
- MV removes most of the descriptive comments
- MV defines f() and g() near the top and together
- MV demand_eqs() is 2 dimensional, defines a value (0) for (1,1) and (2,2) as well
- Some MV functions define their arguments explicitly, while others do not (just use "=")
    - Explicit arguments: f(), g(), v_shifted
    - No explicit arguments: vel_cor, ramp_flux, CorEntryFlux, CorExitFlux, PatchExitFlux, PatchEntryFlux, CorEmissionFlux, TotEmissionFlux
- AV uses a function for EntryFlux() and ExitFlux(), then aggregates in variables \phi_in and \phi_out (for total counts), and finally aggregates for each patch and corridor within the D.(np) and D.(nc) definitions (part of eqs). Meanwhile, MV first calculates CorEntryFlux and CorExitFlux, then uses this to calculate PatchExitFlux and PatchEntryFlux, and finally takes a simple difference to define D.(np) and D.(nc) within eqs.
- MV definition of parameter list is more concise, uses merge, and ModelingToolkit.get_defaults(ode), and Dicts


Things to review (and write in an explanatory document)
- Where do the clock equations come from again? (it's a relationship between sin and cosine, but remember the details)
- Why does MV define some functions with arguments and others not? (Why wasn't it working in my version? Could it be a julia v1.12 thing?)
- How does that ModelingToolkit.get_defaults(ode) part work?


New errors
- Feb 23 everything seemed to run. Feb 25, was prompted to install several packages again. Then, upon run, got error with `+` sign again.
```
ERROR: MethodError: no method matching +(::SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}, ::SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal})

Closest candidates are:
  +(::T, ::Union{Number, AbstractArray{<:Number}, AbstractArray{T}, T}...) where T<:Union{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SafeReal}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}}
   @ SymbolicUtils ~/.julia/packages/SymbolicUtils/OjG8E/src/symbolic_ops/addsub.jl:30
  +(::Any, ::Any, ::Any, ::Any...)
   @ Base operators.jl:587
  +(::Float64, ::T, Any...) where T<:Union{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SafeReal}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}}
   @ SymbolicUtils ~/.julia/packages/SymbolicUtils/OjG8E/src/symbolic_ops/addsub.jl:141
  ...

Stacktrace:
  [1] (::SymbolicUtils.AddWorkerBuffer{SymReal})(terms::Tuple{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}, SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}})
    @ SymbolicUtils ~/.julia/packages/SymbolicUtils/OjG8E/src/symbolic_ops/addsub.jl:80
  [2] add_worker
    @ ~/.julia/packages/SymbolicUtils/OjG8E/src/symbolic_ops/addsub.jl:60 [inlined]
  [3] +
    @ ~/.julia/packages/SymbolicUtils/OjG8E/src/symbolic_ops/addsub.jl:31 [inlined]
  [4] add_sum
    @ ./reduce.jl:24 [inlined]
  [5] _mapreduce(f::typeof(identity), op::typeof(Base.add_sum), ::IndexLinear, A::Matrix{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}})
    @ Base ./reduce.jl:440
  [6] _mapreduce_dim
    @ ./reducedim.jl:367 [inlined]
  [7] mapreduce
    @ ./reducedim.jl:359 [inlined]
  [8] _sum
    @ ./reducedim.jl:1017 [inlined]
  [9] _sum
    @ ./reducedim.jl:1016 [inlined]
 [10] sum(a::Matrix{SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}})
    @ Base ./reducedim.jl:1012
 [11] (::var"#5#6")(i::Int64)
    @ Main ./none:0
 [12] iterate
    @ ./generator.jl:47 [inlined]
 [13] collect(itr::Base.Generator{UnitRange{Int64}, var"#5#6"})
    @ Base ./array.jl:834
 [14] top-level scope
    @ ~/opt/traffic_air_quality_modeling/julia_model/traffic_v2_compact_emissions.jl:64
```