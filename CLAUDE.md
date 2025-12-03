# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains dynamical systems models for simulating transportation systems and their impacts on urban air pollution, with New Delhi, India as the primary case study. The models simulate vehicle flows between patches (locations) through corridors (roads) and compute associated emissions.

## Model Architecture

### Core Concept
- **Patches**: Discrete locations (e.g., residential areas, business districts) where vehicles accumulate
- **Corridors**: Connections between patches (e.g., highways, roads) where vehicles travel
- **Conservation Law**: Total vehicle population is conserved across all patches and corridors

### Key State Variables
- `np` / `P_p*`: Population (vehicle count) in each patch
- `nc` / `P_c*`: Population in each corridor, indexed by direction (e.g., `P_c1__p1p2` = vehicles in corridor 1 heading from patch 1 to patch 2)
- `γ` / `d_p*`: Travel demand (proportion wanting to leave a patch)

### Key Parameters
- `alpha*`: Congestion aversion (higher = fewer travelers leave when congested)
- `beta*` / `kc_jam`: Inverse road capacity (higher = lower capacity, more congestion sensitivity)
- `vff`: Free-flow velocity
- `period`: Diurnal cycle length (typically 24 hours)

### Flow Equations
- **Entry flux** (patch → corridor): Depends on demand, congestion level, and congestion aversion
- **Exit flux** (corridor → patch): Depends on corridor population and average velocity (which decreases with congestion)

### Emissions Calculation
Uses a U-shaped curve relating average vehicle speeds to emission rates (g/km), based on the principle that both very low speeds (stop-and-go) and very high speeds produce higher emissions.

## Implementation Variants

### Julia (Primary)
Location: `julia_model/`

Key files:
- `traffic_v2_compact_emissions.jl`: Latest compact model with emissions calculations
- `traffic_model_v2.jl`: Full model with detailed documentation

Dependencies (add to Julia environment):
```julia
using Pkg
Pkg.add(["ModelingToolkit", "DifferentialEquations", "Plots", "LinearAlgebra", "Interpolations"])
```

Run a model:
```bash
julia julia_model/traffic_v2_compact_emissions.jl
```

### XPPAUT
Location: `xppaut_model/`

Key files:
- `ode_files/travel_4_10212024.ode`: Latest ODE file for XPPAUT simulations

Run with XPPAUT:
```bash
./xppaut_model/xppaut xppaut_model/ode_files/travel_4_10212024.ode
```

## Naming Conventions

### XPPAUT Variable Naming
Due to character limits, XPPAUT uses a specific notation:
- `__` represents "from" (origin)
- `_` represents "to" (destination) or "in" (location)
- Example: `F__p1_c1` = "flux from patch 1 to corridor 1"
- Example: `P_c1__p1p2` = "population in corridor 1 traveling from p1 to p2"

### Julia Variable Naming
Uses more explicit array indexing:
- `np[i]`: Population in patch i
- `nc[i,j,k]`: Population in corridor k heading from patch i to patch j
- `γ[i,j]`: Demand for travel from patch i to patch j
