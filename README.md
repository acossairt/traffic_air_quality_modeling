# traffic_air_quality_modeling
Dynamical systems modeling of transportation systems and their impacts on urban air pollution. Case of interest: New Delhi, India.

The latest version of the model (as of November 13, 2024) is `travel_4_1021.ode` in XPPAUT and `travel_4.ipynb` in Julia. It involves two patches connected by two (bidirectional) corridors.

Syntax:
- P^i : Patch i
- C^k_{ij} : Corridor k, lane from patch i to patch j
- F_i^{ck} : Flux from patch i into corridor k (implies direction from patch i to patch j)
- F_{ck}^i : Flux from corridor k into patch i (implies direction from patch j to patch i)
