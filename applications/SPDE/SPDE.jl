#=
SPDE.jl : Module for simulating a PDE with random coefficents using various
          methods as implemented in MultilevelEstimators.jl.
This file defines a module that contains methods to initialize the estimator.
Currently implemented:
    - MC        : Monte Carlo
    - QMC       : Quasi-Monte Carlo
    - MLMC      : Multilevel Monte Carlo
    - MLQMC     : Multilevel Quasi-Monte Carlo
    - MIMC      : Multi-Index Monte Carlo
    - MIQMC     : Multi-Index Quasi-Monte Carlo
    - AMIMC     : Adaptive Multi-Index Monte Carlo
    - AMIQMC    : Adaptive Multi-Index Quasi-Monte Carlo
    - MGMLMC    : Multigrid Multilevel Monte Carlo
NOTE: Technically, we are simulating a PDE with random coefficients, but this is too
long to write down, so with a slight abuse of notation, this module is called "SPDE".
=#
module SPDE

## dependencies ##
using Interpolations, Reexport, SimpleMultigrid
@reexport using MultilevelEstimators, GaussianRandomFields

## include statements ##
include("SPDE_init.jl")

include("SPDE_data.jl")

include("SPDE_sample.jl")

include("SPDE_multigrid_sample.jl")

end # module SPDE
