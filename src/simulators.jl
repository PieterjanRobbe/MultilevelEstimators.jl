## simulators.jl : denifition of simulators
# Defines a variety of simulators that implement different multilevel
# algorithms.

## Simulator types
"""
`mutable struct Simulator{S}`

Implements a multilevel simulator.

Examples:
```
mlmc_sim = MultilevelMonteCarloSimulator(sampler,tol)
```
"""
mutable struct Simulator{S}
    sampler::Sampler
    tol::Vector
end

struct MonteCarlo end
struct QuasiMonteCarlo end
struct MultiLevelMonteCarlo end
struct MultiLevelQuasiMonteCarlo end
struct MultiIndexMonteCarlo end
struct MultiIndexQuasiMonteCarlo end
struct MultigridMultiLevelQuasiMonteCarlo end

const MonteCarloSimulator = Simulator{MonteCarlo}
const QuasiMonteCarloSimulator = Simulator{QuasiMonteCarlo}
const MultiLevelMonteCarloSimulator = Simulator{MultiLevelMonteCarlo}
const MultiLevelQuasiMonteCarloSimulator = Simulator{MultiLevelQuasiMonteCarlo}
const MultiIndexMonteCarloSimulator = Simulator{MultiIndexMonteCarlo}
const MultiIndexQuasiMonteCarloSimulator = Simulator{MultiIndexQuasiMonteCarlo}
const MultigridMultiLevelQuasiMonteCarloSimulator = Simulator{MultigridMultiLevelQuasiMonteCarlo}

## Useage
# convenience wrapper for simulators
Simulator{S}(sampler::Sampler,tol::T) where {S,T<:AbstractFloat} = Simulator{S}(sampler,tol,".")
Simulator{S}(sampler::Sampler,tol::T,folder::U)  where {S,T<:AbstractFloat,U<:AbstractString} = Simulator{S}(sampler,[tol],folder)

# TODO check arguments of simulator, maybe have something similar to ConvergenceHistory ???
