## sample_methods.jl : collection of probability density functions
#
# Representation of a sampling method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## AbstractSamplemethod ##
"""
```julia
AbstractSampleMethod
```

Represents an abstract sample method.
"""
abstract type AbstractSampleMethod end

## MC ##
"""
```julia
MC()
```

Returns a Monte Carlo sample method.

See also: [`QMC`](@ref)
"""
struct MC <: AbstractSampleMethod end

## QMC ##
"""
```julia
QMC()
```

Returns a Quasi-Monte Carlo sample method.

See also: [`MC`](@ref)
"""
struct QMC <: AbstractSampleMethod end
