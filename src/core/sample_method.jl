## sample_method.jl : types to represent sample methods (MC and QMC)
#
# Representation of a sampling method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

abstract type AbstractSampleMethod end

## MC ##
"""
    MC()

Return a Monte Carlo sample method.

See also: [`QMC`](@ref)
"""
struct MC <: AbstractSampleMethod end

## QMC ##
"""
    QMC()

Return a Quasi-Monte Carlo sample method.

See also: [`MC`](@ref)
"""
struct QMC <: AbstractSampleMethod end
