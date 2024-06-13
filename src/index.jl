## index.jl : representation of a Level and an Index
#
# Representation of levels and multi-dimensional indices. 
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

## Index ##
"""
    Index(i::Integer...)

Return a multi-index.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> index = Index(2, 1)
CartesianIndex(2, 1)
```
See also: [`Level`](@ref)
"""
Index = CartesianIndex

# function __init__()
#     show(io::IO, index::Index) = print(io, string("(", join(index.I, ", "), ")"))
#     show(io::IO, ::Type{Index}) = print(io, "Index")
# end

## Level ##
"""
    Level(l::Integer)

Return a level.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> level = Level(2)
2
```
See also: [`Index`](@ref)
"""
Level = CartesianIndex{1}

# show(io::IO, level::Level) = print(io, level[1])
# show(io::IO, ::Type{Level}) = print(io, "Level")

## utilities ##
function diff(index::Index{d}) where d
    D = Dict{Index{d}, Int}()
    Istart = max(zero(index), index - oneunit(index))
    Iend = index
	R = Istart:Iend
	for I in Base.Iterators.take(R, length(R)-1)
        D[I] = isodd(sum(I - index)) ? -1 : 1  
    end
    return D
end

sum(I::CartesianIndex) = sum(Tuple(I))

# if VERSION < v"1.1.0-"
# 	(:)(I::CartesianIndex{N}, J::CartesianIndex{N}) where N = CartesianIndices(map((i,j) -> i:j, Tuple(I), Tuple(J)))
# end
