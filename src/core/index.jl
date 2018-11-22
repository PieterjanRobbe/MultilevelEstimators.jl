## index.jl : representation of a Level and an Index
#
# Representation of levels and multi-dimensional indices. 
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Index ##
"""
```julia
Index(i...)
```

Returns a multi-index.

# Examples
```jldoctest
julia> index = Index(2,1)
(2,1)

```

See also: [`Level`](@ref)
"""
const Index{d} = NTuple{d,N} where {N<:Integer}

Index(i::Integer...)= all(i .>= 0) ? ntuple(idx -> i[idx], length(i)) : throw(ArgumentError("in Index(i...), arguments i must be larger than or equal to 0"))
Index(i::T...) where {T} = throw(MethodError(Index, i))

## Level ##
"""
```julia
Level(l)
```

Returns a level.

# Examples
```jldoctest
julia> level = Level(2)
(2,)

```

See also: [`Index`](@ref)
"""
const Level = Index{1}

Level(i) = Index(i)

show(io::IO, level::Level) = print(io, level[1])

for f in [:< :> :≤ :≥]
    eval(
         quote
             $f(i::Integer, level::Level) = $f(i, level[1])
             $f(level::Level, i::Integer) = $f(level[1], i)
         end)
end
+(i::Integer, level::Level) = level .+ i
+(level::Level, i::Integer) = level .+ i

## difference ##
function diff(index::Index{d}) where d
    D = Dict{Index{d}, Int}()
    Istart = max.(zero(index), index.-one(index))
    Iend = index
    for I in CartesianIndices(UnitRange.(Istart, Iend))
        if I.I != index
            D[I.I] = isodd(sum(I.I.-index)) ? -1 : 1  
        end 
    end
    return D
end

zero(::Index{d}) where d = ntuple(i->0,d)
one(::Index{d}) where d = ntuple(i->1,d)

## unit ##
unit(i, d) = I[1:d, i]
getindex(u::UniformScaling{Bool}, v::AbstractVector, j::Int) = [u[i,j] for i in v]
