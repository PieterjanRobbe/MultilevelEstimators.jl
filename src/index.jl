## index.jl : representation of a Level and an Index

## Index ##
"""
    Index(i...)

Represents a multi-index.

```jldoctest
julia> index_0_1 = Index(0,1)
(0,1)

```

See also: [`Level`](@ref)
"""
const Index{d} = NTuple{d,N} where {N<:Integer}

Index(i::N...) where {N<:Integer} = all(i .>= 0) ? ntuple(idx -> i[idx], length(i)) : throw(ArgumentError("in Index(i...), arguments i must be larger then or equal to 0"))

## Level ##
"""
    Level(l)

Represents a level.

# Examples
```jldoctest
julia> level_0 = Level(0)
(0,)

```

See also: [`Index`](@ref)
"""
const Level = Index{1}

Level(i) = Index(i)
