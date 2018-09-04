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

Index(i::N...) where {N<:Integer} = all(i .>= 0) ? ntuple(idx -> i[idx], length(i)) : throw(ArgumentError("in Index(i...), arguments i must be larger than or equal to 0"))

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

## difference ##
diff(lvl::Level) = lvl > (0,) ? Dict(lvl.-1 => -1) : Dict{Level,Float64}()

diff(idx::I) where {I<:Index} = begin
    D = Dict{I,Float64}()
    d = length(idx)
    signs = kron([[1,-1] for i=1:d]...)
    for i = 1:length(signs)
        new_idx = idx.-ind2sub(tuple([2 for i = 1:d]...),i).+1
        if new_idx != idx && all(new_idx .> -1)
            D[new_idx] = signs[i]
        end
    end
    return D
end

## unit ##
function unit(T,i,d)
    v = zeros(T,d)
    v[i] = one(T)
    return Index(v...)
end

unit(i,d) = unit(Int64,i,d)
