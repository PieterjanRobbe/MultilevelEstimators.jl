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

## difference ##
diff(lvl::Level) = lvl > (0,) ? Dict(lvl.-1 => -1) : Dict{Level,Float64}()

diff(idx::I) where {I<:Index} = begin
    D = Dict{I,Float64}()
    d = length(idx)
    v = [bit.(1:d,xor(n,n>>1)) for n = 0:2^d-1]
    for i = 1:2^(d-1)
        a = idx.-tuple(v[2*i-1]...)
        b = idx.-tuple(v[2*i]...)
        if all(a.>=0); D[a] = 1; end
        if all(b.>=0); D[b] = -1; end
    end
    return D
end

bit(N,m) = m & (1<<(N-1)) >> (N-1)
