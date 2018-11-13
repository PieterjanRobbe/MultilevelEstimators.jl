## index_sets.jl : representation of an index set
#
# Representation of Multilevel and Multi-Index sets. All indices in the index set for a
# given size parameter `sz` can be computed using `get_index_set()`. This is implemented
# using the `filter`-function from `Base.Iterators`.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## Indexset ##
abstract type AbstractIndexSet{d} end

## SL ##
"""
```julia
SL()
```

Returns single-level index set.

See also: [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct SL{d} <: AbstractIndexSet{d} end

SL() = SL{1}()

filter(::SL,sz) = itr -> sum(itr) == sz

## ML ##
"""
```julia
ML()
```

Returns a multi-level index set.

See also: [`SL`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct ML{d} <: AbstractIndexSet{d} end

ML() = ML{1}()

filter(::ML,sz) = itr -> sum(itr) <= sz

## FT ##
"""
```julia
FT(d, [δ=fill(1, d)])
```

Returns a full tensor index set in `d` dimenions with optional weights `δ`.

See also: [`SL`](@ref), [`ML`](@ref), [`TD`](@ref), [`HC`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct FT{d, W<:AbstractVector} <: AbstractIndexSet{d}
    δ::W
end

function FT(d::Integer; δ::Vector{T}=fill(1, max(0, d))) where {T<:Real}
    check_args(d, δ, "FT")
    return FT{d,Vector{T}}(δ)
end

filter(idxset::FT,sz) = itr -> all(idxset.δ .* itr .<= sz) 

## TD ##
"""
```julia
TD(d, [δ=fill(1, d)])
```

Returns a total degree index set in `d` dimenions with optional weights `δ`.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`HC`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct TD{d, W<:AbstractVector} <: AbstractIndexSet{d}
    δ::W
end

function TD(d::Integer; δ::Vector{T}=fill(1, max(0, d))) where {T<:Real}
    check_args(d, δ, "TD")
    return TD{d,Vector{T}}(δ)
end

filter(idxset::TD,sz) = itr -> sum(idxset.δ .* itr) <= sz

## HC ##
"""
```julia
HC(d, [δ=fill(1, d)])
```

Returns a hyperbolic cross index set in `d` dimenions with optional weights `δ`.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct HC{d, W<:AbstractVector} <: AbstractIndexSet{d}
    δ::W
end

function HC(d::Integer; δ::Vector{T}=fill(1, max(0, d))) where {T<:Real}
    check_args(d, δ, "HC")
    return HC{d,Vector{T}}(δ)
end

filter(idxset::HC,sz) = itr -> prod(idxset.δ .* itr .+ 1) <= sz

function check_args(d, δ, name)
    d > 1 ||
    throw(ArgumentError(string("to use ", name, ", dimension must be greater than 1, got ", d)))
    length(δ) == d ||
    throw(ArgumentError(string("to use weighted ", name, ", weights δ must be of length ", d, ", got ", length(δ))))
    all(δ .> 0) ||
    throw(ArgumentError(string("to use weighted ", name, ", weights δ must be positive")))
    all(isfinite.(δ)) ||
    throw(ArgumentError(string("to use weighted ", name, ", weight must be finite, got ", δ)))
end

## AD ##
"""
```julia
AD(d)
```

Returns an adaptive index set in `d` dimenions.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`MG`](@ref)
"""
struct AD{d} <: AbstractIndexSet{d} end

AD(d::Integer) = d <= 1 ? throw(BoundsError(string("to use AD, dimension must be greater than 1, got ", d))) : AD{d}()

## MG ##
"""
```julia
MG(I<:AbstractIndexSet)
```

Returns a Multigrid wrapper index set for the index set I.

# Examples
```jldoctest
julia> MG(ML())
MG{ML}

julia> MG(TD(2))
MG{TD{2}}
```

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`AD`](@ref)
"""
struct MG{d, I} <: AbstractIndexSet{d}
    idxset::I
end

MG(idxset::I) where {I<:AbstractIndexSet} = MG{ndims(idxset), I}(idxset)

## Union types ##
AbstractML = Union{ML, MG{ML}}
AbstractMI = Union{TD, FT, HC, AD, MG{TD}, MG{FT}, MG{HC}, MG{ML}}
AbstractAD = Union{AD, MG{AD}}

## ndims ##
ndims(::AbstractIndexSet{d}) where {d} = d

## get_index_set ##
"""
```julia
get_index_set(idxset, sz)
```

Returns all indices in the index set `idxset` for a given size paraneter `sz`.

# Examples
```jldoctest
julia> get_index_set(TD(2), 2)
6-element Array{Tuple{Int64,Int64},1}:
(0, 0)
(1, 0)
(2, 0)
(0, 1)
(1, 1)
(0, 2)

```
"""
function get_index_set(idxset::AbstractIndexSet{d}, sz::Integer) where d
    sz >= 0 || throw(ArgumentError("index set size parameter cannot be negative"))
    tensor_grid = Base.product(UnitRange.(0, ntuple(i -> sz+1, d))...)
    filtered_grid = Base.Iterators.filter(filter(idxset, sz), tensor_grid)
    collect(filtered_grid)
end

get_index_set(idxset::MG, sz::Integer) = get_index_set(idxset.idxset, sz)
get_index_set(idxset::AD, sz::Integer) = throw(ArgumentError("get_index_set not implemented for index sets of type AD"))

## is_admissable ##
function is_admissable(idxset::Set{Index{d}}, index::Index{d})  where d
    if in(index, idxset)
        return false
    else
        for i in 1:d
            bw_neighbour = index .- unit(i, d) # backward neighbour
            if !( in(bw_neighbour, idxset) || bw_neighbour[i] < 0 )
                return false
            end
        end
        return true
    end
end

## output formatting ##
show(io::IO, idxset::AbstractIndexSet) = print(io, shortname(idxset))

for name in ["SL" "ML"]
    @eval shortname(::$(Symbol(name))) = $(name)
end

for name in ["TD" "FT" "HC" "AD"]
    @eval shortname(idxset::$(Symbol(name))) = string($(name), "{", ndims(idxset), "}")
end

shortname(idxset::MG) = string("MG{", shortname(idxset.idxset), "}")
