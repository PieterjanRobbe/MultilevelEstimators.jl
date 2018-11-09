## index_set.jl : representation of an IndexSet

## Indexset ##
"""
    IndexSet{d}

Supertype for an index set in `d` dimensions.
"""
abstract type IndexSet{d} end

"""
    SL()

A single-level index set.
"""
struct SL{d} <: IndexSet{d} end

SL() = SL{1}()

filter(::SL,sz) = itr -> sum(itr) == sz

"""
    ML()

A multi-level index set.
"""
struct ML{d} <: IndexSet{d} end

ML() = ML{1}()

filter(::ML,sz) = itr -> sum(itr) <= sz

"""
    FT(d, δ=ones(d))

A full tensor index set in `d` dimenions with weights `δ`.
"""
struct FT{d, W<:AbstractVector} <: IndexSet{d}
    δ::W
end

function FT(d::Integer; δ::Vector{T}=ones(d)) where {T<:Real}
    d > 1 || throw(BoundsError("to use FT, dimension must be greater than 1, got $(d)"))
    length(δ) == d || throw(ArgumentError("to use weighted FT, weights δ must be of length $(d), got $(length(δ))"))
    all(δ.>0) || throw(ArgumentError("to use weighted FT, weights δ must be positive"))
    all(isfinite.(δ)) || throw(ArgumentError("to use weighted FT, weight must be finite, got $(δ)"))
    return FT{d,Vector{T}}(δ)
end

filter(idxset::FT,sz) = itr -> all(idxset.δ.*itr.<=sz) 

"""
    TD(d, δ=ones(d))

A total degree index set in `d` dimenions with weights `δ`.
"""
struct TD{d, W<:AbstractVector} <: IndexSet{d}
    δ::W
end

function TD(d::Integer; δ::Vector{T}=ones(d)) where {T<:Real}
    d > 1 || throw(BoundsError("to use TD, dimension must be greater than 1, got $(d)"))
    length(δ) == d || throw(ArgumentError("to use weighted TD, weights δ must be of length $(d), got $(length(δ))"))
    all(δ.>0) || throw(ArgumentError("to use weighted TD, weights δ must be positive"))
    all(isfinite.(δ)) || throw(ArgumentError("to use weighted TD, weight must be finite, got $(δ)"))
    return TD{d,Vector{T}}(δ)
end

filter(idxset::TD,sz) = itr -> sum(idxset.δ.*itr) <= sz

"""
    HC(d, δ=ones(d))

A hyperbolic cross index set in `d` dimenions with weights `δ`.
"""
struct HC{d, W<:AbstractVector} <: IndexSet{d}
    δ::W
end

function HC(d::Integer; δ::Vector{T}=ones(d)) where {T<:Real}
    d > 1 || throw(BoundsError("to use HC, dimension must be greater than 1, got $(d)"))
    length(δ) == d || throw(ArgumentError("to use weighted HC, weights δ must be of length $(d), got $(length(δ))"))
    all(δ.>0) || throw(ArgumentError("to use weighted HC, weights δ must be positive"))
    all(isfinite.(δ)) || throw(ArgumentError("to use weighted HC, weight must be finite, got $(δ)"))
    return HC{d,Vector{T}}(δ)
end

filter(idxset::HC,sz) = itr -> prod(idxset.δ.*itr.+1) <= sz

"""
    AD(d)

An adaptive index set in `d` dimenions.
"""
struct AD{d} <: IndexSet{d} end

AD(d::Integer) = d <= 1 ? throw(BoundsError("to use AD, dimension must be greater than 1, got $(d)")) : AD{d}()

# Multigrid wrapper type
struct MG{d,I} <: IndexSet{d}
    idxset::I
end

MG(idxset::I) where {I<:IndexSet} = MG{ndims(idxset),I}(idxset)

# utilities
ndims(::IndexSet{d}) where {d} = d

function get_index_set(idxset::IndexSet{d},sz::N) where {d,N<:Integer}
    sz >= 0 || throw(ArgumentError("index set size parameter cannot be negative"))
    tensor_grid = Base.product(UnitRange.(0,ntuple(i->sz+1,d))...)
    filtered_grid = Base.Iterators.filter(filter(idxset,sz),tensor_grid)
    collect(filtered_grid)
end

get_index_set(idxset::MG,sz::N) where {N<:Integer} = get_index_set(idxset.idxset,sz)

function is_valid_index_set(idxset::Vector{T} where {T<:Index{d}}) where {d}
    for idx in idxset
        for i in 1:d
            check_index = Base.setindex(idx,idx[i]-1,i)
            ( in(check_index,idxset) || any(check_index.<0) ) || return false
        end
    end
    return true
end

function is_admissable(idxset::Set{Index{d}}, index::Index{d})  where {d}
    if in(index,idxset)
        return false
    else
        for k in 1:d
            new_index = index .- unit(k,d)
            if !( in(new_index,idxset) || new_index[k] < 0 )
                return false
            end
        end
        return true
    end
end

# output formatting
for i in ["SL" "ML" "TD" "HC" "FT" "AD"]
    ex = :( show(io::IO, ::$(Symbol(i))) = print(io, $(i)) )
    eval(ex)
end

show(io::IO, idxset::MG) = print(io, string("multigrid ", idxset.idxset))
