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

# utilities
ndims(::IndexSet{d}) where {d} = d

function get_index_set(idxset::IndexSet{d},sz::N) where {d,N<:Integer}
    sz >= 0 || throw(ArgumentError("index set size parameter cannot be negative"))
    tensor_grid = Base.product(range.(0,ntuple(i->sz+1,d))...)
    filtered_grid = Base.Iterators.filter(filter(idxset,sz),tensor_grid)
    collect(filtered_grid)
end

function is_valid_index_set(idxset::Vector{T} where {T<:Index{d}}) where {d}
    for idx in idxset
        for i in 1:d
            check_index = Base.setindex(idx,idx[i]-1,i)
            ( in(check_index,idxset) || any(check_index.<0) ) || return false
        end
    end
    return true
end

# output formatting
show(io::IO, ::SL) = print(io, "single-level index set")
show(io::IO, ::ML) = print(io, "multi-level index set")
show(io::IO, ft::FT) = print(io, "$(ndims(ft))-dimensional index set of type FT with weights $(ft.δ)")
show(io::IO, td::TD) = print(io, "$(ndims(td))-dimensional index set of type TD with weights $(td.δ)")
show(io::IO, hc::HC) = print(io, "$(ndims(hc))-dimensional index set of type HC with weights $(hc.δ)")
show(io::IO, ad::AD) = print(io, "$(ndims(ad))-dimensional adaptive index set")
