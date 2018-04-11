## index_sets.jl : representation of index sets

## Indexset ##
"""
IndexSet{d}

Supertype for an index set in `d` dimensions.
"""
abstract type IndexSet{d} end

"""
SL{d}

Representation of a single-level index set.
"""
struct SL{d} <: IndexSet{d}
end

"""
ML{d}

Representatrion of a multilevel index set.
"""
struct ML{d} <: IndexSet{d}
end

"""
FT{d, W<:AbstractVector}

Representation of a full tensor index set in `d` dimenions with weights `W`.
"""
struct FT{d, W<:AbstractVector} <: IndexSet{d}
    weights::W
end

"""
TD{d, W<:AbstractVector}

Representation of a total degree index set in `d` dimenions with weights `W`.
"""
struct TD{d, W<:AbstractVector} <: IndexSet{d}
    weights::W
end

"""
HC{d, W<:AbstractVector}

Representation of a hyperbolic cross index set in `d` dimenions with weights `W`.
"""
struct HC{d, W<:AbstractVector} <: IndexSet{d}
    weights::W
end

"""
AD{d}

Representation of an adaptove index set in `d` dimenions.
"""
struct AD{d} <: IndexSet{d}
end

# utilities
ndims(::IndexSet{d}) where {d} = d

# constructors for index sets
SL() = SL{1}()

ML() = ML{1}()

FT(d::N where N<:Integer) = d <= 1 ? throw(BoundsError("can only use FT with d>1")) : FT{d,Vector{Float64}}(ones(Float64,d))

function FT(weights::Vector{T}) where T<:AbstractFloat
    if ( length(weights) <= 1 ) || !all(weights .> 0)
        throw(ArgumentError("incorrect weights for type FT, must be of length d>1 and positive!"))
    end
    return FT{length(weights),Vector{T}}(weights)
end

TD(d::N where N<:Integer) = d <= 1 ? throw(BoundsError("can only use TD with d>1")) : TD{d,Vector{Float64}}(ones(Float64,d))

function TD(weights::Vector{T}) where T<:AbstractFloat
    if ( length(weights) <= 1 ) || !all(weights .> 0)
        throw(ArgumentError("incorrect weights for type TD, must be of length d>1 and positive!"))
    end
    return TD{length(weights),Vector{T}}(weights)
end

HC(d::N) where N<:Integer = d <= 1 ? throw(BoundsError("can only use HC with d>1")) : HC{d,Vector{Float64}}(ones(Float64,d))

function HC(weights::Vector{T}) where T<:AbstractFloat
    if ( length(weights) <= 1 ) || !all(weights .> 0)
        throw(ArgumentError("incorrect weights for type HC, must be of length d>1 and positive!"))
    end
    return HC{length(weights),Vector{T}}(weights)
end

AD(d::N where N<:Integer) = d <= 1 ? throw(BoundsError("can only use AD with d>1")) : AD{d}()

# output formatting
show(io::IO, sl::SL) = print(io, "Single level index set")

show(io::IO, ml::ML) = print(io, "Multilevel index set")

show(io::IO, ft::FT) = print(io, "$(ndims(ft))-dimensional index set of type FT with weights $(ft.weights)")

show(io::IO, td::TD) = print(io, "$(ndims(td))-dimensional index set of type TD with weights $(td.weights)")

show(io::IO, hc::HC) = print(io, "$(ndims(hc))-dimensional index set of type HC with weights $(hc.weights)")

show(io::IO, ad::AD) = print(io, "$(ndims(ad))-dimensional adaptive index set")

## Main methods for working with indices and index sets ##

# difference operator
function difference(i::Index{d}, current::Index{d}) where {d}
    p = 1
    returnvalue = zero(i)
    while i != -1
        i[p] -= 1
        if (i[p] < current[p]-1) || (!isvalid(i))
            if p == ndims(i)
                returnvalue = -one(i)
                break
            end
            i[p] = current[p]
            p = p + 1
        else
            returnvalue = i
            break
        end
    end
    return returnvalue
end

# drop algorithm for the iterative enumeration of all indices of given kind
# see "M. Holz, Sparse Grid Quadrature in High Dimensions with Applications, Springer, 2010"
function drop(i::Index{d}, L::N, indexset::IndexSet{d}) where {d,N<:Integer}
    p = 1
    returnvalue = zero(i)
    while i != -1 # at most two iterations
        i[p] += 1
        if check_drop(p,i,L,indexset)
            if p == ndims(i)
                returnvalue = -one(i)
                break
            end
            i[p] = 0
            p = p + 1
        else
            returnvalue = i
            break
        end
    end
    return returnvalue
end

check_drop(p::N, i::Index{d}, L::N, indexset::SL{d}) where {d,N<:Integer} = i[1] > L

check_drop(p::N, i::Index{d}, L::N, indexset::ML{d}) where {d,N<:Integer} = i[1] > L

check_drop(p::N, i::Index{d}, L::N, indexset::FT{d}) where {d,N<:Integer} = (i*indexset.weights)[p] > L

check_drop(p::N, i::Index{d}, L::N, indexset::TD{d}) where {d,N<:Integer} = sum(i*indexset.weights) > L

check_drop(p::N, i::Index{d}, L::N, indexset::HC{d}) where {d,N<:Integer} = prod(max(1,i*indexset.weights)) > L

# return index set of given kind for certain parameter L
"""
get_index_set(indexset, L)

Get all indices from a certain index set type `indexset`. The parameter `L` governs the size of the indexset.

# Examples

```jldoctest
julia> get_index_set(TD(2), 3) # get the indices of a 2-dimensional TD index set
Set(MultilevelEstimators.Index{2,Array{Int64,1}}[index [0, 3], index [3, 0], index [0, 1], index [1, 1], index [0, 0], index [2, 0], index [2, 1], index [1, 0], index [0, 2], index [1, 2]])
```
"""
function get_index_set(indexset::IndexSet{d}, L::N) where {d,N<:Integer}
    indices = Set{Index{d,Vector{N}}}()
    if L ≥ 0
        i = zero_index(N,d)
        while i != -1
            push!(indices,copy(i))
            i = drop(i,L,indexset)
        end
    end
    return indices
end

# returns the "boundary" indices of the given index set
"""
get_boundary(indices)

Get the boundary of the indices contained in the set `indices`. An index is part of the boundary if it is not dominated by another index in its maximum dimension.
# Examples

```jldoctest
julia> idcs = get_index_set(TD(2), 3); # get the indices of a 2-dimensional TD index set

julia> get_boundary(idcs)
Set(MultilevelEstimators.Index{2,Array{Int64,1}}[index [1, 2], index [3, 0], index [2, 1], index [0, 3]])

```
"""
function get_boundary(indices::Set{I}) where I<:Index
    boundary = Set{I}()
    for index in indices # could be made more efficient for TD and FT types, but uses same rules now for all types
        index_ = copy(index)
        maxs = indmax(index_) # maximum direction(s)
        isboundary = true
        for m in 1:length(maxs)
            index_[maxs[m]] += 1
            if in(index_,indices)
                isboundary = false
            end
            index_[maxs[m]] -= 1
        end
        if isboundary
            push!(boundary,index_)
        end
    end
    return boundary::Set{I}
end

# admissability for index sets of AD type
"""
is_admissable(set, index)

Check if the index `index` is admissable in the index set `set`.
```jldoctest
julia> idcs = get_index_set(TD(2), 3); # get the indices of a 2-dimensional TD index set

julia> is_admissable(idcs, Index(0,4))
true

julia> is_admissable(idcs, Index(0,0))
false

julia> is_admissable(idcs, Index(4,4))
false

```
"""
function is_admissable(set::Set{Index{d,V}} where V, i::Index{d}) where {d}
    ad = true
    if in(i,set)
        ad = false
    else
        for p in 1:d
            i[p] -= 1
            if !( in(i,set) || i[p] < 0 )
                ad = false
            end
            i[p] += 1
        end
    end
    return ad
end

# sort index sets (for nicer printing)
function sort(index_set::Set{I} where I<:Index)
    sorted_index_set = Index[]
    for index in index_set
        if isempty(sorted_index_set)
            push!(sorted_index_set,index)
        else
            i = 1
            while ( i <= length(sorted_index_set) ) && ( index > sorted_index_set[i] )
                i += 1
            end
            insert!(sorted_index_set,i,index)
        end
    end
    return sorted_index_set
end

# pretty print index set
"""
pretty_print(indexset)

Format the contents of the index set `indexset`.

```jldoctest
julia> idcs = get_index_set(TD(2), 3);

julia> print(pretty_print(idcs))
indexset with 10 indices in 2 dimensions:
**********
*  0  0  *
*  0  1  *
*  1  0  *
*  0  2  *
*  1  1  *
*  2  0  *
*  0  3  *
*  1  2  *
*  2  1  *
*  3  0  *
**********

```
"""
function pretty_print(indexset::Set{Index{d,V}} where V) where {d}
    str = "indexset with $(length(indexset)) indices in $(d) dimensions:\n"
    str *= repeat("*",4+3*d)*"\n"
    for idx in sort(indexset)
        str *= "*  "
        for i in 1:d
            str *= "$(idx[i])  "[1:3]
        end
        str *= "*\n"
    end
    str *= repeat("*",4+3*d)*"\n"
    return str
end