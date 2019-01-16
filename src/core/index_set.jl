## index_set.jl : representation of an index set
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
    SL()

Return a single-level index set.

# Examples
```jldoctest
julia> SL()
SL
```
See also: [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct SL{d} <: AbstractIndexSet{d} end

SL() = SL{1}()

filter(::SL, sz) = itr -> sum(itr) == sz

## ML ##
"""
    ML()

Return a multi-level index set.

# Examples
```jldoctest
julia> ML()
ML
```
See also: [`SL`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct ML{d} <: AbstractIndexSet{d} end

ML() = ML{1}()

filter(::ML, sz) = itr -> sum(itr) ≤ sz

## FT ##
"""
    FT(d::Integer)
    FT(δ::Real...)
    FT(δ::AbstractVector)
    FT(δ::NTuple)

Return a full tensor index set in `d` dimenions with optional weights `δ`. Default weights are all 1's.

# Examples
```jldoctest
julia> FT(2)
FT{2}

julia> FT([2, 1, 1])
FT{3}

julia> print(FT(2), 4)
  ◼ ◼ ◼ ◼ ◼
  ◼ ◼ ◼ ◼ ◼
  ◼ ◼ ◼ ◼ ◼
  ◼ ◼ ◼ ◼ ◼
  ◼ ◼ ◼ ◼ ◼
```
See also: [`SL`](@ref), [`ML`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct FT{d, T} <: AbstractIndexSet{d}
    δ::T
end

filter(idxset::FT, sz) = itr -> all(itr.I ./ idxset.δ .≤ sz) 

## TD ##
"""
    TD(d::Integer)
    TD(δ::Real...)
    TD(δ::AbstractVector)
    TD(δ::NTuple)

# Examples
```jldoctest
julia> TD(2)
TD{2}

julia> TD([2, 1, 1])
TD{3}

julia> print(TD(2), 4)
  ◼
  ◼ ◼
  ◼ ◼ ◼
  ◼ ◼ ◼ ◼
  ◼ ◼ ◼ ◼ ◼
```
Return a total degree index set in `d` dimenions with optional weights `δ`. Default weights are all 1's.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct TD{d, T} <: AbstractIndexSet{d}
    δ::T
end

filter(idxset::TD, sz) = itr -> sum(itr.I ./ idxset.δ) ≤ sz

## HC ##
"""
    HC(d::Integer)
    HC(δ::Real...)
    HC(δ::AbstractVector)
    HC(δ::NTuple)

# Examples
```jldoctest
julia> HC(2)
HC{2}

julia> HC([2, 1, 1])
HC{3}

julia> print(HC(2), 4)
  ◼
  ◼
  ◼
  ◼ ◼
  ◼ ◼ ◼ ◼ ◼
```
Return a hyperbolic cross index set in `d` dimenions with optional weights `δ`. Default weights are all 1's.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct HC{d, T} <: AbstractIndexSet{d}
    δ::T
end

filter(idxset::HC, sz) = itr -> prod(itr.I ./ idxset.δ .+ 1) - 1 ≤ sz

## ZC ##
"""
    ZC(d::Integer)
    ZC(δ::Real...)
    ZC(δ::AbstractVector)
    ZC(δ::NTuple)

# Examples
```jldoctest
julia> ZC(2)
ZC{2}

julia> ZC([2, 1, 1])
ZC{3}

julia> print(ZC(2), 4)
  ◼ ◼ 
  ◼ ◼ 
  ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ ◼ 
  ◼ ◼ ◼ ◼ ◼ 
```
Return a Zaremba cross index set in `d` dimenions with optional weights `δ`. Default weights are all 1's.

See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`AD`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct ZC{d, T} <: AbstractIndexSet{d}
    δ::T
end

filter(idxset::ZC, sz) = itr -> sz == 0 ? all(itr.I .== 0) : prod(max.(1, itr.I ./ idxset.δ)) ≤ sz

## constructors ##
for T in ["FT" "TD" "HC" "ZC"]
    @eval begin
        $(Symbol(T))(d::Integer) = $(Symbol(T))(fill(1, d))

        $(Symbol(T))(δ::Real...) = $(Symbol(T))(promote(δ...))
        
        $(Symbol(T))(δ::AbstractVector{<:Real}) = $(Symbol(T))(Tuple(δ))
        
        $(Symbol(T))(δ::NTuple{d, <:Real} where d) = check_args(δ, $(T)) && $(Symbol(T)){length(δ), typeof(δ)}(δ ./ maximum(δ))
    end
end

## AD ##
"""
    AD(d::Integer)

Return an adaptive index set in `d` dimenions.

# Examples
```jldoctest
julia> AD(2)
AD{2}
```
See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`U`](@ref), [`MG`](@ref)
"""
struct AD{d} <: AbstractIndexSet{d} end

AD(d::Integer) = check_larger_than(AD, d, "dimension", 1) && AD{d}()

## U ##
"""
    U(d::Integer)

Return an unbiased  multi-level or multi-index index set in `d` dimensions.

# Examples
```jldoctest
julia> U(2)
U{2}
```
See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`MG`](@ref)
"""
struct U{d} <: AbstractIndexSet{d} end

U(d::Integer) = check_larger_than(U, d, "dimension", 0) && U{d}()

## MG ##
"""
    MG(I::U)

Return a Multigrid wrapper for the index set `I` of type `U`.

# Examples
```jldoctest
julia> MG(U(1))
MG{U{1}}
```
See also: [`SL`](@ref), [`ML`](@ref), [`FT`](@ref), [`TD`](@ref), [`HC`](@ref), [`ZC`](@ref), [`AD`](@ref), [`U`](@ref)
"""
struct MG{d, I} <: AbstractIndexSet{d}
    idxset::I
end

MG(idxset::I) where {I<:U{d}} where d = MG{d, I}(idxset)

## utilities ##
AbstractMI = Union{TD, FT, HC, ZC, AD, U}

AbstractML = Union{ML, U{1}}

ndims(::AbstractIndexSet{d}) where {d} = d

length(idxset::AbstractIndexSet, sz::Integer) = length(collect(get_index_set(idxset, sz)))

eltype(::Type{<:AbstractIndexSet{d}}) where d = Index{d}

"""
    get_index_set(idxset::AbstractIndexSet, sz::Integer)

Return an iterator over all indices in the index set `idxset` for a given size parameter `sz`.

# Examples
```jldoctest
julia> collect(get_index_set(TD(2), 2))
6-element Array{Tuple{Int64,Int64},1}:
(0, 0)
(1, 0)
(2, 0)
(0, 1)
(1, 1)
(0, 2)
```
"""
get_index_set(idxset::AbstractIndexSet{d}, sz::Integer) where d = Base.Iterators.filter(filter(idxset, sz), CartesianIndices(ntuple(i -> 0:sz + 1, d)))

get_index_set(idxset::MG, sz::Integer) = get_index_set(idxset.idxset, sz)

get_index_set(idxset::AD, sz::Integer) = get_index_set_not_implemented("AD") 

get_index_set(idxset::U, sz::Integer) = get_index_set_not_implemented("U") 

get_index_set_not_implemented(name) = throw(ArgumentError(string("get_index_set not implemented for index sets of type ", name)))

function is_admissable(idxset::Set{Index{d}}, index::Index{d})  where d
    if index ∈ idxset
        return false
    else
        for i in 1:d
            bw_neighbour = index - CartesianIndex(ntuple(j -> j == i, d)) # backward neighbour
            if !( bw_neighbour ∈ idxset || bw_neighbour[i] < 0 )
                return false
            end
        end
        return true
    end
end

## output formatting ##
"""
    print(idxset::AbstractIndexSet{2}, sz::Integer)

Print the indices in the index set `idxset` for a given size parameter `sz`. Only implemented for two-dimensional index sets.

# Examples
```jldoctest
julia> collect(get_index_set(TD(2), 3))
  ◼
  ◼ ◼
  ◼ ◼ ◼
  ◼ ◼ ◼ ◼
```
"""
print(idxset::AbstractIndexSet, sz::Integer) = print(collect(get_index_set(idxset, sz)))

function print(idxset::Vector{<:Index{2}})
	if !isempty(idxset)
		char = "\u25FC"
		n = maximum(idxset).I
		R = CartesianIndices(ntuple(i -> 0:n[i], 2))
		A = map(i -> i ∈ idxset, R)
		str = ""
		for j in n[2]+1:-1:1
			str = string(str, "  ")
			for i in 1:n[1]+1
				str = A[i,j] ? string(str, char, " ") : str
			end
			str = string(str,"\n")
		end
		print(str)
	end
end

print(::Vector{<:Index}) = throw(ArgumentError("print only available for index sets with d = 2"))

show(io::IO, idxset::AbstractIndexSet) = print(io, shortname(idxset))

for T in ["SL" "ML"]
    @eval shortname(::$(Symbol(T))) = $(T)
end

for T in ["TD" "FT" "HC" "ZC" "AD" "U"]
    @eval shortname(idxset::$(Symbol(T))) = string($(T), "{", ndims(idxset), "}")
end

shortname(idxset::MG) = string("MG{", shortname(idxset.idxset), "}")

## input checking ##
function check_args(δ::NTuple{d, <:Real}, name) where d
    check_larger_than(name, d, "dimension", 1)
    check_larger_than(string("weighted ", name), δ, "weights δ", 0)
    check_smaller_than_or_equal_to(string("weighted ", name), δ, "weights δ", 1)
end
