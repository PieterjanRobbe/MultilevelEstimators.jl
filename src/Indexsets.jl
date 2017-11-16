## Indexsets.jl : representation of indices and index sets

## Index ##
struct Index{d,V<:AbstractVector}
  indices::V
end

# utilities
ndims{d}(index::Index{d}) = d
isvalid(index::Index) = all(index.indices.>=0)

# constructor for indices
Index{N<:Integer}(index::Vector{N}) = Index{length(index),Vector{N}}(index)
Index{N<:Integer}(i::N...) = Index{length(i),Vector{N}}([i...]) # slurping and splatting

# create Index form String
function Index{S<:AbstractString}(index_string::S)
  index_array = split(index_string,',')
  if length(index_array) == 1
    idx = tryparse(Int64,index_array[1][2:end-1])
	!isnull(idx) || throw(ArgumentError("could not parse $(index_array[1][2:end-1]) into Int64"))
    return Index(parse(Int64,index_array[1][2:end-1]))
  else
    idx = zeros(Int64,length(index_array))
    tmp = tryparse(Int64,index_array[1][2:end])
    if !isnull(tmp)
      idx[1] = get(tmp)
    else
	  throw(ArgumentError("could not parse $(index_array[1][2:end]) into Int64"))
    end
    tmp = tryparse(Int64,index_array[end][1:end-1])
    if !isnull(tmp)
      idx[end] = get(tmp)
    else
	  throw(ArgumentError("could not parse $(index_array[end][1:end-1]) into Int64"))
    end
    for i = 2:length(index_array)-1
      tmp = tryparse(Int64,index_array[i])
      if !isnull(tmp)
        idx[i] = get(tmp)
      else
		throw(ArgumentError("could not parse $(index_array[i]) into Int64"))
      end
    end
    return Index(idx)
  end
end

# methods
zero{d,V}(::Type{Index{d,V}}) = Index(zeros(eltype(V),d))::Index{d,V}
zero{d,V}(index::Index{d,V}) = Index(zeros(eltype(V),d))::Index{d,V}

one{d,V}(::Type{Index{d,V}}) = Index(ones(eltype(V),d))::Index{d,V}
one{d,V}(index::Index{d,V}) = Index(ones(eltype(V),d))::Index{d,V}

*{d,V}(index1::Index{d,V}, index2::Index{d,V}) = Index(index1.indices.*index2.indices)::Index{d,V} # can only multiply indices of same length
*{d,V,T<:AbstractFloat}(index1::Index{d,V}, vec::Vector{T}) = index1.indices.*vec
*{d,V,T<:AbstractFloat}(vec::Vector{T}, index1::Index{d,V}) = index1.indices.*vec
*{N<:Integer,I<:Index}(s::N, index::I) = Index(s*index.indices)::I
*{I<:Index,N<:Integer}(index::I, s::N) = Index(s*index.indices)::I
*{N<:Number,I<:Index}(s::N, index::I) = s*index.indices
*{I<:Index,N<:Number}(index::I, s::N) = s*index.indices

=={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices == index2.indices )::Bool
=={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices == (s*one(index1)).indices )::Bool

function isless{d}(i1::Index{d},i2::Index{d})
  returnvalue = false
  if isless(sum(i1),sum(i2))
    returnvalue = true
  elseif isequal(i1,i2)
    returnvalue = false
  elseif isequal(sum(i1),sum(i2))
    idx = find(i1.indices .!= i2.indices)[1]
    returnvalue = ( isless(i1[idx],i2[idx]) ? true : false )
  end
  return returnvalue::Bool
end

+{d,V}(index1::Index{d,V}, index2::Index{d,V}) = Index(index1.indices+index2.indices)::Index{d,V}
-{d,V}(index1::Index{d,V}, index2::Index{d,V}) = Index(index1.indices-index2.indices)::Index{d,V}
+{I<:Index,N<:Integer}(index1::I, s::N) = Index(index1.indices+s)::I
-{I<:Index,N<:Integer}(index1::I, s::N) = Index(index1.indices-s)::I
-{I<:Index}(index::I) = Index(-index.indices)::I

getindex{I<:Index,N<:Integer}(index::I, i::N) = index.indices[i]
getindex{I<:Index}(index::I, ::Colon) = index.indices
setindex!{I<:Index}(index::I, value, key) = index.indices[key] = value
maximum{I<:Index}(index::I) = maximum(index.indices)
indmax{I<:Index}(index::I) = indmax(index.indices)
sum{I<:Index}(index::I) = sum(index.indices)
prod{I<:Index}(index::I) = prod(max.(1,index.indices))
diff{I<:Index}(index1::I, index2::I) = count(!,index1.indices.==index2.indices) # returns number of indices that is different
length{I<:Index}(index::I) = length(index.indices)

copy{I}(index::I) = Index(copy(index.indices))
hash(index::Index, h::UInt) = hash(index.indices, hash(:Index, h)) # needed for looking things up in a dict

# output formatting
function show(io::IO, index::Index)
  d = ndims(index)
  print(io, "$d-dimensional index with values [")
  for i = 1:d-1
    print(io, @sprintf("%i,",index[i]))
  end
  print(io, @sprintf("%i]",index[d]))
end

## Indexset ##
abstract type IndexSet{d} end

struct SL{d} <: IndexSet{d}
end

struct ML{d} <: IndexSet{d}
end

struct FT{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

struct TD{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

struct HC{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

struct AD{d} <: IndexSet{d}
end

# utilities
ndims{d}(::IndexSet{d}) = d

# constructors for index sets
SL() = SL{1}()

ML() = ML{1}()

FT{N<:Int}(d::N) = d <= 1 ? throw(BoundsError("can only use FT with d>1")) : FT{d,Vector{Float64}}(ones(Float64,d))

function FT{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
	throw(ArgumentError("incorrect weights for type FT, must be of length d>1 and positive!"))
  end
  return FT{length(weights),Vector{T}}(weights)
end

TD{N<:Int}(d::N) = d <= 1 ? throw(BoundsError("can only use TD with d>1")) : TD{d,Vector{Float64}}(ones(Float64,d))

function TD{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
	throw(ArgumentError("incorrect weights for type TD, must be of length d>1 and positive!"))
  end
  return TD{length(weights),Vector{T}}(weights)
end

HC{N<:Int}(d::N) = d <= 1 ? throw(BoundsError("can only use HC with d>1")) : HC{d,Vector{Float64}}(ones(Float64,d))

function HC{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
	throw(ArgumentError("incorrect weights for type HC, must be of length d>1 and positive!"))
  end
  return HC{length(weights),Vector{T}}(weights)
end

AD{N<:Int}(d::N) = d <= 1 ? throw(BoundsError("can only use AD with d>1")) : AD{d}()

# output formatting
show(io::IO, sl::SL) = print(io, "Single level index set")

show(io::IO, ml::ML) = print(io, "Multilevel index set")

show(io::IO, ft::FT) = print(io, "$(ndims(ft))-dimensional index set of type FT with weights $(ft.weights)")

show(io::IO, td::TD) = print(io, "$(ndims(td))-dimensional index set of type TD with weights $(td.weights)")

show(io::IO, hc::HC) = print(io, "$(ndims(hc))-dimensional index set of type HC with weights $(hc.weights)")

show(io::IO, ad::AD) = print(io, "$(ndims(ad))-dimensional adaptive index set")

## Main methods for working with indices and index sets ##

# difference operator
function difference{d}(i::Index{d}, current::Index{d})
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
function drop{d,N<:Integer}(i::Index{d}, L::N, indexset::IndexSet{d})
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

check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::SL{d}) = i[1] > L

check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::ML{d}) = i[1] > L

check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::FT{d}) = (i*indexset.weights)[p] > L

check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::TD{d}) = sum(i*indexset.weights) > L

check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::HC{d}) = prod(max(1,i*indexset.weights)) > L

# return index set of given kind for certain parameter L
function getIndexSet{d,N<:Integer}(indexset::IndexSet{d}, L::N)
  indices = Set{Index{d,Vector{N}}}()
  if L ≥ 0
    i = Index(zeros(N,d))::Index{d,Vector{N}}
    while i != -1
      push!(indices,copy(i))
      i = drop(i,L,indexset)::Index{d,Vector{N}}
    end
  end
  return indices
end

# returns the "boundary" indices of the given index set
# we define an index as a boundary index when it is not dominated
# by another index in its maximum direction(s)
function getBoundary{d,V}(indices::Set{Index{d,V}}) # must have explicit d in signature
  boundary = Set{Index{d,V}}()
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
  return boundary::Set{Index{d,V}}
end

# admissability for index sets of AD type
function isAdmissable{d,V}(set::Set{Index{d,V}}, i::Index{d,V}) # must have explicit d in signature to ensure same length
  ad = true
  if in(i,set)
    ad = false
  else
    for p in 1:d
      i[p] -= 1
      if !( in(i,set) || i[p] < 0 ) || i[p] >= 15
        ad = false
      end
      i[p] += 1
    end
  end
  return ad
end

# sort index sets (for nicer printing)
function sort{I<:Index}(indexSet::Set{I})
  sortedIndexSet = Index[]
  for index in indexSet
    if isempty(sortedIndexSet)
      push!(sortedIndexSet,index)
    else
      i = 1
      while ( i <= length(sortedIndexSet) ) && ( index > sortedIndexSet[i] )
        i += 1
      end
      insert!(sortedIndexSet,i,index)
    end
  end
  return sortedIndexSet
end

# pretty print index set
function prettyprint{d,N}(indexset::Set{Index{d,Vector{N}}})
	println("indexset with $(length(indexset)) indices in $(d) dimensions:")
	println(repeat("*",4+3*d))
	for idx in sort(indexset)
		str = "*  "
		for i in 1:d
			str *= "$(idx[i])  "[1:3]
		end
		println(str*"*")
	end
	println(repeat("*",4+3*d))
end
