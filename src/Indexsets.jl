# representation of an index
type Index{d,V<:AbstractVector}
  indices::V
end

# utilities
ndims{d}(index::Index{d}) = d

# constructor for indices
Index{N<:Integer}(index::Vector{N}) = Index{length(index),Vector{N}}(index)
Index{N<:Integer}(i::N...) = Index{length(i),Vector{N}}([i...]) # slurping and splatting

# methods
isvalid{I<:Index}(index::I) = all(index.indices .≥ 0)::Bool

zero{d,V}(::Type{Index{d,V}}) = Index(zeros(eltype(V),d))::Index{d,V}
zero{d,V}(index::Index{d,V}) = Index(zeros(eltype(V),d))::Index{d,V}

one{d,V}(::Type{Index{d,V}}) = Index(ones(eltype(V),d))::Index{d,V}
one{d,V}(index::Index{d,V}) = Index(ones(eltype(V),d))::Index{d,V}

*{d,V}(index1::Index{d,V}, index2::Index{d,V}) = Index(index1.indices.*index2.indices)::Index{d,V} # can only multiply indices of same length
*{N<:Integer,I<:Index}(s::N, index::I) = Index(s*index.indices)::I
*{I<:Index,N<:Integer}(index::I, s::N) = Index(s*index.indices)::I
.*{T,I<:Index}(s::Vector{T}, index::I) = (s.*index.indices)::Vector{T}
.*{I<:Index,T}(index::I, s::Vector{T}) = (s.*index.indices)::Vector{T}

=={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices == index2.indices )::Bool
=={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices == (s*one(index1)).indices )::Bool
.=={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices .== (s*one(index1)).indices )

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
prod{I<:Index}(index::I) = prod(max(1,index.indices))
diff{I<:Index}(index1::I, index2::I) = count(!,index1.indices.==index2.indices) # returns number of indices that is different
length{I<:Index}(index::I) = length(index.indices)

copy{I}(index::I) = Index(copy(index.indices))
hash(index::Index, h::UInt) = hash(index.indices, hash(:Index, h)) # needed for looking things up in a dict

function show(io::IO, index::Index)
  d = ndims(index)
  print(io, "$d-dimensional index with values [")
  for i = 1:d-1
    print(io, @sprintf("%i,",index[i]))
  end
  print(io, @sprintf("%i]",index[d]))
end

# enum with index set types
@enum Kind ML FT TD HC AD

# representation of an index set
type IndexSet{K, d, W<:AbstractVector}
  weights::W
end

# utilities
ndims{K,d}(::IndexSet{K,d}) = d
kind{K,d}(::IndexSet{K,d}) = K

# constructor for index sets
IndexSet{T<:AbstractFloat}(K::Kind, weights::Vector{T}) = IndexSet{K,length(weights),Vector{T}}(weights)

# create index set
function createIndexSet{N<:Integer,T<:AbstractFloat}(K::Kind, d::N; weights::Vector{T} = ones(Float64,d))
  d > 0 || error("dimension of index set cannot be negative or zero!")
  K == ML ? ( d == 1 || error("can only perform MLMC when d=1") ) : ( d != 1 || error("cannot perform MIMC simulation with d=1, choose ML or increase dimension!") )
  K == ML ? weights == ones(T,d) || error("cannot assign weights to index set of type ML!") : []
  ( ( length(weights) == d ) && all(weights .> 0) ) || error("incorrect weights specified!")

  return IndexSet(K, weights)
end

# methods
isValid{I<:IndexSet}(indexset::I) = ( typeof(ndims(indexset)) <: Integer && ndims(indexset) > 0 && length(indexset.weights) == ndims(indexset) )

function show(io::IO, indexset::IndexSet)
  print(io, "$(ndims(indexset))-dimensional index set of type $(kind(indexset)) with weights $(indexset.weights)")
end

#
# Main methods for working with indices and index sets
#

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
function drop{d,V,N<:Integer,K}(i::Index{d,V}, L::N, indexset::IndexSet{K,d})
  p = 1
  returnvalue = zero(i)
  while i != -1 # at most two iterations
    i[p] += 1
     if ( K == ML && all(i[:].>L) ) ||
        ( K == FT && (i.*indexset.weights)[p] > L ) ||
        ( K == TD && sum(i.*indexset.weights) > L ) ||
        ( K == HC && prod(max(1,i.*indexset.weights)) > L )     
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

# return index set of given kind for certain parameter L
function getIndexSet{K,d,N<:Integer}(indexset::IndexSet{K,d}, L::N)
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
    maxs = find(index.==maximum(index)) # maximum direction(s)
    isboundary = true
    for m in 1:length(maxs)
      index[maxs[m]] += 1
      if in(index,indices)
        isboundary = false
      end
      index[maxs[m]] -= 1
    end
    if isboundary
      push!(boundary,index)
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