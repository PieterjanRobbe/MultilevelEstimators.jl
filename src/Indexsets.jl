# representation of an "index" in multiple dimensions
type Index{d,V<:AbstractVector}
  indices::V

  #Index(indices::AbstractVector{N}) = new(indices)
end
#Index(d,index::AbstractVector) = Index{d,eltype(index),typeof(index)}(index)

# constructor for indices
Index{N<:Integer}(index::Vector{N}) = Index{length(index),Vector{N}}(index)
Index{N<:Integer}(i::N...) = Index{length(i),Vector{N}}([i...]) # slurping and splatting

ndims{d}(index::Index{d}) = d

zero{d,V}(::Type{Index{d,V}}) = Index(zeros(eltype(V),d))
zero{d,V}(index::Index{d,V}) = Index(zeros(eltype(V),d))

one{d,V}(::Type{Index{d,V}}) = Index(ones(eltype(V),d))
one{d,V}(index::Index{d,V}) = Index(ones(eltype(V),d))

*{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices.*index2.indices) # can only multiply indices of same length
*{N<:Integer,I<:Index}(s::N, index::I) = Index(s*index.indices)
*{I<:Index,N<:Integer}(index::I, s::N) = Index(s*index.indices)
*{T,I<:Index}(s::Vector{T}, index::I) = s.*index.indices
*{I<:Index,T}(index::I, s::Vector{T}) = s.*index.indices

=={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices == index2.indices )
=={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices == (s*one(index1)).indices )
.=={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices .== (s*one(index1)).indices )

hash(index::Index, h::UInt) = hash(index.indices, hash(:Index, h))
isequal(index1::Index, index2::Index) = isequal(hash(index1),hash(index2))

!={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices != index2.indices )
!={I<:Index,N<:Integer}(index1::I, s::N) = ( index1.indices != (s*one(index1)).indices )

+{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices+index2.indices)
-{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices-index2.indices)
+{I<:Index,N<:Integer}(index1::I, s::N) = Index(index1.indices+s)
-{I<:Index,N<:Integer}(index1::I, s::N) = Index(index1.indices-s)
-{I<:Index}(index::I) = Index(-index.indices)

getindex{I<:Index,N<:Integer}(index::I, i::N) = index.indices[i]
getindex{I<:Index}(index::I, ::Colon) = index.indices
setindex!{I<:Index}(index::I, value, key) = index.indices[key] = value
maximum{I<:Index}(index::I) = maximum(index.indices)
indmax{I<:Index}(index::I) = indmax(index.indices)
sum{I<:Index}(index::I) = sum(index.indices)
prod{I<:Index}(index::I) = prod(max(1,index.indices))
diff{I<:Index}(index1::I, index2::I) = count(!,index1.indices.==index2.indices) # returns number of indices that is different

copy{I}(index::I) = Index(copy(index.indices))

function show(io::IO, index::Index)
  d = ndims(index)
  print(io, "$d-dimensional index with values [")
  for i = 1:d-1
    print(io, @sprintf("%i,",index[i]))
  end
  print(io, @sprintf("%i]",index[d]))
end

isvalid{I<:Index}(index::I) = all(index.indices .≥ 0)

# enum with index set types
@enum Kind ML FT TD HC AD

# representation of an index set
type IndexSet{d, W<:AbstractVector}
  kind::Kind
  weights::W
end

ndims{d}(index::Index{d}) = d

# create index set
function createIndexSet{N<:Int,T<:AbstractFloat}(kind::Kind, d::N; weights::Vector{T} = ones(Float64,d))
  d > 0 || error("dimension of index set cannot be negative or zero!")
  (d > 1) $ (kind == ML) || error("cannot perform MIMC simulation with d=1, choose ML or increase dimension!")
  kind == ML ? weights == ones(T,d) || error("cannot assign weights to index set of type ML!") : []
  ( ( length(weights) == d ) && all(weights .> 0) ) || error("incorrect weights specified!")

  return IndexSet{length(weights),Vector{T}}(kind, weights)::IndexSet{length(weights),Vector{T}}
end

# convenience methods for index sets
ndims{d}(indexSet::IndexSet{d}) = d

function show(io::IO, indexset::IndexSet)
  d = ndims(indexset)
  k = indexset.kind
  w = indexset.weights
  print(io, "$d-dimensional index set of type $k with weights $w")
end

# difference operator
function difference{N<:Integer,d}(p::N, i::Index{d}, current::Index{d})
  while i != -1
    i[p] -= 1
    if (i[p] < current[p]-1) || (!isvalid(i))
      if p == ndims(i)
        return -one(i)
      end
      i[p] = current[p]
      p = p + 1
    else
      p = 1
      return i
    end
  end
end

# drop algorithm for the iterative enumeration of all indices of given kind
# see "M. Holz, Sparse Grid Quadrature in High Dimensions with Applications, Springer, 2010"
function drop{N<:Integer,d,T<:AbstractFloat}(p::N, i::Index{d}, K::N, kind::Kind, weights::Vector{T})
  while i != -1 # at most two iterations
    i[p] += 1
    if ( kind == ML && all(i[:].>K) ) ||
      ( kind == FT && (i*weights)[p] > K ) ||
      ( kind == TD && sum(i*weights) > K ) ||
      ( kind == HC && prod(max(1,i*weights)) > K )     
      if p == ndims(i)
        return -one(i)
      end
      i[p] = 0
      p = p + 1
    else
      p = 1
      return i
    end
  end
end

# return index set of given kind for certain parameter K
function getIndexSet{S<:IndexSet,N<:Integer}(indexset::S, K::N)
  d = ndims(indexset)
  indices = Index{d,Vector{N}}[]
  if K ≥ 0
    i = Index(zeros(N,d))
    p = 1
    while i != -1
      push!(indices,copy(i))
      i = drop(p,i,K,indexset.kind,indexset.weights)
    end
  end
  return Set(indices)
end

# returns the "boundary" indices of the given index set
# we define an index as a boundary index when it is not dominated
# by another index in its maximum direction(s)
function getBoundary{S<:IndexSet,N<:Integer}(indexset::S, K::N)
  indices = getIndexSet(indexset,K)
  return getBoundary(indices)
end

# returns the "boundary" indices of an adaptive index set
function getBoundary{d,V}(indices::Set{Index{d,V}}) # must have explicit d in signature
  boundary = Set{Index{d,V}}()
  for index in indices
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
