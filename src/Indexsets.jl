# representation of an "index" in multiple dimensions
type Index{d}
  indices::Array{Int32,1}
end

# constructor for indices
Index{N<:Integer}(index::Array{N,1}) = Index{length(index)}(index)
Index{N<:Integer}(i::N...) = Index{length(i)}([i...]) # slurping and splatting

ndims{d}(index::Index{d}) = d

zero{d}(::Type{Index{d}}) = Index(zeros(Int32,d))
zero{d}(index::Index{d}) = Index(zeros(Int32,d))

one{d}(::Type{Index{d}}) = Index(ones(Int32,d))
one{d}(index::Index{d}) = Index(ones(Int32,d))

*{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices.*index2.indices)
*{d,N}(s::N, index::Index{d}) = Index(s*index.indices)
*{d,N}(index::Index{d}, s::N) = Index(s*index.indices)
*{d,T}(s::Vector{T}, index::Index{d}) = s.*index.indices
*{d,T}(index::Index{d}, s::Vector{T}) = s.*index.indices

=={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices == index2.indices )
=={d,N<:Integer}(index1::Index{d}, s::N) = ( index1.indices == (s*one(index1)).indices )
.=={d,N<:Integer}(index1::Index{d}, s::N) = ( index1.indices .== (s*one(index1)).indices )

hash(index::Index, h::UInt) = hash(index.indices, hash(:Index, h))
isequal(index1::Index, index2::Index) = isequal(hash(index1),hash(index2))

!={d}(index1::Index{d}, index2::Index{d}) = ( index1.indices != index2.indices )
!={d,N}(index1::Index{d}, s::N) = ( index1.indices != (s*one(index1)).indices )

+{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices+index2.indices)
-{d}(index1::Index{d}, index2::Index{d}) = Index(index1.indices-index2.indices)
+{d,N}(index1::Index{d}, s::N) = Index(index1.indices+s)
-{d,N}(index1::Index{d}, s::N) = Index(index1.indices-s)
-{d}(index::Index{d}) = Index(-index.indices)

getindex{d}(index::Index{d}, i::Integer) = index.indices[i]
getindex{d}(index::Index{d}, ::Colon) = index.indices
setindex!{d}(index::Index{d}, value, key) = index.indices[key] = value
maximum{d}(index::Index{d}) = maximum(index.indices)
indmax{d}(index::Index{d}) = indmax(index.indices)
sum{d}(index::Index{d}) = sum(index.indices)
prod{d}(index::Index{d}) = prod(max(1,index.indices))
diff{d}(index1::Index{d}, index2::Index{d}) = count(!,index1.indices.==index2.indices) # returns number of indices that is different

copy{d}(index::Index{d}) = Index(copy(index.indices))

function show(io::IO, index::Index)
  d = ndims(index)
  print(io, "$d-dimensional index with values [")
  for i = 1:d-1
    print(io, @sprintf("%i,",index[i]))
  end
  print(io, @sprintf("%i]",index[d]))
end

isvalid{d}(index::Index{d}) = all(index.indices .≥ 0)

# enum with index set types
@enum Kind ML FT TD HC AD

# representation of an index set
type IndexSet{d, T<:AbstractFloat}
  kind::Kind
  weights::Vector{T}
end

# create index set
function createIndexSet{N<:Int,T<:AbstractFloat}(kind::Kind, d::N; weights::Array{T,1} = ones(Float64,d))
  d > 0 || error("dimension of index set cannot be negative or zero!")
  (d > 1) $ (kind == ML) || error("cannot perform MIMC simulation with d=1, choose ML or increase dimension!")
  kind == ML ? weights == ones(T,d) || error("cannot assign weights to index set of type ML!") : []
  ( ( length(weights) == d ) && all(weights .> 0) ) || error("incorrect weights specified!")
  
  return IndexSet{d,T}(kind, weights)
end

# convenience methods for index sets
ndims{d,T}(indexSet::IndexSet{d,T}) = d

function show(io::IO, indexset::IndexSet)
  d = ndims(indexset)
  k = indexset.kind
  w = indexset.weights
  print(io, "$d-dimensional index set of type $k with weights $w")
end

# difference operator
function difference{d,N}(p::N, i::Index{d}, current::Index{d})
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
function drop{d,N,T}(p::N, i::Index{d}, K::N, kind::Kind, weights::Array{T,1})
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
function getIndexSet{d,N,T}(indexset::IndexSet{d,T}, K::N)
  indices = Index{d}[]
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
function getBoundary{d,N,T}(indexset::IndexSet{d,T}, K::N)
  indices = getIndexSet(indexset,K)
  return getBoundary(indices)
end

# returns the "boundary" indices of an adaptive index set
function getBoundary{d}(indices::Set{Index{d}})
  boundary = Set{Index{d}}()
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
  return boundary
end

# admissability for index sets of AD type
function isAdmissable{d}(set::Set{Index{d}}, i::Index{d})
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
