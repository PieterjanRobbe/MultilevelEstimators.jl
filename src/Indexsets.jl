# representation of an index
type Index{d,V<:AbstractVector}
  indices::V
end

# utilities
ndims{d}(index::Index{d}) = d

# constructor for indices
Index{N<:Integer}(index::Vector{N}) = Index{length(index),Vector{N}}(index)
Index{N<:Integer}(i::N...) = Index{length(i),Vector{N}}([i...]) # slurping and splatting

# create Index form String
function Index{S<:AbstractString}(index_string::S)
  index_array = split(index_string,',')
  if length(index_array) == 1
    idx = tryparse(Int64,index_array[1][2:end-1])
    !isnull(idx) || error("could not parse $(index_array[1][2:end-1]) into Int64")
    return Index(parse(Int64,index_array[1][2:end-1]))
  else
    idx = zeros(Int64,length(index_array))
    tmp = tryparse(Int64,index_array[1][2:end])
    if !isnull(tmp)
      idx[1] = get(tmp)
    else
      error("could not parse $(index_array[1][2:end]) into Int64")
    end
    tmp = tryparse(Int64,index_array[end][1:end-1])
    if !isnull(tmp)
      idx[end] = get(tmp)
    else
      error("could not parse $(index_array[end][1:end-1]) into Int64")
    end
    for i = 2:length(index_array)-1
      tmp = tryparse(Int64,index_array[i])
      if !isnull(tmp)
        idx[i] = get(tmp)
      else
        error("could not parse $(index_array[i]) into Int64")
      end
    end
    return Index(idx)
  end
end

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

# representation of an index set
abstract IndexSet{d}

type SL{d} <: IndexSet{d}
end

type ML{d} <: IndexSet{d}
end

type FT{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

type TD{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

type HC{d, W<:AbstractVector} <: IndexSet{d}
  weights::W
end

type AD{d} <: IndexSet{d}
end

# utilities
ndims{d}(::IndexSet{d}) = d

# constructor for index sets
function SL()
  return SL{1}()
end

function ML()
  return ML{1}()
end

function FT{N<:Int}(d::N)
  if d <= 1
    error("can only use FT with d>1!")
  end
  return FT{d,Vector{Float64}}(ones(Float64,d))
end

function FT{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
    error("incorrect weights for type FT, must be of length d>1 and positive!")
  end
  return FT{length(weights),Vector{T}}(weights)
end

function TD{N<:Int}(d::N)
  if d <= 1
    error("can only use TD with d>1!")
  end
  return TD{d,Vector{Float64}}(ones(Float64,d))
end

function TD{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
    error("incorrect weights for type TD, must be of length d>1 and positive!")
  end
  return TD{length(weights),Vector{T}}(weights)
end

HC{T<:AbstractFloat}(weights::Vector{T}) = HC{length(weights),Vector{T}}(weights)
function HC{N<:Int}(d::N)
  if d <= 1
    error("can only use HC with d>1!")
  end
  return HC{d,Vector{Float64}}(ones(Float64,d))
end

function HC{T<:AbstractFloat}(weights::Vector{T})
  if ( length(weights) <= 1 ) || !all(weights .> 0)
    error("incorrect weights for type HC, must be of length d>1 and positive!")
  end
  return HC{length(weights),Vector{T}}(weights)
end

function AD{N<:Int}(d::N)
  if d<=1
    error("can only use AD with d>1!")
  end
  return AD{d}()
end

# methods
isValid{I<:IndexSet}(indexset::I) = ( typeof(ndims(indexset)) <: Integer && ndims(indexset) > 0 && length(indexset.weights) == ndims(indexset) )

function show(io::IO, sl::SL)
  print(io, "Single level index set")
end

function show(io::IO, ml::ML)
  print(io, "Multilevel index set")
end

function show(io::IO, ft::FT)
  print(io, "$(ndims(ft))-dimensional index set of type FT with weights $(ft.weights)")
end

function show(io::IO, td::TD)
  print(io, "$(ndims(td))-dimensional index set of type TD with weights $(td.weights)")
end

function show(io::IO, hc::HC)
  print(io, "$(ndims(hc))-dimensional index set of type HC with weights $(hc.weights)")
end

function show(io::IO, ad::AD)
  print(io, "$(ndims(ad))-dimensional adaptive index set")
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

function check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::SL{d})
  return i[1] > L
end

function check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::ML{d})
  return i[1] > L
end

function check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::FT{d})
  return (i.*indexset.weights)[p] > L
end

function check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::TD{d})
  return sum(i.*indexset.weights) > L
end

function check_drop{d,N<:Integer}(p::N, i::Index{d}, L::N, indexset::HC{d})
  return prod(max(1,i.*indexset.weights)) > L
end

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