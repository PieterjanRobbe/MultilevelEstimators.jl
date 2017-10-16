## indices.jl : representation of indices

## Index ##
"""
Index{d,V<:AbstractVector}

Representation of a multi-index in `d` dimensions.

```jldoctest
julia> Index([2,1,3])
index [2, 1, 3]

julia> Index(0,0)
index [0, 0]

julia> Index("[2,3,1]")
index [2, 3, 1]

```
"""
struct Index{d,V<:AbstractVector}
    indices::V
end

# utilities
ndims(index::Index{d}) where {d} = d
isvalid(index::Index) = all(index.indices.>=0)

# constructor for indices
Index(index::Vector{N}) where {N<:Integer} = Index{length(index),Vector{N}}(index)
Index(i::N...) where {N<:Integer} = Index([i...]) # slurping and splatting

# create index form String
function Index(index_string::S) where {S<:AbstractString}
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

# zero index
zero(::Type{Index{d,V}}) where {d,V} = Index(zeros(eltype(V),d))::Index{d,V}
zero(index::Index{d,V}) where {d,V} = Index(zeros(eltype(V),d))::Index{d,V}
zero_index(d) = Index(zeros(Int64,d))::Index{d,Vector{Int64}}
zero_index(::Type{N},d) where N<:Integer = Index(zeros(N,d))::Index{d,Vector{N}}

# one index
one(::Type{Index{d,V}}) where {d,V} = Index(ones(eltype(V),d))::Index{d,V}
one(index::Index{d,V}) where {d,V} = Index(ones(eltype(V),d))::Index{d,V}

# multiplication
*(index1::Index{d,V}, index2::Index{d,V}) where {d,V} = Index(index1.indices.*index2.indices)::Index{d,V} 
*(index1::Index{d,V} where {d,V}, vec::Vector{T} where {T<:AbstractFloat}) = index1.indices.*vec
*(vec::Vector{T} where T<:AbstractFloat, index1::Index{d,V} where {d,V}) = index1.indices.*vec
*(s::N where N<:Integer, index::I) where I<:Index = Index(s*index.indices)::I
*(index::I, s::N where N<: Integer) where I<:Index = Index(s*index.indices)::I
*(s::N where N<:Number, index::I where I<:Index) = s*index.indices
*(index::I where I<:Index, s::N where N<:Number) = s*index.indices

# comparison
==(index1::I, index2::I) where I<:Index = ( index1.indices == index2.indices )::Bool
==(index1::I where I<:Index, s::N where N<:Integer) = ( index1.indices == (s*one(index1)).indices )::Bool

function isless(i1::I,i2::I) where {I<:Index}
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

# addition and substraction
+(index1::I, index2::I) where I<:Index = Index(index1.indices+index2.indices)::I
-(index1::I, index2::I) where I = Index(index1.indices-index2.indices)::I
+(index1::I, s::N where N<:Integer) where I<:Index = Index(index1.indices+s)::I
-(index1::I, s::N where N<:Integer) where I<:Index = Index(index1.indices-s)::I
-(index::I) where I<:Index = Index(-index.indices)::I

# other usefull methods
getindex(index::I where I<:Index, i::N where N<:Integer) = index.indices[i]
getindex(index::I, ::Colon) where I<:Index = index.indices
getindex(A::Array{T,d} where T,i::Index{d}) where d = A[i+1]
setindex!(index::I where I<:Index, value, key) = index.indices[key] = value
maximum(index::I where I<:Index) = maximum(index.indices)
indmax(index::I where I<:Index) = indmax(index.indices)
sum(index::I where I<:Index) = sum(index.indices)
prod(index::I where I<:Index) = prod(max.(1,index.indices))
diff(index1::I, index2::I) where I<:Index = count(!,index1.indices.==index2.indices) # returns number of indices that is different
length(index::I where I<:Index) = length(index.indices)
copy(index::I where I<:Index) = Index(copy(index.indices))
hash(index::I where I<:Index, h::UInt) = hash(index.indices, hash(:Index, h)) # needed for looking things up in a dict

# output formatting
function show(io::IO, index::Index)
    d = ndims(index)
    if d == 1
        print(io, "level $(index.indices)")
    else
        print(io, "index $(index.indices)")
    end
end
