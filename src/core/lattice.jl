## lattice.jl : implementation of rank-1 lattice rules
#
# Basic implementation of rank-1 lattice rules to generate points in [0,1] 
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2019

abstract type AbstractLatticeRule{s} <: AbstractRNG end

## LatticeRule32 ##
mutable struct LatticeRule32{s} <: AbstractLatticeRule{s}
    z::Vector{UInt32}
    n::UInt32
end

"""
    LatticeRule32(z, s, n)
    LatticeRule32(z, s)
    LatticeRule32(z)

Returns a rank-1 lattice rule in `s` dimensions with generating vector `z` and at most `n` points.

Technically, we return an extensible lattice sequence where the `k`-th point is transformed using the radical inverse function. This has the advantage that we can add points to the set without changing the already computed points. When no maximum number of points `n` is provided, we assume `n=2^32`. When no number of dimensions `s` is provided, we assume `s=length(z)`. 

# Examples
```jldoctest
julia> lattice_rule = LatticeRule32([0x00000001,0x00000022], 2, 0x00000037) # Fibonacci lattice rule
LatticeRule32{2}

julia> get_point(lattice_rule, 2)
2-element Array{Float64,1}:
 0.75
 0.5

julia> get_point(lattice_rule, 56) # returns random points beyond n
2-element Array{Float64,1}:
[...]
```

"""
LatticeRule32(z::Vector{UInt32}) = LatticeRule32(z, length(z))
LatticeRule32(z::Vector{UInt32}, s::Integer) = LatticeRule32(z, s, typemax(UInt32))
LatticeRule32(z::Vector{UInt32}, s::Integer, n::Integer) = check_args(z, s, n) && LatticeRule32{s}(view(z, 1:s), convert(UInt32, n))

"""
    LatticeRule32(str, s, n)
    LatticeRule32(str, s)
    LatticeRule32(str)

Returns a rank-1 lattice rule in `s` dimensions with generating vector `z` read from the file `str` and at most `n` points.

Technically, we return an extensible lattice sequence where the `k`-th point is transformed using the radical inverse function. This has the advantage that we can add points to the set without changing the already computed points. When no maximum number of points `n` is provided, we assume `n=2^32`. When no number of dimensions `s` is provided, we assume `s=length(z)`. The file `str` should contain a vector of integer points.

# Examples
```jldoctest
julia> import MultilevelEstimators

julia> z_file = joinpath(dirname(pathof(MultilevelEstimators)), "core", "generating_vectors", "K_3600_32.txt")
[...]

julia> lattice_rule = LatticeRule32(z_file, 16)
LatticeRule32{16}

julia> get_point(lattice_rule, 123)
16-element Array{Float64,1}:
 0.3828125
 0.28125  
 0.15625  
 0.34375  
 0.6015625
 0.1171875
 0.4296875
 0.4609375
 0.2265625
 0.65625  
 0.9609375
 0.6796875
 0.0546875
 0.21875  
 0.1328125
 0.1640625

```

"""
LatticeRule32(s::Int) = check_arg(s) && LatticeRule32(joinpath(@__DIR__(), "generating_vectors", "K_3600_32.txt"), s)
LatticeRule32(str::String, s::Integer) = LatticeRule32(vec(readdlm(str, UInt32)), s)
LatticeRule32(str::String, s::Integer, n::Integer) = LatticeRule32(vec(readdlm(str, UInt32)), s, n)

show(io::IO, lattice_rule::LatticeRule32{s}) where s = print(io, string("LatticeRule32{", s, "}"))

## ShiftedLatticeRule ##
struct ShiftedLatticeRule{L<:AbstractLatticeRule, V<:AbstractVector{<:Real}}
    lattice_rule::L
    shift::V
end

ShiftedLatticeRule(lattice_rule::AbstractLatticeRule{s}) where s = ShiftedLatticeRule(lattice_rule, rand(s))

show(io::IO, s::ShiftedLatticeRule) = print(io, string("ShiftedLatticeRule{", s.lattice_rule, "}"))

## get_point ##
@inline function get_point!(x::Vector{<:AbstractFloat}, lattice_rule::LatticeRule32{s}, k::UInt32) where s
    ϕ = convert(UInt32, reversebits(xor(k, k >> 1)))
	x .= ϕ * 2.0^(-32) * lattice_rule.z .% 1
end

@inline get_point(lattice_rule::LatticeRule32{s}, k::UInt32) where s = k > lattice_rule.n ? rand(Float32, s) : get_point!(Vector{Float32}(undef, s), lattice_rule, k)
@inline get_point(lattice_rule::LatticeRule32{s}, k::Int64) where s = k > lattice_rule.n ? rand(s) : get_point!(Vector{Float64}(undef, s), lattice_rule, convert(UInt32, k))

@inline get_point(shifted_lattice_rule::ShiftedLatticeRule, k::Integer) = get_point(shifted_lattice_rule.lattice_rule, k) .+ shifted_lattice_rule.shift

## reverse_bits ##
reversebits(n::U) where U<:Unsigned = parse(U, reverse(bitstring(n)), base=2)

## input checking ##
function check_args(z::Vector{UInt32}, s::Integer, n::Integer)
	check_larger_than(LatticeRule32, s, "number of dimensions s", 0)
	check_ordered(LatticeRule32, s, length(z)+1, "number of dimensions s", "or equal to the length of the generating vector z")
	check_ordered(LatticeRule32, n, 2^32, "maximum number of points n", "or equal to 2^32, consider implementing a LatticeRule64 type")
end

check_arg(s::Integer) = check_ordered(LatticeRule32, s, 3601, "number of dimensions s", "or equal to 3600, please supply your own generating vector z to proceed")
