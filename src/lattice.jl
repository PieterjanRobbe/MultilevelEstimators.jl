## lattice.jl : implementation of rank-1 lattice rules
#
# Basic implementation of rank-1 lattice rules to generate points in [0,1] 
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

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

When no maximum number of points `n` is provided, we assume `n=2^32`. When no number of dimensions `s` is provided, we assume `s=length(z)`. 

!!! info

    Technically, we return an extensible lattice sequence where the `k`-th point is transformed using the radical inverse function. This has the advantage that we can add points to the set without changing the already computed points.

More generating vectors can be found online [here](https://web.maths.unsw.edu.au/~fkuo/lattice/index.html) or [here](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/).

# Examples
```jldoctest; setup = :(using MultilevelEstimators; import Random; Random.seed!(1))
julia> lattice_rule = LatticeRule32([0x00000001,0x00000022], 2, 0x00000037) # Fibonacci lattice rule
LatticeRule32{2}

julia> get_point(lattice_rule, 2)
2-element Array{Float64,1}:
 0.75
 0.5

julia> get_point(lattice_rule, 56) # returns random points beyond n
2-element Array{Float64,1}:
 0.23603334566204692
 0.34651701419196046
```
See also: [`get_point`](@ref), [`ShiftedLatticeRule`](@ref)
"""
LatticeRule32(z::Vector{UInt32}) = LatticeRule32(z, length(z))
LatticeRule32(z::Vector{UInt32}, s::Integer) = LatticeRule32(z, s, typemax(UInt32))
LatticeRule32(z::Vector{UInt32}, s::Integer, n::Integer) = check_args(z, s, n) && LatticeRule32{s}(view(z, 1:s), convert(UInt32, n))

"""
    LatticeRule32(str, s, n)
    LatticeRule32(str, s)
    LatticeRule32(str)

Returns a rank-1 lattice rule in `s` dimensions with generating vector `z` read from the file `str` and at most `n` points.

# Examples
```jldoctest; setup = :(using MultilevelEstimators)
julia> z_file = joinpath(dirname(pathof(MultilevelEstimators)), "generating_vectors", "K_3600_32.txt");

julia> lattice_rule = LatticeRule32(z_file, 16)
LatticeRule32{16}

julia> get_point(lattice_rule, 123)
16-element Array{Float64,1}:
 0.3828125
 0.2890625
 0.1484375
 0.3515625
 0.6015625
 0.1171875
 0.4296875
 0.4609375
 0.2265625
 0.6484375
 0.9609375
 0.6796875
 0.0546875
 0.2265625
 0.1328125
 0.1640625

```
See also: [`get_point`](@ref), [`ShiftedLatticeRule`](@ref)
"""
LatticeRule32(s::Int) = check_arg(s) && LatticeRule32(joinpath(@__DIR__(), "generating_vectors", "K_3600_32.txt"), s)
LatticeRule32(str::String, s::Integer) = LatticeRule32(vec(readdlm(str, UInt32)), s)
LatticeRule32(str::String, s::Integer, n::Integer) = LatticeRule32(vec(readdlm(str, UInt32)), s, n)

show(io::IO, lattice_rule::LatticeRule32{s}) where s = print(io, string("LatticeRule32{", s, "}"))

## ShiftedLatticeRule ##
"""
    ShiftedLatticeRule(lat::AbstractLatticeRule)

Returns a shifted rank-1 lattice rule based on the lattice rule `lat`.

# Examples
```jldoctest; setup = :(using MultilevelEstimators; import Random; Random.seed!(1))
julia> z_file = joinpath(dirname(pathof(MultilevelEstimators)), "generating_vectors", "K_3600_32.txt");

julia> lattice_rule = LatticeRule32(z_file, 16)
LatticeRule32{16}

julia> shifted_lattice_rule = ShiftedLatticeRule(lattice_rule)
ShiftedLatticeRule{LatticeRule32{16}}

julia> get_point(shifted_lattice_rule, 0)
16-element Array{Float64,1}:
 0.23603334566204692
 0.34651701419196046
 0.3127069683360675
 0.00790928339056074
 0.4886128300795012
 0.21096820215853596
 0.951916339835734
 0.9999046588986136
 0.25166218303197185
 0.9866663668987996
 0.5557510873245723
 0.43710797460962514
 0.42471785049513144
 0.773223048457377
 0.2811902322857298
 0.20947237319807077

```
See also: [`LatticeRule32`](@ref), [`get_point`](@ref)
"""
struct ShiftedLatticeRule{L<:AbstractLatticeRule, V<:AbstractVector{<:Real}}
    lattice_rule::L
    shift::V
end

ShiftedLatticeRule(lattice_rule::AbstractLatticeRule{s}) where s = ShiftedLatticeRule(lattice_rule, rand(s))

show(io::IO, s::ShiftedLatticeRule) = print(io, string("ShiftedLatticeRule{", s.lattice_rule, "}"))

## get_point ##
@inline function get_point!(x::Vector{<:AbstractFloat}, lattice_rule::LatticeRule32, k::UInt32)
    ϕ = convert(UInt32, reversebits(xor(k, k >> 1)))
	copyto!(x, ϕ * 2.0^(-32) * lattice_rule.z .% 1)
end

"""
    get_point(lat::AbstractLatticeRule, k::Integer)

Get the `k`-th point of the lattice rule `lat`.

```jldoctest; setup = :(using MultilevelEstimators)
julia> lattice_rule = LatticeRule32([0x00000001,0x00000022], 2, 0x00000037) # Fibonacci lattice rule
LatticeRule32{2}

julia> get_point(lattice_rule, 2)
2-element Array{Float64,1}:
 0.75
 0.5

```
See also: [`LatticeRule32`](@ref), [`ShiftedLatticeRule`](@ref)
"""
@inline get_point(lattice_rule::LatticeRule32{s}, k::UInt32) where s = k > lattice_rule.n ? rand(Float32, s) : get_point!(Vector{Float32}(undef, s), lattice_rule, k)

@inline get_point(lattice_rule::LatticeRule32{s}, k::Int) where s = k > lattice_rule.n ? rand(s) : get_point!(Vector{Float64}(undef, s), lattice_rule, convert(UInt32, k))

@inline function get_point!(x::Vector{<:AbstractFloat}, shifted_lattice_rule::ShiftedLatticeRule, k::UInt32)
    ϕ = convert(UInt32, reversebits(xor(k, k >> 1)))
	copyto!(x, (ϕ * 2.0^(-32) * shifted_lattice_rule.lattice_rule.z .+ shifted_lattice_rule.shift) .% 1)
end

@inline get_point(shifted_lattice_rule::ShiftedLatticeRule{<:LatticeRule32{s}}, k::UInt32) where s = k > shifted_lattice_rule.lattice_rule.n ? rand(Float32, s) : get_point!(Vector{Float32}(undef, s), shifted_lattice_rule, k)

@inline get_point(shifted_lattice_rule::ShiftedLatticeRule{<:LatticeRule32{s}}, k::Int) where s = k > shifted_lattice_rule.lattice_rule.n ? rand(s) : get_point!(Vector{Float64}(undef, s), shifted_lattice_rule, convert(UInt32, k))

## reverse_bits ##
reversebits(n::U) where U<:Unsigned = parse(U, reverse(bitstring(n)), base=2)

## input checking ##
function check_args(z::Vector{UInt32}, s::Integer, n::Integer)
	check_larger_than(LatticeRule32, s, "number of dimensions s", 0)
	check_ordered(LatticeRule32, s, length(z) + 1, "number of dimensions s", "or equal to the length of the generating vector z")
	check_ordered(LatticeRule32, n, typemax(UInt32) + 1, "maximum number of points n", "or equal to 2^32, consider implementing a LatticeRule64 type")
end

check_arg(s::Integer) = check_ordered(LatticeRule32, s, 3601, "number of dimensions s", "or equal to 3600, please supply your own generating vector z to proceed")
