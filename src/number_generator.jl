# number_generators.jl : collection of number generators

## Distribution ##
abstract type Distribution end

# Uniform
struct Uniform{V<:Vector{T} where {T<:Real}} <: Distribution
    a::V
    b::V
end

Uniform(s) = Uniform(zeros(s),ones(s))

# Normal
struct Normal{V<:Vector{T} where {T<:Real}} <: Distribution 
    μ::V
    σ::V
end

Normal(s) = Normal(zeros(s),ones(s))

# TruncatedNormal
struct TruncatedNormal{V<:Vector{T} where {T<:Real}} <: Distribution
    μ::V
    σ::V
    a::V
    b::V
end

TruncatedNormal(s) = TruncatedNormal(zeros(s),ones(s),-3*ones(s),3*ones(s))

## PointSet ##
abstract type PointSet{s} end

# MC
struct MC{s} <: PointSet{s} end

# QMC
struct QMC{s,q,R} <: PointSet{s}
    generator::R
end

## NumberGenerator ##
struct NumberGenerator{D<:Distribution,P<:PointSet}
    distribution::D
    pointset::P
end

## Type aliases ##
"""
    UniformMCGenerator(s)
    UniformMCGenerator(a,b)

An `s`-dimensional uniform Monte Carlo point set generator with lower bound `a` and upper bound `b`. The default is `a=zeros(s)` and `b=ones(s)`. When the lower and upper bounds are provided, the number of dimensions is extracted from their respective lengths.

# Examples
```jldoctest
julia> UniformMCGenerator(100)
100-dimensional uniform Monte Carlo point generator

julia> UniformMCGenerator(-ones(25),ones(25))
25-dimensional uniform Monte Carlo point generator

```
See also: [UniformQMCGenerator](@ref)
"""
const UniformMCGenerator = NumberGenerator{Uniform,MC}

"""
    UniformQMCGenerator(s,q)
    UniformQMCGenerator(q,a,b)
    UniformQMCGenerator(rwr,a,b)

An `s`-dimensional uniform Quasi-Monte Carlo point set generator with lower bound `a` and upper bound `b` that uses `q` random shifts. The default is `a=zeros(s)` and `b=ones(s)`. When the lower and upper bounds are provided, the number of dimensions is extracted from their respective lengths. Unless otherwise specified with the `rwr` keyword, a default randomized lattice rule is used.

# Examples
```jldoctest
julia> UniformQMCGenerator(100,16)
100-dimensional uniform Quasi-Monte Carlo point generator with 16 random shifts

julia> UniformQMCGenerator(32,-ones(25),ones(25))
25-dimensional uniform Quasi-Monte Carlo point generator with 32 random shifts

julia> rwr = RandWrapper(LatSeq(125),20)
125-dimensional randomized sequence with 20 shifts

julia> UniformQMCGenerator(rwr,zeros(125),ones(125))
125-dimensional uniform Quasi-Monte Carlo point generator with 20 random shifts

```
See also: [UniformMCGenerator](@ref)
"""
const UniformQMCGenerator = NumberGenerator{Uniform,QMC}

"""
    NormalMCGenerator(s)
    NormalMCGenerator(μ,σ)

An `s`-dimensional standard normal Monte Carlo point set generator with mean `μ` and standard devitation `σ`. The default is `μ=zeros(s)` and `σ=ones(s)` When the mean and standard deviation are provided, the number of dimensions is extracted from their respective lengths.

# Examples
```jldoctest
julia> NormalMCGenerator(100)
100-dimensional normal Monte Carlo point generator

julia> NormalMCGenerator(3*ones(25),ones(25))
25-dimensional normal Monte Carlo point generator

```
See also: [NormalQMCGenerator](@ref)
"""
const NormalMCGenerator = NumberGenerator{Normal,MC}

"""
    NormalQMCGenerator(s,q)
    NormalQMCGenerator(q,μ,σ)
    NormalQMCGenerator(rwr,μ,σ)

An `s`-dimensional standard normal Quasi-Monte Carlo point set generator with mean `μ` and standard deviation `σ` that uses `q` random shifts. The default is `μ=zeros(s)` and `σ=ones(s)`. When the mean and standard deviation are provided, the number of dimensions is extracted from their respective lengths. Unless otherwise specified with the `rwr` keyword, a default randomized lattice rule is used. An inverse mapping is used to map the QMC points to ``\\mathbb{R}^s``.

# Examples
```jldoctest
julia> NormalQMCGenerator(100,16)
100-dimensional normal Quasi-Monte Carlo point generator with 16 random shifts

julia> NormalQMCGenerator(32,3*ones(25),ones(25))
25-dimensional normal Quasi-Monte Carlo point generator with 32 random shifts

julia> rwr = RandWrapper(LatSeq(125),20)
125-dimensional randomized sequence with 20 shifts

julia> NormalQMCGenerator(rwr,zeros(125),ones(125))
125-dimensional normal Quasi-Monte Carlo point generator with 20 random shifts

```
See also: [NormalMCGenerator](@ref)
"""
const NormalQMCGenerator = NumberGenerator{Normal,QMC}

"""
    TruncatedNormalMCGenerator(s)
    TruncatedNormalMCGenerator(μ,σ,a,b)

An `s`-dimensional truncated standard normal Monte Carlo point set generator with mean `μ`, standard devitation `σ`, lower bound `a` and upper bound `b`. The default is `μ=zeros(s)`, `σ=ones(s)`, `a=-3*ones(s)` and `b=3*ones(s)`. When the mean, standard deviation, lower and upper bound are provided, the number of dimensions is extracted from their respective lengths.

# Examples
```jldoctest
julia> TruncatedNormalMCGenerator(100)
100-dimensional truncated normal Monte Carlo point generator

julia> TruncatedNormalMCGenerator(3*ones(25),ones(25))
25-dimensional truncated normal Monte Carlo point generator

```
See also: [TruncatedNormalQMCGenerator](@ref)
"""
const TruncatedNormalMCGenerator = NumberGenerator{TruncatedNormal,MC}

"""
    TruncatedNormalQMCGenerator(s,q)
    TruncatedNormalQMCGenerator(q,μ,σ,a,b)
    TruncatedNormalQMCGenerator(rwr,μ,σ,a,b)

    An `s`-dimensional truncated standard normal Quasi-Monte Carlo point set generator with mean `μ`, standard deviation `σ`, lower bound `a` and upper bound `b` that uses `q` random shifts. The default is `μ=zeros(s)`, `σ=ones(s)`, `a=-3*ones(s)` and `b=3*ones(s)`. When the mean, standard deviation, lower and upper bounds are provided, the number of dimensions is extracted from their respective lengths. Unless otherwise specified with the `rwr` keyword, a default randomized lattice rule is used. An inverse mapping is used to map the QMC points to ``\\mathbb{R}^s``.

# Examples
```jldoctest
julia> TruncatedNormalQMCGenerator(100,16)
100-dimensional truncated normal Quasi-Monte Carlo point generator with 16 random shifts

julia> TruncatedNormalQMCGenerator(32,3*ones(25),ones(25),ones(25),5*ones(25))
25-dimensional truncated normal Quasi-Monte Carlo point generator with 32 random shifts

julia> rwr = RandWrapper(LatSeq(125),20)
125-dimensional randomized sequence with 20 shifts

julia> TruncatedNormalQMCGenerator(rwr,zeros(125),ones(125),-2*ones(125),2*ones(125))
125-dimensional truncated normal Quasi-Monte Carlo point generator with 20 random shifts

```
See also: [TruncatedNormalMCGenerator](@ref)
"""
const TruncatedNormalQMCGenerator = NumberGenerator{TruncatedNormal,QMC}

## constructors ##

# input checks
assert_length(::Type{Uniform},params) = length(params) == 2 || throw(ArgumentError("to use Uniform points, supply a lower and an upper bound, got $(length(params)) inputs"))
assert_length(::Type{Normal},params) = length(params) == 2 || throw(ArgumentError("to use Normal points, supply a mean and a standard deviation, got $(length(params)) inputs"))
assert_length(::Type{TruncatedNormal},params) = length(params) == 4 || throw(ArgumentError("to use TruncatedNormal points, supply a mean, standard deviation, lower and upper bound, got $(length(params)) inputs"))

assert_bounds(::Type{Uniform},params) = all(params[1].<params[2]) || throw(ArgumentError("upper bound must be larger than lower bound"))
assert_bounds(::Type{Normal},params) = all(params[2].>=0) || throw(ArgumentError("to use normal random numbers, stadard deviation must be positive")) 
assert_bounds(::Type{TruncatedNormal},params) = ( all(params[3].<params[4]) || throw(ArgumentError("upper bound must be larger than lower bound")) ) && ( all(params[2].>=0) || throw(ArgumentError("to use truncated normal random numbers, standard deviation must be positive")) )

check_ndims(s) = (s > 0) || throw(BoundsError("number of dimenions must be positive, got $(s)"))
check_nshifts(q) = (q > 0 ) || throw(BoundsError("number of shifts must be positive, got $(q)"))

assert_lengths(params,s) = all(length.(params).==s) || throw(DimensionMismatch("inputs must have same length"))

# MC generator, specify number of dimensions
function NumberGenerator{D,P}(s::N where {N<:Integer}) where {D<:Distribution,P<:MC}
    check_ndims(s)
    NumberGenerator(D(s),MC{s}())
end

# MC generator, specify lower and upper bounds
function NumberGenerator{D,P}(params::Vector...) where {D<:Distribution,P<:MC}
    assert_length(D,params)
    s = length(params[1])
    check_ndims(s)
    assert_lengths(params,s)
    assert_bounds(D,params)
    NumberGenerator(D{promote_type(typeof.(params)...)}(promote(params...)...),MC{s}())
end


# QMC generator, specify number of dimensions and number of shifts
function NumberGenerator{D,P}(s::N where {N<:Integer},q::N where {N<:Integer}) where {D<:Distribution,P<:QMC}
    check_ndims(s)
    check_nshifts(q)
    lat = LatSeq(s)
    rwr = RandWrapper(lat,q)
    NumberGenerator{D,P}(rwr)
end

# QMC generator, specify point generator
NumberGenerator{D,P}(rwr::R) where {R<:RandWrapper,D<:Distribution,P<:QMC} = NumberGenerator(D(ndims(rwr)),QMC{ndims(rwr),nshifts(rwr),R}(rwr))

# QMC generator, specify lower and upper bounds and number of shifts
function NumberGenerator{D,P}(q::N where {N<:Integer},params::Vector{T}...) where {T<:Real,D<:Distribution,P<:QMC}
    s = length(params[1])
    check_ndims(s)
    check_nshifts(q)
    lat = LatSeq(s)
    rwr = RandWrapper(lat,q)
    NumberGenerator{D,P}(rwr,params...)
end

# QMC generator, specify lower and upper bounds and point generator
function NumberGenerator{D,P}(rwr::R,params::Vector{T}...) where {R<:RandWrapper,T<:Real,D<:Distribution,P<:QMC}
    assert_length(D,params)
    s = ndims(rwr)
    assert_lengths(params,s)
    assert_bounds(D,params)
    NumberGenerator(D{promote_type(typeof.(params)...)}(promote(params...)...),QMC{s,nshifts(rwr),R}(rwr))
end

## get_point
function get_point(nb::NumberGenerator{D,P} where {D,P},k::N where {N<:Integer})
    point =  get_point(nb.pointset,k)
    transform(nb.distribution,point)
end

get_point(::MC{s},k::N where {N<:Integer}) where {s} = rand(s)
get_point(pointset::QMC,k::N where {N<:Integer}) = getPoint(pointset.generator,k)

transform(u::Uniform,point) = u.a .+ diagm(u.b - u.a) * point
transform(n::Normal,point) = n.μ .+ diagm(n.σ) * Φ⁻¹.(point)
function transform(t::TruncatedNormal,point)
    α = (t.a-t.μ)./t.σ
    β = (t.b-t.μ)./t.σ
    t.μ .+ diagm(t.σ) * Φ⁻¹.(Φ.(α) .+ diagm(Φ.(β)-Φ.(α)) * point )
end

Φ⁻¹(x::T where {T<:Real}) = √2*erfinv(2*x-1)
Φ(x::T where {T<:Real}) = 1/2*(1+erf(x/√2))

## new
new(nb::NumberGenerator{D,P} where {D,P<:MC}) = nb
function new(nb::NumberGenerator{D,P} where {D,P<:QMC})
    nb2 = deepcopy(nb)
    nb2.pointset.generator = RandWrapper(nb.pointset.generator.generator,nshifts(nb.pointset.generator))
end

## show methods
show(io::IO,nb::NumberGenerator{D,MC{s}}) where {D,s} = print(io,"$(s)-dimensional $(nb.distribution) Monte Carlo point generator")
show(io::IO,nb::NumberGenerator{D,QMC{s,q,R}}) where {D,s,q,R} = print(io,"$(s)-dimensional $(nb.distribution) Quasi-Monte Carlo point generator with $(q) random shifts")
show(io::IO,::Uniform) = print(io,"uniform")
show(io::IO,::Normal) = print(io,"normal")
show(io::IO,::TruncatedNormal) = print(io,"truncated normal")
