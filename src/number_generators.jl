# number_generators.jl : collection of number generators

abstract type NumberGenerator{s} end

abstract type MCgenerator{s} <: NumberGenerator{s} end

abstract type QMCgenerator{s} <: NumberGenerator{s} end

## uniform random number generators

# uniform Monte Carlo generator
struct UniformMCgenerator{s,V<:AbstractVector} <: MCgenerator{s}
    lb::V # lower bound for uniform uncertainties
    ub::V # upper bound for uniform uncertainties
end

# constructors
function UniformMCgenerator(s::N where N<:Integer)
    return UniformMCgenerator{s,Vector{Float64}}(zeros(s),ones(s))
end

function UniformMCgenerator(s::N where N<:Integer,lb::Vector{T},ub::Vector{T}) where T<:AbstractFloat
    @assert length(ub) == s && length(lb) == s
    return UniformMCgenerator{s,Vector{T}}(lb,ub)
end

# utilities
nshifts(::UniformMCgenerator) = 1
reset(::UniformMCgenerator) = Void

# methods
is_valid(U::UniformMCgenerator{s}) where s = length(U.lb) == s && length(U.ub) == s

function show{s}(io::IO,U::UniformMCgenerator{s})
    @assert is_valid(U)
    print(io,"$(ndims(U))-dimensional uniform Monte Carlo sampler with \n")
    print(io," - lower bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.lb)
    print(io," - upper bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.ub)
end

function getPoint(generator::UniformMCgenerator{s},k::N where N<:Integer) where s
    return generator.lb .+ (generator.ub - generator.lb).*rand(s)
end

# uniform Quasi-Monte Carlo generator
struct UniformQMCgenerator{s,q,R<:RandWrapper,V<:AbstractVector} <: QMCgenerator{s}
    generator::R
    lb::V # lower bound for uniform unvertainties
    ub::V # upper bound for uniform uncertainties
end

# constructors
function UniformQMCgenerator(s::N, q::N) where N<:Integer
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return UniformQMCgenerator{s,q,typeof(randlat),Vector{Float64}}(randlat,zeros(s),ones(s))
end

function UniformQMCgenerator(randlat::RandWrapper)
    return UniformQMCgenerator{QMC.ndims(s),nshifts(q),typeof(randlat),Vector{Float64}}(randlat,zeros(s),ones(s))
end

function UniformQMCgenerator(s::N, q::N,lb::Vector{T},ub::Vector{T}) where {N<:Integer,T<:AbstractFloat}
    @assert length(ub) == s && length(lb) == s
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return UniformQMCgenerator{s,q,typeof(randlat),Vector{T}}(randlat,lb,ub)
end

function UniformQMCgenerator(randlat::RandWrapper,lb::Vector{T},ub::Vector{T}) where T<:AbstractFloat
    @assert length(ub) == QMC.ndims(randlat) && length(lb) == QMC.ndims(randlat)
    return UniformQMCgenerator{QMC.ndims(randlat),q = nshifts(randlat),typeof(randlat),T,Vector{T}}(randlat,1.,lb,ub)
end

function UniformQMCgenerator{T<:AbstractFloat}(randlat::RandWrapper,λ::T,lb::Vector{T},ub::Vector{T})
    @assert length(ub) == s && length(lb) == s
    return UniformQMCgenerator{QMC.ndims(randlat),q = nshifts(randlat),typeof(randlat),T,Vector{T}}(randlat,λ,lb,ub)
end

# utilities
nshifts{s,q}(generator::UniformQMCgenerator{s,q}) = q
reset(generator::UniformQMCgenerator) = reset(generator.generator)

# methods
isValid{s}(U::UniformQMCgenerator{s}) = G.λ >= 0.5 && length(lb) == s && length(ub) == s

function show{s,q}(io::IO,U::UniformQMCgenerator{s,q})
    @assert isValid(U)
    str = "$(nshifts(U)) x $(ndims(U))-dimensional uniform Quasi-Monte Carlo sampler with λ= $(U.λ) and \n"
    str *= " - lower bound = $(U.lb)\n"
    str *= " - upper bound = $(U.lb)\n"
    print(io,str)
end

function getPoint{s,q,N<:Integer}(generator::UniformQMCgenerator{s,q},k::N)
    return generator.lb .+ (generator.ub - generator.lb).*getPoint(generator.generator,k)
end

#
# Gaussian random number generators
#

# random number generator
struct GaussianMCgenerator{s,T<:AbstractFloat} <: MCgenerator{s}
    λ::T
end

# constructors
function GaussianMCgenerator{N<:Integer}(s::N)
    return GaussianMCgenerator{s,Float64}(0.5)
end

# utilities
nshifts{s}(generator::GaussianMCgenerator{s}) = 1
reset{s}(generator::GaussianMCgenerator{s}) = Void

# methods
isValid(G::GaussianMCgenerator) = G.λ == 0.5

function show{s}(io::IO,G::GaussianMCgenerator{s})
    @assert isValid(G)
    str = "$(ndims(G))-dimensional Gaussian Monte Carlo sampler with λ= $(G.λ)"
    print(io,str)
end

getPoint{s,N<:Integer}(generator::GaussianMCgenerator{s},k::N) = randn(s)

# randomised lattice rule generator
struct GaussianQMCgenerator{s,q,R<:RandWrapper,T<:AbstractFloat} <: QMCgenerator{s}
    generator::R
    λ::T
end

# constructors
function GaussianQMCgenerator{N<:Integer}(s::N, q::N)
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return GaussianQMCgenerator{s,q,typeof(randlat),Float64}(randlat,1.)
end

function GaussianQMCgenerator(randlat::RandWrapper)
    return GaussianQMCgenerator{QMC.ndims(randlat),nshifts(randlat),typeof(randlat),Float64}(randlat,1.)
end

# utilities
nshifts{s,q}(generator::GaussianQMCgenerator{s,q}) = q
ndims{s}(generator::NumberGenerator{s}) = s

# methods
isValid(G::GaussianQMCgenerator) = G.λ >= 0.5

function show{s,q}(io::IO,G::GaussianQMCgenerator{s,q})
    @assert isValid(G)
    str = "$(nshifts(G)) x $(ndims(G))-dimensional Gaussian Quasi-Monte Carlo sampler with λ= $(G.λ)"
    print(io,str)
end

function getPoint{s,q,N}(generator::GaussianQMCgenerator{s,q},k::N)
    return sqrt(2)*erfinv.(2*getPoint(generator.generator,k)-1)
end

# reset(generator::GaussianQMCgenerator) = reset(generator.generator)














