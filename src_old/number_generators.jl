# number_generators.jl : collection of number generators

"""
NumberGenerator{s}

Supertype for a number generator in `s` dimensions.
"""
abstract type NumberGenerator{s} end

"""
MCgenerator{s}

Supertype for a Monte Carlo number generator in `s` dimensions. Subtype of `NumberGenerator{s}`.
"""
abstract type MCgenerator{s} <: NumberGenerator{s} end

"""
QMCgenerator{s}

Supertype for a Quasi-Monte Carlo number generator in `s` dimensions. Subtype of `NumberGenerator{s}`.
"""
abstract type QMCgenerator{s} <: NumberGenerator{s} end

## uniform random number generators

# uniform Monte Carlo generator
"""
UniformMCgenerator{s,V}

Generator for uniformly distributed Monte Carlo points in an `s`-dimensional hypercube. 
"""
struct UniformMCgenerator{s,V<:AbstractVector} <: MCgenerator{s}
    lb::V # lower bound for uniform uncertainties
    ub::V # upper bound for uniform uncertainties
end

# constructors
"""
UniformMCgenerator(s)

Construct a generator for uniformly distributed Monte Carlo points in an `s`-dimensional hypercube. 

# Examples

```jldoctest
julia> UniformMCgenerator(1400)
1400-dimensional uniform Monte Carlo sampler with 
- lower bound = 
1400-element Array{Float64,1}:
0.0
0.0
⋮  
0.0
0.0
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformMCgenerator(s::N where N<:Integer)
    s <= 0 && throw(ArgumentError("number of dimensions s must be positive, got $(s)"))
    return UniformMCgenerator{s,Vector{Float64}}(zeros(s),ones(s))
end

"""
UniformMCgenerator(s,lb,ub)

Construct a generator for uniformly distributed Monte Carlo points in an `s`-dimensional hypercube, scaled by a lower bound `lb` and upper bound `ub`. 

# Examples

```jldoctest
julia> s = 1400
1400

julia> UniformMCgenerator(s,-ones(s),ones(s))
1400-dimensional uniform Monte Carlo sampler with 
- lower bound = 
1400-element Array{Float64,1}:
-1.0
-1.0
⋮  
-1.0
-1.0 
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformMCgenerator(s::N where N<:Integer,lb::Vector{T},ub::Vector{T}) where T<:AbstractFloat
    length(lb) != s && throw(ArgumentError("lower bound lb must be of length s = $(s), got $(length(lb))"))
    length(ub) != s && throw(ArgumentError("upper bound ub must be of length s = $(s), got $(length(ub))"))
    any(lb.>ub) && throw(ArgumentError("upper bound lb must be larger than lower bound lb"))
    return UniformMCgenerator{s,Vector{T}}(lb,ub)
end

# utilities
nshifts(::UniformMCgenerator) = 1

# methods
is_valid(U::UniformMCgenerator{s}) where s = length(U.lb) == s && length(U.ub) == s

function show(io::IO,U::UniformMCgenerator{s} where s)
    print(io,"$(ndims(U))-dimensional uniform Monte Carlo sampler with \n")
    print(io," - lower bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.lb)
    print(io,"\n - upper bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.ub)
end

function getPoint(generator::UniformMCgenerator{s},k::N where N<:Integer) where s
    return generator.lb .+ (generator.ub - generator.lb).*rand(s)
end

# uniform Quasi-Monte Carlo generator
"""
UniformQMCgenerator{s,V}

Generator for Quasi-Monte Carlo points in an `s`-dimensional hypercube. 
"""
struct UniformQMCgenerator{s,q,R<:RandWrapper,V<:AbstractVector} <: QMCgenerator{s}
    generator::R
    lb::V # lower bound for uniform unvertainties
    ub::V # upper bound for uniform uncertainties
end

# constructors
"""
UniformQMCgenerator(s, q)

Construct a generator for Quasi-Monte Carlo points with `q` shifts in an `s`-dimensional hypercube. 

# Examples

```jldoctest
julia> UniformQMCgenerator(1400,16)
1400-dimensional Quasi-Monte Carlo sampler based on a 1400-dimensional randomized sequence with 16 shifts that generates points in a hypercube with 
- lower bound = 
1400-element Array{Float64,1}:
0.0
0.0
⋮  
0.0
0.0
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformQMCgenerator(s::N, q::N) where N<:Integer
    s <= 0 && throw(ArgumentError("number of dimensions s must be positive, got $(s)"))
    q <= 0 && throw(ArgumentError("number of shifts q must be positive, got $(q)"))
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return UniformQMCgenerator{s,q,typeof(randlat),Vector{Float64}}(randlat,zeros(s),ones(s))
end

"""
UniformQMCgenerator(randlat)

Construct a generator for Quasi-Monte Carlo points using the randomly shifted lattice rule `randlat`. 

# Examples

```jldoctest
julia> lat = LatSeq(1400)
1400-dimensional lattice sequence

julia> randlat = RandWrapper(lat,16)
1400-dimensional randomized sequence with 16 shifts

julia> UniformQMCgenerator(randlat)
1400-dimensional Quasi-Monte Carlo sampler based on a 1400-dimensional randomized sequence with 16 shifts that generates points in a hypercube with 
- lower bound = 
1400-element Array{Float64,1}:
0.0
0.0
⋮  
0.0
0.0
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformQMCgenerator(randlat::RandWrapper)
    s = QMC.ndims(randlat)
    q = nshifts(randlat)
    return UniformQMCgenerator{s,q,typeof(randlat),Vector{Float64}}(randlat,zeros(s),ones(s))
end

"""
UniformQMCgenerator(s, q, lb, ub)

Construct a generator for Quasi-Monte Carlo points with `q` shifts in an `s`-dimensional hypercube with lower bound `lb` and upper bound `ub`. 

# Examples

```jldoctest
julia> s = 1400
1400

julia> UniformQMCgenerator(s,16,-ones(s),ones(s))
1400-dimensional Quasi-Monte Carlo sampler based on a 1400-dimensional randomized sequence with 16 shifts that generates points in a hypercube with 
- lower bound = 
1400-element Array{Float64,1}:
-1.0
-1.0
⋮  
-1.0
-1.0
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformQMCgenerator(s::N, q::N,lb::Vector{T},ub::Vector{T}) where {N<:Integer,T<:AbstractFloat}
    s <= 0 && throw(ArgumentError("number of dimensions s must be positive, got $(s)"))
    q <= 0 && throw(ArgumentError("number of shifts q must be positive, got $(q)"))
    length(lb) != s && throw(ArgumentError("lower bound lb must be of length s = $(s), got $(length(lb))"))
    length(ub) != s && throw(ArgumentError("upper bound ub must be of length s = $(s), got $(length(ub))"))
    any(lb.>ub) && throw(ArgumentError("upper bound lb must be larger than lower bound lb"))
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return UniformQMCgenerator{s,q,typeof(randlat),Vector{T}}(randlat,lb,ub)
end

"""
UniformQMCgenerator(randlat, lb, ub)

Construct a generator for Quasi-Monte Carlo points using the randomly shifted lattice rule `randlat` and scale the points between `lb` and `ub`. 

# Examples

```jldoctest
julia> s = 1400
1400

julia> lat = LatSeq(s)
1400-dimensional lattice sequence

julia> randlat = RandWrapper(lat,16)
1400-dimensional randomized sequence with 16 shifts

julia> UniformQMCgenerator(randlat, -ones(s), ones(s))
1400-dimensional Quasi-Monte Carlo sampler based on a 1400-dimensional randomized sequence with 16 shifts that generates points in a hypercube with 
- lower bound = 
1400-element Array{Float64,1}:
-1.0
-1.0
⋮  
-1.0
-1.0
- upper bound = 
1400-element Array{Float64,1}:
1.0
1.0
⋮  
1.0
1.0

```
"""
function UniformQMCgenerator(randlat::RandWrapper,lb::Vector{T},ub::Vector{T}) where T<:AbstractFloat
    s = QMC.ndims(randlat)
    length(lb) != s && throw(ArgumentError("lower bound lb must be of length s = $(s), got $(length(lb))"))
    length(ub) != s && throw(ArgumentError("upper bound ub must be of length s = $(s), got $(length(ub))"))
    any(lb.>ub) && throw(ArgumentError("upper bound lb must be larger than lower bound lb"))
    return UniformQMCgenerator{s,nshifts(randlat),typeof(randlat),Vector{T}}(randlat,lb,ub)
end

# utilities
nshifts(generator::UniformQMCgenerator{s,q} where s) where q = q

# methods
is_valid(U::UniformQMCgenerator{s}) where s = length(U.lb) == s && length(U.ub) == s

function show(io::IO,U::UniformQMCgenerator)
    print(io,"$(ndims(U))-dimensional Quasi-Monte Carlo sampler based on a $(U.generator) that generates points in a hypercube with \n")
    print(io," - lower bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.lb)
    print(io,"\n - upper bound = \n")
    show(IOContext(io,limit=true),"text/plain",U.ub)
end

function getPoint{s,q,N<:Integer}(generator::UniformQMCgenerator{s,q},k::N)
    return generator.lb .+ (generator.ub - generator.lb).*getPoint(generator.generator,k)
end

## Gaussian random number generators

# Gaussian Monte Carlo generator
"""
GaussianMCgenerator{s,V}

Generator for normal distributed Monte Carlo points in \mathbb{R}^`s`. 
"""
struct GaussianMCgenerator{s} <: MCgenerator{s}
end

# constructors
"""
GaussianMCgenerator(s)

Construct a generator for normal distributed Monte Carlo points in \mathbb{R}^`s`. 

# Examples

```jldoctest
julia> GaussianMCgenerator(1400)
1400-dimensional Gaussian Monte Carlo sampler  

```
"""
function GaussianMCgenerator(s::N where N<:Integer)
    s <= 0 && throw(ArgumentError("number of dimensions s must be positive, got $(s)"))
    return GaussianMCgenerator{s}()
end

# utilities
nshifts(generator::GaussianMCgenerator{s} where s) = 1

# methods
is_valid(::GaussianMCgenerator) = true

function show(io::IO,G::GaussianMCgenerator{s} where s)
    str = "$(ndims(G))-dimensional Gaussian Monte Carlo sampler"
    print(io,str)
end

getPoint{s,N<:Integer}(generator::GaussianMCgenerator{s},k::N) = randn(s)

# Gaussian Quasi-Monte Carlo generator
"""
GaussianQMCgenerator{s,q,R}

Generator for normal distributed Monte Carlo points with `q` shifts in \mathbb{R}^`s`, using the QMC points generator of type `R`. The points are mapped from the unit cube using the inverse normal CDF. 
"""
struct GaussianQMCgenerator{s,q,R<:RandWrapper} <: QMCgenerator{s}
    generator::R
end

# constructors
"""
GaussianQMCgenerator(s, q)

Construct a generator for normal distributed Quasi-Monte Carlo points with `q` shifts in \mathbb{R}^`s`.

# Examples

```jldoctest
julia> UniformQMCgenerator(1400,16)
1400-dimensional Gaussian Quasi-Monte Carlo sampler with 16 shifts

```
"""
function GaussianQMCgenerator{N<:Integer}(s::N, q::N)
    s <= 0 && throw(ArgumentError("number of dimensions s must be positive, got $(s)"))
    q <= 0 && throw(ArgumentError("number of shifts q must be positive, got $(q)"))
    lat = LatSeq(s)
    randlat = RandWrapper(lat,q)
    return GaussianQMCgenerator{s,q,typeof(randlat)}(randlat)
end

"""
GaussianQMCgenerator(randlat)

Construct a generator for normal distributed Quasi-Monte Carlo points in \mathbb{R}^`s` that draws points from the QMC point set generator `randlat`.

# Examples

```jldoctest
julia> lat = LatSeq(16)
16-dimensional lattice sequence

julia> randlat = RandWrapper(lat,8)
16-dimensional randomized sequence with 8 shifts

julia> UniformQMCgenerator(randlat)
1400-dimensional Gaussian Quasi-Monte Carlo sampler with 16 shifts

```
"""
function GaussianQMCgenerator(randlat::RandWrapper)
    return GaussianQMCgenerator{QMC.ndims(randlat),nshifts(randlat),typeof(randlat)}(randlat)
end

# utilities
nshifts(generator::GaussianQMCgenerator{s,q} where s) where q = q
ndims(generator::NumberGenerator{s}) where s = s

# methods
is_valid(::GaussianQMCgenerator) = true

function show(io::IO,G::GaussianQMCgenerator{s,q} where {s,q})
    str = "$(QMC.ndims(G))-dimensional Gaussian Quasi-Monte Carlo sampler based on a $(G.generator)"
    print(io,str)
end

function getPoint{s,q,N}(generator::GaussianQMCgenerator{s,q},k::N)
    return sqrt(2)*erfinv.(2*getPoint(generator.generator,k)-1)
end
