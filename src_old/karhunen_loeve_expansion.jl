## karhunen_loeve_expansion.jl : Karhunen-Loeve expansion of Gaussian random fields

## Karhunen-Loeve (KL) expansion ##
"""
KarhunenLoeveExpansion

Representation of a Karhunen-Loeve expansion.
"""
struct KarhunenLoeveExpansion{C<:CovarianceFunction,N,V<:AbstractVector,F<:AbstractVector} <: GaussianRandomFieldSampler
    cov::C # covariance function
    nterms::N # number of terms
    eigenval::V # KL eigenvalues
    eigenfunc::F # KL eigenfunctions
end

# utilities
ndims(K::KarhunenLoeveExpansion{C{d}}) = d
nterms(K::KarhunenLoeveExpansion) = K.nterms

# methods
function show(io::IO, K::KarhunenLoeveExpansion)
    print(io, "Karhunen-Lo\'eve expansion with $(K.nterms) terms for a $(K.cov)")
end

## Main methods for creating and composing the KL expansion ##

# For separable kernels, we assume that the 1d points are given.
# When composing the random field, we take a tensor product of the 1d eigenfunctions.
# By specifying which eigenfunctions must be tensorized, we can do MIMC.
# In the special case of the SeparableExponentialCovarianceFunction with p=1,
# we have analytic expressions for eigenvalues and eigenfunctions.
# TODO how to do non-separable covariance in MIMC case? => interpolate
# TODO or have x represented by a matrix instead of a vector...
# also TODO is all of this necessary? We can just generate the random field at the finest grid
# and do restriction... => not if #terms is mem intensive...

"""
KarhunenLoeveExpansion(cov, x, nterms)

Compute the Karhunen-Lo\'eve expansion for the covariance function `cov` in all points `x`, where `nterms` terms are used on each level. The points `x` must be explicitly given at each level. For separable covariance functions, this must be the 1d points. For non-separable covariance functions, this must be an `n`-by-`d` matrix. The number of terms can be an integer, or a Vector with integers that specifies how many terms should be used on each level. The number of terms can be an integer, or a Vector with integers that specifies how many terms should be used on each level.  

"""
# TODO add example
function KarhunenLoeveExpansion(cov::C where C<:CovarianceFunction, x::Vector{AbstractArray{T}} where T, nterms::Vector{N} where N)
    in(:apply,fieldnames(cov)) || throw(ArgumentError("no function apply in covariance function $(cov)"))
    # TODO more error handling here...

    if cov.p == 1 && ( typeof(cov) == SeparableMaternCovarianceFunction || ( typeof(cov) == MaternCovarianceFunction && ndims(cov) == 1 ) )
        return compute_kl_analytical(cov,x,nterms)
    elseif typeof(C) <: SeparableCovarianceFunction
        return compute_separable_kl(cov,x,nterms)
    else
        return compute_kl(cov,x,nterms)
    end
end

# alias for integer nterms
KarhunenLoeveExpansion(cov::M where M<:CovarianceFunction, x::Vector{Vector{T}} where T, nterms::N where N<:Number) = KarhunenLoeveExpansion(cov,x,nterms*ones(length(x)))


""" 
compose(kl, xi, index)

Compose the Karhunen-Lo\'eve expansion `kl` at index `index` based on the sample `xi`.

"""
# TODO add examples
# TODO switch between separable/non-separable
function compose(kl::KarhunenLoeveExpansion{d} where d, xi::Vector{T}, index::Index)
    nterms = length(xi)

    # first term for type-stability of k
    size(kl.eigenfunc)
    size(kl.eigenfunc[index_[1]+1])
    v = kl.eigenfunc[index_[1]+1][kl.eigenval[1][1][1],:]::Array{T,1}
    for p = 2:d
        @inbounds v = kron(kl.eigenfunc[index_[p]+1][kl.eigenval[1][1][p],:],v)::Array{T,1}
    end
    k = xi[1]*sqrt(kl.eigenval[1][2])*v
    # rest of the terms
    for i in 2:nterms
        @inbounds v = kl.eigenfunc[index_[1]+1][kl.eigenval[i][1][1],:]::Array{T,1}
        for p = 2:d
            @inbounds v = kron(kl.eigenfunc[index_[p]+1][kl.eigenval[i][1][p],:],v)::Array{T,1}
        end
        @inbounds k += xi[i]*sqrt(kl.eigenval[i][2])*v
    end
    return k
end

## internals ##
function compute_kl_analytical(cov::C where C<:CovarianceFunction, x::Vector{AbstractArray{T}} where T,nterms::Vector{N} where N)
    ω = findroots(cov.λ, maximum(nterms)+1)
    θ = 2*cov.σ^2*cov.λ./(cov.λ^2*ω.^2+1)
    eigenval = find_d_dim_eigenvalues(d,θ,nterms)

    # calculate 1d eigenfunctions
    eigenfunc = Array{T,2}[]
    ω = ω[1:maximum(nterms)]
    n = sqrt(2)/2*sqrt.(1./ω.*(cov.λ^2*ω.^2.*cos.(ω).*sin.(ω)+cov.λ^2*ω.^3-2*cov.λ*ω.*cos.(ω).^2-cos.(ω).*sin.(ω)+ω)+2*cov.λ)
    for i in 1:length(x)
        push!(eigenfunc,diagm(1./n)*( sin.(ω*collect(x[i])') + kernel.λ*diagm(ω)*cos.(ω*collect(x[i])') ))
    end

    return KarhunenLoeveExpansion{typeof(cov),Vector{N},typeof(eigenval),typeof(eigenfunc)}(cov,nterms,eigenval,eigenfunc)
end

function compute_separable_kl(cov::C where C<:SeparableCovarianceFunction, x::Vector{AbstractArray{T}} where T,nterms::Vector{N} where N)
    d = ndims(cov)

    pts = collect(x[end])
    θ, f = nystrom(kernel, pts, maximum(nterms))	
    eigenval = find_d_dim_eigenvalues(d,θ,nterms)

    # calculate 1d eigenfunctions
    eigenfunc = Array{T,2}[]
    for i in 1:length(x)
        ipts = x[i]
        m = length(ipts)
        ar = zeros(T,mkl,m)
        for r = 1:nterms[i]
            ip = interpolate((pts,), normalize_ef(f[r,:],pts), Gridded(Linear()))
            ar[r,:] = ip[ipts]
        end
        push!(eigenfunc, ar)
    end

    return KLExpansion{typeof(cov),Vector{N},typeof(eigenval),typeof(eigenfunc)}(cov,nterms,eigenval,eigenfunc)
end

function compute_kl(cov::C where C<:SeparableCovarianceFunction, x::Vector{AbstractArray{T}} where T,nterms::Vector{N} where N)
end

# find the d-dimensional eigenvalues by adaptive search
function find_d_dim_eigenvalues(d::N, θ::Vector{T}, nterms::N) where {T<:AbstractFloat,N<:Integer}
    old = zeros(nterms,d+1)
    active = Dict{Index{d,Vector{N}},T}()
    active[one(Index{d,Vector{N}})] = θ[1]^d
    ptr = 0
    while ptr < nterms
        ptr += 1
        index = collect(keys(active))[indmax(values(active))] # find maximum eigenvalue in active set
        old[ptr,1:d] = index
        old[ptr,d+1] = active[index]
        for p = 1:d # check for new admissable indices
            index[p] = index[p] + 1
            if ~haskey(active,index) # add to admissable set, if not already added
                active[copy(index)] = prod(θ[index])
            end
            index[p] = index[p] - 1
        end
        delete!(active,index)
    end
    return old
end


function normalize_ef{T}(f::Vector{T},pts)
    A = sqrt(trapezium_rule(pts,f.^2))
    return vec(f/A)
end

function trapezium_rule(pts, f)
    return 1/2*diff(pts)'*(f[2:end]+f[1:end-1])
end

# helper function to compute nodes and weights of Gauss quadrature on [0,1]
function gauss_legendre_0_1(m)
    nodes, weights = gausslegendre(m)
    nodes = 0.5 + nodes/2 # scale roots
    weights = 0.5*weights # scale weights

    return nodes, weights
end

# nodes and weights on [0,1]^n
function gauss_legendre_0_1_n(m)
    nodes, weights = gauss_legendre_0_1(m)    

    Nd = collect(Base.product([1:m for i = 1:n]...))[:]
    N = zeros(length(Nd),n)
    for i = 1:length(Nd)
        N[i,:] = flipdim(nodes[Nd[i]...],1)
    end
    W = kron([weights for i = 1:n]...)

    return N,W
end

# scale function for nodes
function scale_to_0_1(nodes::Matrix{T} where T)
    mx = maximum(nodes,1)
    mn = minimum(nodes,1)

    for n in 1:size(nodes,2)
        nodes[:,n] = (nodes[:,n] - mn[n])/(mx[n]-mn[n])
    end 

    return nodes
end

# nystrom method
function nystrom(cov::CovarianceFunction,x::AbstractArray{T},nterms::N) where {T,N}
    # compute covariance matrix
    nodes, weights = gauss_legendre_0_1_n(nterms) # TODO how many terms ???
    K = cov.apply(nodes,nodes) # TODO check implementation of apply function

    # obtain eigenvalues and eigenfunctions
    D = diagm(weights)
    Dsqrt = sqrt.(D)
    Z = Symmetric(triu(Dsqrt*K*Dsqrt))
    eigenval, eigenfunc = eig(Z)
    eigenval = eigenval*kernel.σ^2
    #eigenval = real(eigenval)
    idx = sortperm(eigenval,rev=true)
    sort!(eigenval,rev=true)
    eigenfunc = eigenfunc[:,idx]

    f = weights.*(diagm(sqrt.(1./weights))*eigenfunc)
    lambda = 1./eigenval

    # nystrom method
    m = length(x)
    eigenfunc = zeros(nterms,m)
    #x = 1/2/m:1/m:1-1/2/m
    K = applyKernel(kernel,x,nodes')*kernel.σ^2
    for j in 1:nterms
        @inbounds eigenfunc[j,:] = lambda[j]*K*f[:,j]
    end

    return eigenval, eigenfunc
end

# find all positive (>0) zeros of the transcendental function tan(ω) = 2*λ*ω/(λ^2*ω^2-1)
function findroots{T<:AbstractFloat,N<:Integer}(λ::T, n::N)
    λ > 0 || error("λ must be positive (λ ≥ 0)!")
    n > 0 || error("requested number of roots must be larger than zero!")

    # define the transcendental function
    f(ω) = (λ^2*ω^2-1)*sin(ω)-2*λ*ω*cos(ω)

    # find range around singularity 1/λ
    left_point_of_range = (2*floor(1/(π*λ)-1/2))*π/2 # left odd multiple of π/2
    right_point_of_range = (2*ceil(1/(π*λ)-1/2)+1)*π/2 # right odd multiple of π/2

    # find roots before 1/λ, if any
    if left_point_of_range ≠ π/2
        roots = zeros(min(n,round(UInt64,floor(abs(1/λ/π-1/2)))))
        left_point = π/2
        right_point = 1*π
        for i = 1:length(roots)
            @inbounds roots[i] = bisect_root(f,left_point,right_point)[1]
            right_point = roots[i] + π
            left_point = left_point + π
        end
    else
        roots = zeros(0)
    end

    # find roots inside range around 1/λ
    ( length(roots) ≥ n || floor(1/(π*λ)-1/2) < 0 ) || 
    push!(roots,bisect_root(f,left_point_of_range+eps(T),1/λ)[1]) # first intersection point
    ( length(roots) ≥ n || ceil(1/(π*λ)-1/2) < 0 ) || 
    push!(roots,bisect_root(f,1/λ,right_point_of_range)[1]) # second intersection point

    # if the first root is zero, cut it off
    roots[1] == 0 ? shift!(roots) : [] # empty expression

    # find roots after 1/λ
    startindex = 1 + length(roots)
    if n-length(roots) > 0
        roots = [roots; zeros(n-length(roots))]
        left_point = (2*ceil(1/(π*λ)-1/2)+2)*π/2
        right_point = (2*ceil(1/(π*λ)-1/2)+3)*π/2
        for i = startindex:length(roots)
            @inbounds roots[i] = bisect_root(f,left_point,right_point)[1]
            right_point = roots[i] + π
            left_point = left_point + π
        end
    end
    return roots
end

# bissection method to find the zeros of a function in a particular interval [x1,x2]
function bisect_root{T<:AbstractFloat}(fn::Function, x1::T, x2::T)
    xm = middle(x1, x2)
    s1 = sign(fn(x1))
    s2 = sign(fn(x2))
    while x1 < xm < x2
        sm =  sign(fn(xm))

        if s1 != sm
            x2 = xm
            s2 = sm
        else
            x1 = xm
            s1 = sm
        end

        xm = middle(x1, x2)
    end

    return x1, x2
end

# helper function to find the real "mid point" of two given floating point numbers
function middle(x1::Float64, x2::Float64)
    # use the usual float rules for combining non-finite numbers
    if !isfinite(x1) || !isfinite(x2)
        return x1 + x2
    end

    # always return 0.0 when inputs have opposite sign
    if sign(x1) != sign(x2) && x1 != 0.0 && x2 != 0.0
        return 0.0
    end

    negate = x1 < 0.0 || x2 < 0.0

    x1_int = reinterpret(UInt64, abs(x1)) # TODO: find automatic conversion to UInt based on T
    x2_int = reinterpret(UInt64, abs(x2))
    unsigned = reinterpret(Float64, (x1_int + x2_int) >> 1)

    return negate ? -unsigned : unsigned
end
