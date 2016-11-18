# abstract Gaussian field sampler
abstract GaussianFieldSampler

# Empty sampler type
type EmptySampler <: GaussianFieldSampler

end

# KL expansion
type KLExpansion{d,N<:Integer,A<:AbstractVector,B<:AbstractVector,C<:AbstractVector} <: GaussianFieldSampler
  mkl::N # number of terms in KL expansion
  eigenval::A # d-dimensional KL eigenvalues (Array of tupples with index and value)
  eigenfunc::B # 1d KL eigenfunctions
  ad::Bool # true if adaptive in number of terms
  s::C
end

# utilities
ndims{d}(G::KLExpansion{d}) = d
inttype{d,N}(G::KLExpansion{d,N}) = N
floattype(G::KLExpansion) = eltype(eltype(G.eigenfunc))

# methods
isValid(G::EmptySampler) = true

function show(io::IO, G::EmptySampler)
  print(io, "emtyp Gaussian field sampler")
end

function isValid(G::KLExpansion)
  d = ndims(G)
  N = inttype(G)
  T = floattype(G)
  return G.mkl > 0 && typeof(G.eigenval) == Vector{Tuple{Index{d,Vector{N}},T}} && typeof(G.eigenfunc) == Vector{Array{T,2}}
end

function show(io::IO, G::KLExpansion)
  str = G.ad ? "level-adaptive " : ""
  print(io, str*"KL-expansion with $(G.mkl) terms")
end

# covariance kernel
abstract Kernel

type MaternKernel{T<:AbstractFloat} <: Kernel
  λ::T # correlation length
  σ::T # variance
  ν::T # smoothness parameter
  p::T # p-norm

  function MaternKernel(λ, σ, ν, p)
    λ > 0 || error("correlation length of random field cannot be negative or zero!")
    σ > 0 || error("variance of random field cannot be negative or zero!")
    ν > 0 || error("smoothness of random field cannot be negative or zero!")
    p >= 1 || error("p-norm needs p>1!")

    return new(λ, σ, ν, p)
  end
end

MaternKernel{T<:AbstractFloat}(λ::T,σ::T,ν::T,p::T) = MaternKernel{T}(λ,σ,ν,p)

# methods
isValid(M::MaternKernel) = ( λ > 0 && σ > 0 && ν > 0 && p >= 1 )

function applyKernel{T<:AbstractFloat}(K::MaternKernel,x::AbstractArray{T},y::AbstractArray{T})
  cov = zeros(T,size(x.-y))
  for i in 1:length(x)
    for j in 1:length(y)
      @inbounds cov[i,j] = K.σ^2*2^(1-K.ν)/gamma(K.ν)*(sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ).^K.ν.*besselk(K.ν,sqrt(2*K.ν)*norm(x[i]-y[j],K.p)/K.λ)
    end
  end
  cov[(x.-y).==zero(T)]=one(T)

  return cov
end

function show(io::IO, M::MaternKernel)
  print(io, "Matern kernel with correlation length λ = $(M.λ), variance σ^2 = $((M.σ)^2), smoothness ν = $(M.ν) and $(M.p)-norm.")
end

#
# Main methods for creating and composing the KL expansion
#

# constructor for separable kernels (non-separable is not implemented yet...)
function KLExpansion{T,N}(kernel::Kernel, d::N, mkl::N; m0::N = 4, maxL::N = 10, ad::Bool = false, x::Vector{Vector{T}} = Vector{Vector{T}}(), s::Vector{N} = Vector{N}())
  d > 0 || error("dimension cannot be negative or zero!")
  mkl > 0 || error("number of KL-terms must be a positive integer!")
  maxL >= 0 || error("maximum indexset indicator must be postitive")
  isempty(x) || ( length(x) == maxL+1 || error("supply as many points at each level as maxL+1") )
  if isempty(s)
    s = mkl*ones(maxL+1)
  end

  # calculate eigenvalues
  if typeof(kernel) <: MaternKernel && kernel.ν == 0.5
    ω = findroots(kernel.λ, mkl+1)
    θ = 2*kernel.σ^2*kernel.λ./(kernel.λ^2*ω.^2+1)
  else
    θ,dummy = nystrom(kernel,0,mkl+1)
  end

  # compose the d-dimensional eigenvalues by adaptive search
  O = Vector{Tuple{Index{d,Vector{N}},T}}() # could be oredered dict as well
  A = Dict{Index{d,Vector{N}},T}()
  A[one(Index{d,Vector{N}})] = θ[1]^d
  converged = false
  max_mkl = 1 # store maximum number of 1d eigenvalues
  while length(O) < mkl
    idx = collect(keys(A))[indmax(values(A))] # find minimum of values
    max_mkl = maximum([max_mkl maximum(idx)])
    push!(O,(idx,A[idx]))
    for p = 1:d # check each new admissable index
      idx = copy(idx)
      idx[p] = idx[p] + 1
      if ~haskey(A,idx) # add to admissable set, if not already added
        A[copy(idx)] = prod(θ[idx.indices])
      end
      idx = copy(idx)
      idx[p] = idx[p] - 1
    end
    delete!(A,idx)
  end

  # calculate 1d eigenfunctions
  eigenfunc = Array{T,2}[]
  if typeof(kernel) <: MaternKernel && kernel.ν == 0.5
    ω = ω[1:max_mkl] # cut off ω
    n = sqrt(2)/2*sqrt(1./ω.*(kernel.λ^2*ω.^2.*cos(ω).*sin(ω)+kernel.λ^2*ω.^3-2*kernel.λ*ω.*cos(ω).^2-cos(ω).*sin(ω)+ω)+2*kernel.λ)
    for L in 0:maxL
      if isempty(x)
        m = m0*2^L
        pts = collect(1/2/m:1/m:1-1/2/m) # finite volume method, cell centered values
      else
        pts = collect(x[L+1]) # other values if specified
      end
      push!(eigenfunc,diagm(1./n)*( sin(ω*pts') + kernel.λ*diagm(ω)*cos(ω*pts') ))
    end
  else
    for L in 0:maxL
      if isempty(x)
        m = m0*2^L
      else
        m = length(x[L+1])
      end
      dummy,temp = nystrom(kernel,m,max_mkl) # HACK for now use points on [0,1]
      push!(eigenfunc, temp)
    end
  end

  return KLExpansion{d,N,typeof(O),typeof(eigenfunc),typeof(s)}(mkl,O,eigenfunc,ad,s)
end

# function to compose KL expansion
function compose{K<:KLExpansion,T<:AbstractFloat,t,V}(kl::K, xi::Vector{T}, index::Index{t,V})
  N = inttype(kl)
  d = ndims(kl)
  nterms = length(kl.eigenval)
  if kl.ad
    nterms = kl.s[index.indices[end]+1]
    index_ = Index(index.indices[1:end-1])::Index{t-1,V}
  elseif ndims(index) < d # FIX for multilevel case
    index_ = Index(repeat([index[1]],inner=[d]))::Index{d,V}
  else
    index_ = index
  end

  # first term for type-stability of k
  v = kl.eigenfunc[index_[1]+1][kl.eigenval[1][1][1],:]::Array{T,2}
  for p = 2:d
    v = kron(kl.eigenfunc[index_[p]+1][kl.eigenval[1][1][p],:],v)::Array{T,2}
  end
  k = xi[1]*sqrt(kl.eigenval[1][2])*v

  # rest of the terms
  for i in 2:nterms
    v = kl.eigenfunc[index_[1]+1][kl.eigenval[i][1][1],:]::Array{T,2}
    for p = 2:d
      v = kron(kl.eigenfunc[index_[p]+1][kl.eigenval[i][1][p],:],v)::Array{T,2}
    end
    k += xi[i]*sqrt(kl.eigenval[i][2])*v
  end
  return k
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

# helper function to compute nodes and weights of Gauss quadrature on [0,1]
function gauss_legendre_0_1(m)
  nodes, weights = gausslegendre(m)
  nodes = 0.5 + nodes/2 # scale roots
  weights = 0.5*weights # scale weights

  return nodes, weights
end

# nystrom method
function nystrom{N<:Integer}(kernel::Kernel,m::N,nterms::N)
  # compute covariance matrix
  nodes, weights = gauss_legendre_0_1(nterms)
  K = applyKernel(kernel,nodes,nodes')

  # obtain eigenvalues and eigenfunctions
  D = diagm(weights)
  Dsqrt = sqrt(D)
  Z = Symmetric(triu(Dsqrt*K*Dsqrt))
  eigenval, eigenfunc = eig(Z)
  #eigenval = real(eigenval)
  idx = sortperm(eigenval,rev=true)
  sort!(eigenval,rev=true)
  eigenfunc = eigenfunc[:,idx]

  f = weights.*(diagm(sqrt(1./weights))*eigenfunc)
  lambda = 1./eigenval

  # nystrom method
  eigenfunc = zeros(nterms,m)
  x = 1/2/m:1/m:1-1/2/m
  K = applyKernel(kernel,x,nodes')
  for j in 1:nterms
    @inbounds eigenfunc[j,:] = lambda[j]*K*f[:,j]
  end

  return eigenval, eigenfunc
end
