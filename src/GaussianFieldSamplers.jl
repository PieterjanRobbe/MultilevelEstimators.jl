# abstract Gaussian field sampler
abstract GaussianFieldSampler

# KL expansion
type KLexpansion{d,T,N} <: GaussianFieldSampler
  mkl::N # number of terms in KL expansion
  eigenval::Array{Array{T,1},1} # d-dimensional KL eigenvalues
  eigenfunc::Array{Array{T,2},1} # 1d KL eigenfunctions
end

function createKLexpansion{T,N}(d::N, λ::T, σ::T, ν::T, mkl::N; m0::N = 4, maxK::N = 15, cov::Function = nothing)
  d > 0 || error("dimension cannot be negative or zero!")
  λ > 0 || error("correlation length of random field cannot be negative or zero!")
  σ > 0 || error("variance of random field cannot be negative or zero!")
  ν > 0 || error("smoothness of random field cannot be negative or zero!")
  mkl > 0 || error("number of KL-terms must be a positive integer!")

  # calculate eigenvalues
  if ν == 0.5
    ω = findroots(λ, mkl+1)
    θ = 2*σ^2*λ./(λ^2*ω.^2+1)
  else
    θ,dummy = nystrom(cov,0,mkl+1)
  end

  # compose the d-dimensional eigenvalues by adaptive search
  O = Array{T,1}[] #Dict{Index{d},T}()
  A = Dict{Index{d},T}()
  A[one(Index{d})] = θ[1]^d
  converged = false
  max_mkl = 1 # store maximum number of 1d eigenvalues
  while length(O) < mkl
    idx = collect(keys(A))[indmax(values(A))] # find minimum of values
    max_mkl = maximum([max_mkl maximum(idx)])
    push!(O,[idx.indices;A[idx]])
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
  ev = zeros(length(O))
  for i in 1:mkl
    ev[i] = O[i][end]
  end

  # calculate 1d eigenfunctions
  if ν == 0.5
    ω = ω[1:max_mkl] # cut off ω
    eigenfunc = Array{T,2}[]
    n = sqrt(2)/2*sqrt(1./ω.*(λ^2*ω.^2.*cos(ω).*sin(ω)+λ^2*ω.^3-2*λ*ω.*cos(ω).^2-cos(ω).*sin(ω)+ω)+2*λ)
    for K in 0:maxK
      m = m0*2^K
      x = collect(1/2/m:1/m:1-1/2/m)
      push!(eigenfunc,diagm(1./n)*( sin(ω*x') + λ*diagm(ω)*cos(ω*x') ))
    end
  else
    eigenfunc = Array{T,2}[]
    for K in 0:maxK
      dummy,temp = nystrom(cov,m0*2^K,max_mkl)
      push!(eigenfunc, temp)
    end
  end

  return KLexpansion{d,T,N}(mkl,O,eigenfunc)
end

# function to compose KL expansion
function compose{d,T,N}(kl::KLexpansion{d,T,N}, xi::Vector{T}, index::Index)
  index = ndims(index) < d ? Index(repeat([index[1]],inner=[d])) : index # FIX for multilevel case v
  k = 0
  for i in 1:length(kl.eigenval)
    v = kl.eigenfunc[index[1]+1][convert(N,kl.eigenval[i][1]),:]
    for p = 2:d
      v = kron(kl.eigenfunc[index[p]+1][convert(N,kl.eigenval[i][p]),:],v)
    end
    k += xi[i]*sqrt(kl.eigenval[i][d+1])*v
  end
  return k
end

# find all positive (>0) zeros of the transcendental function tan(ω) = 2*λ*ω/(λ^2*ω^2-1)
function findroots{T<:Real,N<:Integer}(λ::T, n::N)
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
      roots[i] = bisect_root(f,left_point,right_point)[1]
      right_point = roots[i] + π
      left_point = left_point + π
    end
  else
    roots = zeros(0)
  end

  # find roots inside range around 1/λ
  ( length(roots) ≥ n || floor(1/(π*λ)-1/2) < 0 ) || 
    push!(roots,bisect_root(f,left_point_of_range+eps(Float64),1/λ)[1]) # first intersection point
  ( length(roots) ≥ n || ceil(1/(π*λ)-1/2) < 0 ) || 
    push!(roots,bisect_root(f,1/λ,right_point_of_range)[1]) # second intersection point

  # if the first root is zero, cut it off
  roots[1] == 0 ? shift!(roots) : "" # empty expression

  # find roots after 1/λ
  startindex = 1 + length(roots)
  if n-length(roots) > 0
    roots = [roots; zeros(n-length(roots))]
    left_point = (2*ceil(1/(π*λ)-1/2)+2)*π/2
    right_point = (2*ceil(1/(π*λ)-1/2)+3)*π/2
    for i = startindex:length(roots)
      roots[i] = bisect_root(f,left_point,right_point)[1]
      right_point = roots[i] + π
      left_point = left_point + π
    end
  end
  return roots
end

# bissection method to find the zeros of a function in a particular interval [x1,x2]
function bisect_root{T<:Real}(fn::Function, x1::T, x2::T)
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

# helper function to find the real "mid point" of two given Float64 numbers
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

  x1_int = reinterpret(UInt64, abs(x1))
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

# matern kernel
function matern{T<:AbstractFloat}(λ::T,σ::T,ν::T,p::T,x::AbstractArray{T},y::AbstractArray{T})
  cov = σ^2*2^(1-ν)/gamma(ν)*(sqrt(2*ν)*norm(x.-y,1)/λ).^ν.*besselk(ν,sqrt(2*ν)*norm(x.-y,1)/λ)
  cov[(x.-y).==zero(T)]=one(T)

  return cov
end

# nystrom method
function nystrom{N<:Integer}(covariance::Function,m::N,nterms::N)
  # compute covariance matrix
  nodes, weights = gauss_legendre_0_1(nterms)
  K = covariance(nodes,nodes')

  # obtain eigenvalues and eigenfunctions
  D = diagm(weights)
  Dsqrt = sqrt(D)
  eigenval, eigenfunc = eig(Dsqrt*K*Dsqrt)
  eigenval = real(eigenval)
  idx = sortperm(eigenval,rev=true)
  sort!(eigenval,rev=true)
  eigenfunc = eigenfunc[:,idx]

  f = weights.*(diagm(sqrt(1./weights))*eigenfunc)
  lambda = 1./eigenval

  # nystrom method
  eigenfunc = zeros(nterms,m)
  x = 1/2/m:1/m:1-1/2/m
  K = covariance(x,nodes')
  for j in 1:nterms
    @inbounds eigenfunc[j,:] = lambda[j]*K*f[:,j]
  end

  return eigenval, eigenfunc
end