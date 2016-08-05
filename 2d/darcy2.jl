# sample function for parametrized PDE with point evaluation as QoI
function parametrizedPDEpointEvaluation{T<:AbstractFloat,d,S<:Settings}(xi::Vector{T},j::Index{d},settings::S)
  kl = compose(settings.gaussianFieldSampler,xi,j)::Matrix{T} # KL expansion
  m = settings.m0.*2.^j.indices # grid sizes
  mx, my = length(j.indices) < 2 ? tuple(repeat(m,inner=[2])...) : tuple(m...)
  k = exp(reshape(kl,(mx,my)))::Matrix{T} # take exponential
  p = epde2(k)::Matrix{T} # solve the deterministic PDE (flow cell geometry)
  vx = 1/2/mx:1/mx:1-1/2/mx
  vy = 1/2/my:1/my:1-1/2/my
  myInterpolator = interpolate((vx,vy), p, Gridded(Linear()))
  z = [0.25,0.5,0.75]
  return myInterpolator[z,z]::Matrix{T}
end

# sample function for parametrized PDE with effective conductivity as QoI
function parametrizedPDEEffectiveConductivity{T<:AbstractFloat,d,S<:Settings}(xi::Vector{T},j::Index{d},settings::S)
  kl = compose(settings.gaussianFieldSampler,xi,j)::Matrix{T} # KL expansion
  mx, my = length(j.indices) < 2 ? tuple(repeat(settings.m0.*2.^j.indices,inner=[2])...) : tuple(settings.m0.*2.^j.indices...)
  k = exp(reshape(kl,(mx,my)))::Matrix{T} # take exponential
  p = epde2d(k)::Matrix{T} # solve the deterministic PDE (flow cell geometry)
  return trapz(2.0*mx*squeeze(k[mx,:],1).*squeeze(p[mx,:],1),my)::T
  #return trapz(1/3.0*mx*squeeze(k[mx,:],1).*(9*squeeze(p[mx,:],1)-squeeze(p[mx-1,:],1)),my) # 2nd order approx
end

#
# helper functions
#

# trapezoidal rule in 1D
function trapz{T<:Real,N<:Integer}(y::Array{T,1}, mx::N)
  m = length(y)
  r = 0.0
  for i in 2:m
    r += y[i] + y[i-1]
  end
  return r/mx/2.0
end