# simulate(sampler, absTOL)
simulate{T<:AbstractFloat}(sampler::Sampler, absTOL::T) = simulate(sampler, absTOL=absTOL)

# simulate(sampler, absTOL, failProb)
simulate{T<:AbstractFloat}(sampler::Sampler, absTOL::T, failProb::T) = simulate(sampler, absTOL=absTOL, failProb=failProb)

# simulate(sampler, absTOL=..., relTOL=..., failProb=...)
function simulate{T<:AbstractFloat}(sampler::Sampler; absTOL::T=Inf, relTOL::T=Inf, failProb::T=0.1)

  # checks on inputs
  if absTOL < 0 || relTOL < 0
    error("supplied tolerances (absolute and/or relative) must be positive!")
  elseif failProb < 0 || failProb > 1
    error("failure probability failProb must be between 0 and 1!")
  elseif isinf(absTOL) && isinf(relTOL)
    error("supply a relative or absolute tolerance!")
  end

  t = Vector{T}()

  # run mimc or cmimc
  if sampler.continuate
    absTOL_ = 2^(sampler.nTOL-1)*absTOL
    relTOL_ = 2^(sampler.nTOL-1)*relTOL
    sampler.continuate = false # run initial hierarchy
    mimc(sampler,absTOL_,relTOL_,failProb)
    sampler.continuate = true
    for i in 1:sampler.nTOL-1
      absTOL_ = 2^(sampler.nTOL-1-i)*absTOL
      relTOL_ = 2^(sampler.nTOL-1-i)*relTOL
      push!(t, @elapsed mimc(sampler,absTOL_,relTOL_,failProb))
    end
  else
      push!(t, @elapsed mimc(sampler,absTOL,relTOL,failProb)) # plain mimc
  end

  return cumsum(t) # return timing results

end

# mimc simulation
function mimc{d,T<:AbstractFloat}(sampler::Sampler{d}, absTOL::T, relTOL::T, failProb::T)

  # print some info to screen
  if ( sampler.showInfo )
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
    @printf(sampler.ioStream,"%s","*** MultilevelEstimators.jl @$(now())\n")
    @printf(sampler.ioStream,"%s","*** Simulating $(sampler.sampleFunction)\n")
    isCont = sampler.continuate ? "" : "not"
    idxSet = ( isa(sampler.indexSet,ML) || isa(sampler.indexSet,SL) || isa(sampler.indexSet,AD) ) ? "$(sampler.indexSet)" : "$(sampler.indexSet)"[1:34]
    @printf(sampler.ioStream,"%s","*** Using a $(idxSet), $isCont continuating \n")
    @printf(sampler.ioStream,"%s",
    @sprintf("*** absTOL = %0.3e / relTOL = %0.3e (failure probability of %0.2f)\n",absTOL,relTOL,failProb))
    @printf(sampler.ioStream,"%s","--------------------------------------------------------------------------------\n")
  end

  # loop variables
  converged = false
  L = 0

  # aliases
  N = Int64
  λ = sampler.numberGenerator.λ
  q = nshifts(sampler.numberGenerator)
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1

  # variable definition
  F = zeros(T,sampler.Z); G = zeros(T,sampler.Z)
  E = Dict{Index{d,Vector{N}},Vector{T}}() # dict for expected values
  V = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variances
  Vf = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of the function
  W = Dict{Index{d,Vector{N}},T}() # dict for costs
  S = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples
  Vest = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of estimator
  C = sqrt(2)*erfcinv(failProb)
  oldindexset = Set{Index{d,Vector{N}}}()
  old = Set{Index{d,Vector{N}}}()
  profit = Dict{Index{d,Vector{N}},T}()
  profit[Index(zeros(N,d))] = 0.

  while !converged
    # get index set in current iteration
    if typeof(sampler.indexSet) <: AD
      # find index with largest gain (sparse addition)
      i = collect(keys(profit))[indmin(collect(values(profit)))]
      delete!(profit,i) # delete from profit set
      push!(old,i) # add to old set and check admissables
      for p in 1:d 
        j = copy(i)::Index{d,Vector{N}}
        j[p] += 1
        if isAdmissable(union(old, Set(keys(profit))),j)
          profit[j] = 0.
        end
      end
      newindexset = union(old, Set(keys(profit)))::Set{Index{d,Vector{N}}}
    else
      newindexset = getIndexSet(sampler.indexSet,L)::Set{Index{d,Vector{N}}}
    end

    # get boundary in current iteration
    boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
    if typeof(sampler.indexSet) == ML
      push!(boundary,Index(sampler.L*ones(N,d)))
    elseif typeof(sampler.indexSet) <: AD
      boundary = union(boundary,getBoundary(newindexset)::Set{Index{d,Vector{N}}})
    else
      boundary = union(boundary,getBoundary(getIndexSet(sampler.indexSet, L))::Set{Index{d,Vector{N}}})
    end

    # new indices to add
    indicesToAdd = setdiff( newindexset, oldindexset )

    # print some info to the screen
    !sampler.showInfo || print_with_color(:blue, sampler.ioStream, "*** currently running at L = $(L)...\n")

    # take initial number of samples
    if !sampler.continuate
      for index::Index{d,Vector{N}} in indicesToAdd
        if !haskey(sampler.times,index)
          time = @elapsed sample(sampler, sampler.Nstar, index)
          setindex!(sampler.times, time, index) # cumulative time
        end
      end
    end
    (E,Vf,V,W,Vest,profit) = updateDicts(sampler,newindexset,E,Vf,V,W,Vest,profit)

    # determine actual absolute tolerance
    minTOL = min(absTOL,relTOL*maximum(sum(values(E))))

    # estimate initial bias and splitting
    G = zero(T)
    for index::Index{d,Vector{N}} in boundary
      G += maximum(abs(E[index]))
    end
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("   *** initial bias = %0.6e \n", G ) )
    splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
    !sampler.showInfo || print_with_color(:green, sampler.ioStream, @sprintf("   *** splitting = %0.6e \n", splitting ) )

    # calculate optimal number of samples
    mySum = zeros(T,sampler.Z)
    for index::Index{d,Vector{N}} in newindexset
      mySum += ( ( Vf[index].*W[index].^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index::Index{d,Vector{N}} in newindexset
      Nopt = ( ( C/(splitting*minTOL) )^2 * 1/q * (Vf[index]./W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,max(3.,maximum(Nopt)))
    end

    # take additional samples at each level
    for index::Index{d,Vector{N}} in newindexset
      samplesToTake = haskey(sampler.samples,index) ? max( 0, S[index]-size(sampler.samples[index],dir) ) : S[index]
      if samplesToTake > 0
        time = @elapsed sample(sampler, samplesToTake, index )
        sampler.times[index] = haskey(sampler.times,index) ? sampler.times[index] + time : time # cumulative time
      end
    end
    (E,Vf,V,W,Vest,profit) = updateDicts(sampler,newindexset,E,Vf,V,W,Vest,profit)

    # safety
    F = C*sqrt(sum(values(Vest)))
    G = zero(T)
    for index::Index{d,Vector{N}} in boundary
      G += maximum(abs(E[index]))
    end
    splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
    if sampler.safety
      while maximum(F) > splitting*minTOL
        # double number of samples on index with max ratio Vest/W
        maxratio = G = zero(T)
        maxindex = Index(zeros(N,d))::Index{d,Vector{N}}
        for index::Index{d,Vector{N}} in newindexset
          ratio = maximum(Vest[index])/W[index]
          if ratio > maxratio
            maxratio = ratio
            maxindex = index::Index{d,Vector{N}}
          end
        end
        n = size(sampler.samples[maxindex],dir)
        time = @elapsed sample(sampler, nextpow2(n+1)-n, maxindex ) # round to nearest power of two
        sampler.times[maxindex] += time # cumulative time
        (E,Vf,V,W,Vest,profit) = updateDicts(sampler,newindexset,E,Vf,V,W,Vest,profit)

        # recompute statistical error and splitting
        F = C*sqrt(sum(values(Vest)))
        G = zero(T)
        for index::Index{d,Vector{N}} in boundary
          G += maximum(abs(E[index]))
        end
        splitting = ( G < minTOL/2 ) ? 1 - G/minTOL : 0.5
      end
    end

    # check for convergence
    if ( L > 1 )
      Etemp = Dict{Index{d,Vector{N}},Vector{T}}()
      Vtemp = Dict{Index{d,Vector{N}},Vector{T}}()
      for index::Index{d,Vector{N}} in newindexset # compute actual mean and variance also in continuation case
        Etemp[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
        Vtemp[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      end
      converged = ( maximum(F + G)::T < minTOL )

      # print some info to the screen
      !sampler.showInfo || print_with_color(:magenta, sampler.ioStream, 
        @sprintf("*** error estimate is %0.6e + %0.6e = %0.6e \n", maximum(F), maximum(G), maximum(F+G) ) )
      μ = maximum(sum(values(Etemp)))
      σ = sqrt(maximum(sum(values(Vtemp))))
      !sampler.showInfo || ( !converged || print_with_color(:magenta, sampler.ioStream, 
        @sprintf("*** result is %0.6e ±(%0.6e, %0.6e, %0.6e) \n", μ, σ, 2*σ, 3*σ ) ) )
    end

    # print warning if no convergence
    if !converged
      if L == sampler.maxL
        warn("maximum level reached and no convergence yet, sorry! :(\n")
        converged = true
      end
    end

    # update L
    L += 1
    oldindexset = newindexset
  end

  #
  # PRINT OVERVIEW TABLE
  #
  if ( sampler.showInfo )
    @printf(sampler.ioStream,"%s","-------------------------------------------------------------------------------- \n")
    itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
    @printf(sampler.ioStream,"%s","  "*itype*"       E              V               N               W               \n")
    @printf(sampler.ioStream,"%s","-------------------------------------------------------------------------------- \n")
    for index::Index{d,Vector{N}} in sort(Set(collect(keys(E))))
      str = ("  $(index.indices)            "[1:13])
      str *= @sprintf("%12.5e",maximum(E[index]))
      str *= @sprintf("    %0.6e",maximum(Vf[index]))
      str *= @sprintf("    %d               ",prod(size(sampler.samples[index])[1:2]))[1:16]
      str *= @sprintf("    %0.6e \n",W[index])
      @printf(sampler.ioStream,"%s",str)
    end
  end
end

# refactored function that updates E/Vf/V/W/Vest/profit from current sampler
function updateDicts{d,N,T}(sampler::Sampler, indices::Set{Index{d,Vector{N}}}, 
  E::Dict{Index{d,Vector{N}},Vector{T}}, Vf::Dict{Index{d,Vector{N}},Vector{T}}, 
    V::Dict{Index{d,Vector{N}},Vector{T}}, W::Dict{Index{d,Vector{N}},T}, 
      Vest::Dict{Index{d,Vector{N}},Vector{T}}, profit::Dict{Index{d,Vector{N}},T})

  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1

  # update all
  for index::Index{d,Vector{N}} in indices
    if haskey(sampler.samples,index)
      E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
      Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
      nsamples = size(sampler.samples[index],dir)
      W[index] =  sampler.useTime ? sampler.times[index]/nsamples : prod(2.^(sampler.γ.*index))
      Vest[index] = V[index]/nsamples
      if typeof(sampler.indexSet) <: AD && haskey(profit,index) # compute gains when adaptive
        profit[index] = maximum(abs(E[index])./(Vf[index].*W[index]))
      end
    end
  end

  # correction when doing continuation
  if sampler.continuate
    (A,α,B,β,Bacc,βacc,C,γ) = estimateProblemParameters(sampler)  # ! fits only max E_\ell, V_\ell
    Vftilde = bayesianUpdateVariance(sampler,A,α,B,β,E,Vf)
    Vtilde = bayesianUpdateVariance(sampler,A,α,Bacc,βacc,E,V)
    for index::Index{d,Vector{N}} in indices # use model fit on all indices
      E[index] = [A*prod(2.^(-α.*index))]
      Vf[index] = haskey(Vftilde,index) ? [Vftilde[index]] : [B*prod(2.^(-β.*index))]
      V[index] = haskey(Vtilde,index) ? [Vtilde[index]] : [Bacc*prod(2.^(-βacc.*index))]
      W[index] = C*prod(2.^(γ.*index))
      Vest[index] = haskey(sampler.samples,index) ? V[index]/size(sampler.samples[index],dir) : [zero(T)]
      if typeof(sampler.indexSet) == AD && haskey(profit,index)
        profit[index] = abs(E[index])./sqrt(Vf[index].*W[index])
      end
    end
  end

  return (E::Dict{Index{d,Vector{N}},Vector{T}},Vf::Dict{Index{d,Vector{N}},Vector{T}},
    V::Dict{Index{d,Vector{N}},Vector{T}},W::Dict{Index{d,Vector{N}},T},
      Vest::Dict{Index{d,Vector{N}},Vector{T}},profit::Dict{Index{d,Vector{N}},T})
end

# estimate problem parameters A, α, B, β, Bacc, βacc, C, γ and Vtilde in the continuation problem
function estimateProblemParameters{d}(sampler::Sampler{d})

  # for now, let A, α, B, β, Bacc, βacc, C and γ follow from a fit through the available E_\ell, V_\ell and W_\ell
  E = Dict{Index{d,Vector{Int64}},Float64}()
  Vf = Dict{Index{d,Vector{Int64}},Float64}()
  V = Dict{Index{d,Vector{Int64}},Float64}()
  W = Dict{Index{d,Vector{Int64}},Float64}()
  for index in keys(sampler.samples)
    E[index] = maximum(squeeze(mean(sampler.samples[index],(1,2)),(1,2)))
    Vf[index] = maximum(squeeze(mean(var(sampler.samples[index],2),1),(1,2)))
    V[index] = maximum(squeeze(var(mean(sampler.samples[index],1),2),(1,2)))
    W[index] = sampler.useTime ? sampler.times[index] : prod(2.^(index.*sampler.γ))
  end

  # least-squares fit
  (A,α) = leastSquaresFit(E)
  (B,β) = leastSquaresFit(Vf)
  (Bacc,βacc) = leastSquaresFit(V)
  (C,γ) = leastSquaresFit(W)

  # println("I have fitted A = $(A) and α= $(α)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       E              E_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(E))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(E[index])))
  #   str *= @sprintf("    %0.6e",A*prod(2.^(-α.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(A*prod(2.^(-α.*index.indices))-maximum(abs(E[index])))/(maximum(abs(E[index]))))
  #   println(str)
  # end

  # println("I have fitted B = $(B) and β= $(β)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       Vf             Vf_model        error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(Vf))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(Vf[index])))
  #   str *= @sprintf("    %0.6e",B*prod(2.^(-β.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(B*prod(2.^(-β.*index.indices))-maximum(abs(Vf[index])))/(maximum(abs(Vf[index]))))
  #   println(str)
  # end

  # println("I have fitted Bacc = $(Bacc) and βacc = $(βacc)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       V              V_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(V))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(V[index])))
  #   str *= @sprintf("    %0.6e",Bacc*prod(2.^(-βacc.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(Bacc*prod(2.^(-βacc.*index.indices))-maximum(abs(V[index])))/(maximum(abs(V[index]))))
  #   println(str)
  # end

  # println("I have fitted C = $(C) and γ= $(γ)")
  # println("--------------------------------------------------------------------------------")
  # itype = ndims(sampler.indexSet) == 1 ? "level" : "index"
  # println("  "*itype*"       W              W_model         error                  ")
  # println("--------------------------------------------------------------------------------")
  # for index in sort(Set(collect(keys(W))))
  #   str = ("  $(index.indices)            "[1:13])
  #   str *= @sprintf("%12.5e",maximum(abs(W[index])))
  #   str *= @sprintf("    %0.6e",C*prod(2.^(γ.*index.indices)))
  #   str *= @sprintf("    %0.6e",abs(C*prod(2.^(γ.*index.indices))-maximum(abs(W[index])))/(maximum(abs(W[index]))))
  #   println(str)
  # end

  # # index set from which to derive the parameters
  # L = max(1,sampler.L-sampler.bayesianProperties.ell_0)
  # I = setdiff(getIndexSet(sampler.indexset,sampler.L),getIndexSet(sampler.indexset,L))
  # numerator = 0.
  # denominator = 0.
  # for i in I
  #   numerator += prod(2^(-(α-β).*i.indices))*N[i]*E[i]
  #   denominator += N[i]*prod(2^(-(2*α-β).*i.indices))
  # end
  # A_ = numerator/denominator

  # numerator = 0.
  # denominator = 0.
  # for i in I
  #   numerator += prod(2^(β.*i.indices))*      N[i]*E[i]
  #   denominator += N[i]*prod(2^(-(2*α-β).*i.indices))
  # end
  # B_ = numerator/denominator

  return (A::Float64, α::Vector{Float64}, B::Float64, β::Vector{Float64}, Bacc::Float64, βacc::Vector{Float64}, C::Float64, γ::Vector{Float64})
end

# Bayesian update of variances
function bayesianUpdateVariance{d,N1<:Integer,T<:AbstractFloat}(sampler::Sampler,A::T,α::Vector{T},B::T,β::Vector{T}, E::Dict{Index{d,Vector{N1}},Vector{T}}, V::Dict{Index{d,Vector{N1}},Vector{T}})
  dir = isa(sampler.numberGenerator,MCgenerator) ? 2 : 1
  Vtilde = Dict{Index{d,Vector{N1}},T}()

  for index in keys(E)
    N = size(sampler.samples[index],dir)
    μ = A*prod(2.^(-α.*index))
    λ = 1/B*prod(2.^(β.*index))
    Γ3 = sampler.k[2]*λ + N/2
    Γ4 = sampler.k[2] + 0.5*(N-1)*maximum(V[index]) + sampler.k[1]*N*(maximum(E[index])-μ)^2/(2*(sampler.k[1]+N))
    Vtilde[index] = Γ4/Γ3
  end
  zero_idx = Index(zeros(N1,d))
  Vtilde[zero_idx] = maximum(V[zero_idx]) # ! correction for level / index 0

  return Vtilde::Dict{Index{d,Vector{N1}},T}
end

# do a least-squares fit of the indices and data provided in the dict
function leastSquaresFit{d,N<:Integer,T<:AbstractFloat}(dict::Dict{Index{d,Vector{N}},T})
  keyz = sort(Set(collect(keys(dict))))
  m = length(keyz)
  data = zeros(m-1,d+1)
  rhs = zeros(m-1)
  for i in 1:m-1 # ! correction for level / index 0
    data[i,1:d] = keyz[i+1].indices
    data[i,d+1] = 1
    rhs[i] = log2(abs(dict[keyz[i+1]][1]))
  end
  coeffs = data\rhs # least-squares system
  X = 2^coeffs[d+1]
  ξ = abs(coeffs[1:d])

  return (X,ξ)
end
