# MIMC simulation
function simulate{d,N<:Integer,T<:AbstractFloat}(sampler::Sampler{d,N}, TOL::T)

  # checks on inputs
  @assert TOL ≥ 0

  # failure probability 
  epsilon = 0.1

  # algorithm details, loop variables etc.
  Nstar = isa(sampler.settings.numberGenerator,MCgenerator) ? 16 : 16

  # Loop variables
  converged = false
  sampler.K = 0

  A = zeros(T,sampler.settings.Z); B = zeros(T,sampler.settings.Z) # variable definition for use outside while-loop
  E     = Dict{Index{d,Vector{N}},Vector{T}}() # dict for expected values
  V     = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variances
  Vf    = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of the function
  W     = Dict{Index{d,Vector{N}},T}() # dict for costs
  S     = Dict{Index{d,Vector{N}},N}() # dict for optimal number of samples
  Vest  = Dict{Index{d,Vector{N}},Vector{T}}() # dict for variance of estimator
  C     = sqrt(2)*erfcinv(epsilon)
  λ     = sampler.settings.numberGenerator.λ
  q     = nshifts(sampler.settings.numberGenerator)
  splitting = sampler.settings.splitting
  dir = isa(sampler.settings.numberGenerator,MCgenerator) ? 2 : 1
  oldindexset = Set{Index{d,Vector{N}}}()
  old = Set{Index{d,Vector{N}}}()
  profit = Dict{Index{d,Vector{N}},T}()
  profit[Index(zeros(N,d))] = 0.

  while !converged
    # get index set in current iteration
    if kind(sampler) == AD
      # find index with largest gain
      i = collect(keys(profit))[indmin(collect(values(profit)))]
      delete!(profit,i) # delete from profit set
      push!(old,i) # add to old set and check admissables
      # full addition
      full = createIndexSet(FT,d)::IndexSet{d,Vector{T}}
      FTSet = getIndexSet(full,1)::Set{Index{d,Vector{N}}}
      for p in FTSet
        j = copy(p+i)::Index{d,Vector{N}}
        if isAdmissable(union(old, Set(keys(profit))),j)
          profit[j] = 0.
        end
      end
      # sparse addition
      # for p in 1:d
      #   j = copy(i)
      #   j[p] += 1
      #   if isAdmissable(union(old, Set(keys(profit))),j)
      #     profit[j] = 0.
      #   end
      # end
      newindexset = union(old, Set(keys(profit)))::Set{Index{d,Vector{N}}}
    else
      newindexset = getIndexSet(sampler.settings.indexset,sampler.K)::Set{Index{d,Vector{N}}}
    end

    # get boundary in current iteration
    boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
    if kind(sampler) == ML
      push!(boundary,Index(sampler.K*ones(N,d)))
    elseif kind(sampler) == AD
      boundary = union(boundary,getBoundary(newindexset)::Set{Index{d,Vector{N}}})
    else
      boundary = union(boundary,getBoundary(sampler.settings.indexset, sampler.K)::Set{Index{d,Vector{N}}})
    end

    # new indices to add
    indicesToAdd = setdiff( newindexset, oldindexset )

    # print some info to the screen
    !sampler.settings.showInfo || print_with_color(:blue, "*** currently running at K = $(sampler.K)...\n")

    # take initial number of samples
    for index in indicesToAdd
      time = @elapsed sample(sampler, Nstar, index)
      E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
      V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
      Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      W[index] =  sampler.settings.useTime ? time/16 : prod(sampler.settings.m0.*2.^(index.indices))^sampler.settings.gamma
      Vest[index] = V[index]/Nstar
      if kind(sampler) == AD && haskey(profit,index) # compute gains when adaptive
        profit[index] = maximum(abs(E[index])./(Vf[index].*W[index]))
      end
    end

    # estimate bias and splitting parameter
    B = zeros(T,sampler.settings.Z)
    for index in boundary
      B += abs(E[index])
    end
    !sampler.settings.showInfo || print_with_color(:green, "   >>> B = $(maximum(B))\n")
    splitting = ( maximum(B) < TOL/2 ) ? 1 - maximum(B)/TOL : 0.5
    !sampler.settings.showInfo || print_with_color(:green, "   >>> splitting = $(splitting)\n")

    # calculate optimal number of samples
    mySum = zeros(T,sampler.settings.Z)
    for index::Index{d,Vector{N}} in newindexset
      mySum += ( ( Vf[index].*W[index].^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index in newindexset
      Nopt = ( ( C/(splitting*TOL) )^2 * 1/q * (Vf[index]./W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,maximum(Nopt))
    end

    # take additional samples at each level
    for index in newindexset
      samplesToTake = max( 0, S[index]-size(sampler.samples[index],dir) )
      if samplesToTake > 0
        sample(sampler, samplesToTake, index )
        E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
        V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
        Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
        Vest[index] = V[index]/length(sampler.samples[index])
        if kind(sampler) == AD && haskey(profit,index) # compute gains when adaptive
        profit[index] = maximum(abs(E[index])./(Vf[index].*W[index]))
        end
      end
    end

    # safety
    A = C*sqrt(sum(values(Vest)))
    while maximum(A) > splitting*TOL # 0.9 is safety
      # double number of samples on index with max ratio V/W
      maxratio = 0.
      maxindex = Index(zeros(N,d))
      for index in newindexset
        ratio = maximum(Vest[index])/W[index] # shorter version not guaranteed when V/W?
        if ratio > maxratio
          maxratio = ratio
          maxindex = index::Index{d,Vector{N}}
        end
      end
      n = size(sampler.samples[maxindex],dir)
      sample(sampler, nextpow2(n+1)-n, maxindex ) # round to nearest power of two
      E[maxindex] = squeeze(mean(sampler.samples[maxindex],(1,2)),(1,2))
      V[maxindex] = squeeze(var(mean(sampler.samples[maxindex],1),2),(1,2))
      Vf[maxindex] = squeeze(mean(var(sampler.samples[maxindex],2),1),(1,2))
      Vest[maxindex] = V[maxindex]/length(sampler.samples[maxindex])
      if kind(sampler) == AD && haskey(profit,maxindex) # compute gains when adaptive
        profit[maxindex] = maximum(abs(E[maxindex])./(Vf[maxindex].*W[maxindex]))
      end
      A = C*sqrt(sum(values(Vest)))
    end

    # check for convergence
    if ( sampler.K > 0 )
      B = zeros(T,sampler.settings.Z)
      for index in boundary
        B += abs(E[index])
      end
      Etemp = sum(values(E))
      Vtemp = sum(values(Vf))
      converged = ( maximum(A + B)::T < TOL )

      # print some info to the screen
      !sampler.settings.showInfo || print_with_color(:magenta, "*** error estimate is $(maximum(A)) + $(maximum(B)) = $(maximum(A+B))\n")
      !sampler.settings.showInfo || ( !converged || print_with_color(:magenta, "*** result is $(maximum(Etemp)) with a variance of $(maximum(Vtemp))\n" ) )
    end

    # print warning if no convergence
    if !converged
      if ( kind(sampler) != AD && sampler.K == sampler.settings.maxK )
        warn("maximum level reached and no convergence yet, sorry! :(\n")
      end
    end

    # update K and get new indices to add
    sampler.K += 1
    oldindexset = newindexset
  end

  return (sum(values(E))::Vector{T}, sum(values(Vf))::Vector{T}, A::Vector{T}, B::Vector{T}, splitting::T, S::Dict{Index{d,Vector{N}},N})
  
end
