# MIMC simulation
function simulate{d,T,N}(sampler::Sampler{d,T,N}, TOL::T)
  # checks on inputs
  @assert TOL ≥ 0

  # failure probability 
  epsilon = 0.1

  # algorithm details, loop variables etc.
  Nstar = isa(sampler.settings.numberGenerator,MCgenerator) ? 16 : 16

  # Loop variables
  converged = false
  sampler.K = 0

  A = 0.; B = 0. # variable definition for use outside while-loop
  E     = Dict{Index{d},Array{T,1}}() # dict for expected values
  V     = Dict{Index{d},Array{T,1}}() # dict for variances
  Vf    = Dict{Index{d},Array{T,1}}() # dict for variance of the function
  W     = Dict{Index{d},T}() # dict for costs
  S     = Dict{Index{d},N}() # dict for optimal number of samples
  Vest  = Dict{Index{d},Array{T,1}}() # dict for variance of estimator
  C     = sqrt(2)*erfcinv(epsilon)
  λ     = sampler.settings.numberGenerator.λ
  q     = nshifts(sampler.settings.numberGenerator)
  splitting = sampler.settings.splitting
  dir = isa(sampler.settings.numberGenerator,MCgenerator) ? 2 : 1
  oldindexset = Set{Index{d}}()

  old = Set{Index{d}}()
  profit = Dict{Index{d},T}()
  profit[zero(Index{d})] = 0.

  while !converged
    # get index set in current iteration
    if kind(sampler) == AD
      # find index with largest gain
      i = collect(keys(profit))[indmin(collect(values(profit)))]
      delete!(profit,i) # delete from profit set
      push!(old,i) # add to old set and check admissables
      # full addition
      FTSet = createIndexSet(FT,d)
      FTSet = getIndexSet(FTSet,1)
      for p in FTSet
        j = copy(p+i)
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
      newindexset = union(old, Set(keys(profit)))
    else
      newindexset = getIndexSet(sampler.settings.indexset,sampler.K)
    end

    # get boundary in current iteration
    if kind(sampler) == ML
      boundary = Set{Index{d}}()
      push!(boundary,Index(sampler.K*ones(N,d)))
    elseif kind(sampler) == AD
      boundary = getBoundary(newindexset)
    else
      boundary = getBoundary(sampler.settings.indexset, sampler.K)
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
    B = 0.
    for index in boundary
      B += abs(E[index])
    end
    !sampler.settings.showInfo || print_with_color(:green, "   >>> B = $(maximum(B))\n")
    splitting = ( maximum(B) < TOL/2 ) ? 1 - maximum(B)/TOL : 0.5
    !sampler.settings.showInfo || print_with_color(:green, "   >>> splitting = $(splitting)\n")

    # calculate optimal number of samples
    mySum = 0.
    for index in newindexset
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
      maxindex = NaN
      for index in newindexset
        ratio = maximum(Vest[index])/W[index] # shorter version not guaranteed when V/W?
        if ratio > maxratio
          maxratio = ratio
          maxindex = index
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
      B = 0.
      for index in boundary
        B += abs(E[index])
      end
      Etemp = sum(values(E))
      Vtemp = sum(values(Vf))
      converged = ( maximum(A + B) < TOL )

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

  return (sum(values(E)), sum(values(Vf)), A, B, splitting, S)
  
end
