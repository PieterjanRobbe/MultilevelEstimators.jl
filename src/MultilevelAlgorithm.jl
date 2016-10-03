# MIMC simulation
function simulate{d,N<:Integer,T<:AbstractFloat}(sampler::Sampler{d,N}, TOL::T, epsilon::T)

  # checks on inputs
  @assert TOL ≥ 0
  @assert epsilon ≥ 0

  # print some info to screen
  if ( sampler.settings.showInfo )
    println("--------------------------------------------------------------------------------")
    println("*** MultilevelEstimators.jl @$(now())")
    println("*** Simulating $(sampler.settings.sampleFunction)")
    println(@sprintf("*** TOL = %0.3e (failure probability of %0.2f)",TOL,epsilon))
    println("--------------------------------------------------------------------------------")
  end

  # Loop variables
  converged = false
  sampler.L = 0

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
      # full = createIndexSet(FT,d)::IndexSet{d,Vector{T}}
      # FTSet = getIndexSet(full,1)::Set{Index{d,Vector{N}}}
      # for p in FTSet
      #   j = copy(p+i)::Index{d,Vector{N}}
      #   if isAdmissable(union(old, Set(keys(profit))),j)
      #     profit[j] = 0.
      #   end
      # end
      # sparse addition
      for p in 1:d
        j = copy(i)::Index{d,Vector{N}}
        j[p] += 1
        if isAdmissable(union(old, Set(keys(profit))),j)
          profit[j] = 0.
        end
      end
      newindexset = union(old, Set(keys(profit)))::Set{Index{d,Vector{N}}}
    else
      newindexset = getIndexSet(sampler.settings.indexset,sampler.L)::Set{Index{d,Vector{N}}}
    end

    # get boundary in current iteration
    boundary = Set{Index{d,Vector{N}}}()::Set{Index{d,Vector{N}}} # empty set
    if kind(sampler) == ML
      push!(boundary,Index(sampler.L*ones(N,d)))
    elseif kind(sampler) == AD
      boundary = union(boundary,getBoundary(newindexset)::Set{Index{d,Vector{N}}})
    else
      boundary = union(boundary,getBoundary(getIndexSet(sampler.settings.indexset, sampler.L))::Set{Index{d,Vector{N}}})
    end

    # new indices to add
    indicesToAdd = setdiff( newindexset, oldindexset )

    # print some info to the screen
    !sampler.settings.showInfo || print_with_color(:blue, "*** currently running at L = $(sampler.L)...\n")

    # take initial number of samples
    for index::Index{d,Vector{N}} in indicesToAdd
      if !haskey(sampler.times,index)
        time = @elapsed sample(sampler, sampler.settings.Nstar, index)
        sampler.times[index] = time/16
      end
      E[index] = squeeze(mean(sampler.samples[index],(1,2)),(1,2))
      V[index] = squeeze(var(mean(sampler.samples[index],1),2),(1,2))
      Vf[index] = squeeze(mean(var(sampler.samples[index],2),1),(1,2))
      W[index] =  sampler.settings.useTime ? sampler.times[index] : prod((sampler.settings.m0.*2.^(index.indices)).^sampler.settings.γ)
      Vest[index] = V[index]/sampler.settings.Nstar
      if kind(sampler) == AD && haskey(profit,index) # compute gains when adaptive
        profit[index] = maximum(abs(E[index])./(Vf[index].*W[index]))
      end
    end

    # estimate bias and splitting parameter
    B = zeros(T,sampler.settings.Z)
    for index::Index{d,Vector{N}} in boundary
      B += abs(E[index])
    end
    !sampler.settings.showInfo || print_with_color(:green, @sprintf("   *** B = %0.6e \n", maximum(B) ) )
    splitting = ( maximum(B) < TOL/2 ) ? 1 - maximum(B)/TOL : 0.5
    !sampler.settings.showInfo || print_with_color(:green, @sprintf("   *** splitting = %0.6e \n", epsilon ) )

    # calculate optimal number of samples
    mySum = zeros(T,sampler.settings.Z)
    for index::Index{d,Vector{N}} in newindexset
      mySum += ( ( Vf[index].*W[index].^(2*λ) ).^(1/(2*λ+1) ) )
    end
    for index::Index{d,Vector{N}} in newindexset
      Nopt = ( ( C/(splitting*TOL) )^2 * 1/q * (Vf[index]./W[index]).^(2*λ/(2*λ+1)) .* mySum ).^(1/(2*λ))
      S[index] = ceil(N,maximum(Nopt))
    end

    # take additional samples at each level
    for index::Index{d,Vector{N}} in newindexset
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
    if !(typeof(sampler.settings.numberGenerator) <: MCgenerator && !sampler.settings.safety)
      while maximum(A) > splitting*TOL # 0.9 is safety
        # double number of samples on index with max ratio V/W
        maxratio = 0.
        maxindex = Index(zeros(N,d))::Index{d,Vector{N}}
        for index::Index{d,Vector{N}} in newindexset
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
    end

    # check for convergence
    if ( sampler.L > 1 )
      B = zeros(T,sampler.settings.Z)
      for index::Index{d,Vector{N}} in boundary
        B += abs(E[index])
      end
      Etemp = sum(values(E))
      Vtemp = sum(values(Vf))
      converged = ( maximum(A + B)::T < TOL )

      # print some info to the screen
      !sampler.settings.showInfo || print_with_color(:magenta, @sprintf("*** error estimate is %0.6e + %0.6e = %0.6e \n", maximum(A), maximum(B), maximum(A+B) ) )
      !sampler.settings.showInfo || ( !converged || print_with_color(:magenta, @sprintf("*** result is %0.6e with a variance of %0.6e \n", maximum(Etemp), maximum(Vtemp) ) ) )
    end

    # print warning if no convergence
    if !converged
      if ( kind(sampler) != AD && sampler.L == sampler.settings.maxL )
        warn("maximum level reached and no convergence yet, sorry! :(\n")
          converged = true
      end
    end

    # update L and get new indices to add
    sampler.L += 1
    oldindexset = newindexset
  end

  #
  # PRINT OVERVIEW TABLE
  #
  if ( sampler.settings.showInfo )
    println("--------------------------------------------------------------------------------")
    itype = ndims(sampler.settings.indexset) == 1 ? "level" : "index"
    println("  "*itype*"       E              V               N               W                  ")
    println("--------------------------------------------------------------------------------")
    for index::Index{d,Vector{N}} in sort(Set(collect(keys(E))))
      str = ("  $(index.indices)            "[1:13])::ASCIIString
      str *= @sprintf("%12.5e",maximum(E[index]))
      str *= @sprintf("    %0.6e",maximum(Vf[index]))
      str *= @sprintf("    %d               ",prod(size(sampler.samples[index])[1:2]))[1:16]
      str *= @sprintf("    %0.6e",W[index])
      println(str)
    end
  end

  return (sum(values(E))::Vector{T}, sum(values(Vf))::Vector{T}, A::Vector{T}, B::Vector{T}, splitting::T, S::Dict{Index{d,Vector{N}},N})
  
end
