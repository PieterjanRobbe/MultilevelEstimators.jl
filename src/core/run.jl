## run.jl : entry point method for simulating an Estimator
#
# Runs the Estimator for the given tolerance, and returns a History object.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

"""
```julia
run(estimator, ϵ)
```
Run the estimator and compute the expected value of the quantity of interest up to the given tolerance(s) `ϵ`.

# Examples
```jldoctest
julia>

```
"""
function run(estimator::Estimator, tols::AbstractVector{<:Real})

    # input checking
    all(tols.>0) || throw(ArgumentError(string("supplied tolerance(s) must be positive, got ", tols)))

    # make history
    history = History()

    # run the sequence of tolerances
    for tol in tols
        _run(estimator,tol)
        push!(history, estimator, tol) # log the results in history
        clear(estimator) # prepare new run
    end

    return history
end

run(estimator::Estimator, tol::Real) = run(estimator, get_tols(estimator, tol))

## Main routine ##
function _run(estimator::Estimator{T, <:MC}, ϵ::Real) where T<:Union{SL, ML, MI}
    # print status
    verbose(estimator) && print_header(estimator, ϵ)

    # start level is 0
    L = 0

    # MSE splitting parameter
    θ = min_splitting(estimator)

    # main loop
    while L ≤ 2 || !converged(estimator, ϵ, θ)

        # print level
        verbose(estimator) && print_level(estimator, L)

        # update index set
        index_set = boundary(estimator, L)
        verbose(estimator) && print_index_set(estimator, index_set)

        # obtain initial variance estimate
        ns = regress_nb_of_samples(estimator, index_set, ϵ, θ, L) 
        for index in index_set
            if !contains_samples_at_index(estimator, index)
                sample!(estimator, index, ns[index])
            end
        end

        # add new indices to the index set
        for index in index_set
            push!(estimator, index)
        end
        set_sz(estimator, L)

        # print status
        verbose(estimator) && print_status(estimator)

        # value of the MSE splitting parameter
        θ = do_mse_splitting(estimator) ? compute_splitting(estimator, ϵ) : min_splitting(estimator)

        # evaluate optimal number of samples
        n_opt = Dict(τ => optimal_nb_of_samples(estimator, τ, ϵ, θ) for τ in keys(estimator))

        # print optimal number of samples
        verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

        # take additional samples
        for τ in keys(estimator)
            n_due = n_opt[τ] - nb_of_samples(estimator, τ)
            n_due > 0 && sample!(estimator, τ, n_due)
        end

        if T <: SL
            break
        else
            # show status
            verbose(estimator) && print_status(estimator)
            verbose(estimator) && L ≥ 2 && print_mse_analysis(estimator, ϵ, θ)

            # check if the new level exceeds the maximum level
            if !converged(estimator, ϵ, θ) && ( sz(estimator) ≥ max_index_set_param(estimator) ) 
                verbose(estimator) && warn_max_level(estimator)
                break
            end

            # update level
            L += 1
        end
    end

    # print convergence status
    verbose(estimator) && print_convergence(estimator, converged(estimator, ϵ, θ))
end

## Routine for N samples on the coarsest mesh ##
function run2(estimator::Estimator, Ns::AbstractVector{<:Integer})

    # input checking
    all(Ns.>0) || throw(ArgumentError(string("supplied number of samples must be positive, got ", Ns)))

    # make history
    history = History()

    # run the sequence of tolerances
    for N in Ns
        _run(estimator,N)
        push!(history, estimator, N) # log the results in history
        clear(estimator) # prepare new run
    end

    return history
end

# TODO in powers of 2
run2(estimator::Estimator, N::Integer) = run2(estimator, fill(N, 1))

function run2(estimator::Estimator{T, <:MC}, N::Integer) where T<:AbstractIndexSet
    # print status
    verbose(estimator) && print_header(estimator, N)

    # steps
    step = 1000

    # index set is maximum index set allowed
    all_indices = get_index_set(estimator, max_index_set_param(estimator))

    # take samples at boundary (warm-up)
    if isempty(nb_of_samples(estimator))
        for index in all_indices
            add_index(estimator, index)
        end
    end

    # push all levels/indices to the index set
    for index in all_indices
        push!(estimator, index)
    end
    set_sz(estimator, max_index_set_param(estimator))

    # accumulator
    acc = Vector{Float64}(undef, 0)

    # "step" samples at a time
    for m in 1:div(N, step)+1
        n = N > step*m ? step : N % step 

        # print status
        verbose(estimator) && print_status(estimator)
        verbose(estimator) && print_index_set(estimator, all_indices)

        # compute r param
        r = 1/2*(β(estimator) + γ(estimator))
        r = map(rᵢ->isnan(rᵢ) || rᵢ ≤ 0 ? 1.5 : rᵢ, r)
        P = Dict(index=>p(r, index) for index in keys(estimator))
        display(P)
        println("")
        W = Dict(index=>w(r, index) for index in keys(estimator))
        display(W)
        println("")
        
        # compute indices where to take samples
        n_opt = Dict(i=>0 for i in keys(estimator))
        while sum(values(n_opt)) < n
            idx = Index(floor.(Int, randexpr.(r))...)
            if idx ∈ keys(estimator)
                n_opt[idx] += 1
            end
        end

        # subtract
        if T <: MG
            for index in keys(estimator)
                R = CartesianIndices(UnitRange.(zero(eltype(T)), index.-1))
                for I in R
                    n_opt[Index(I)] = max(0, n_opt[Index(I)] -  n_opt[index])
                end
            end
        end

        # print optimal number of samples
        verbose(estimator) && print_optimal_nb_of_samples(estimator, n_opt)

        # take samples and update accumulator
        for τ in keys(estimator)
            # take samples
            n_opt[τ] > 0 && sample!(estimator, τ, n_opt[τ])
            # update accumulator
            append!(acc, zeros(n_opt[τ]))
            ϕ = T <: MG ? zero(eltype(T)) : τ
            R = CartesianIndices(UnitRange.(ϕ, τ))
            for I in R
                acc[end-n_opt[τ]+1:end] += W[Index(I)]*samples_diff(estimator, 1, Index(I))[end-n_opt[τ]+1:end]
            end
        end

        # show status
        verbose(estimator) && print_mse_analysis(estimator, NaN, NaN)

        ve = sum(varest(estimator, index) for index in keys(estimator) if nb_of_samples(estimator, index) > 1)

        # print var est
        println("*****************")
        println("*****************")
        println("*****************")
        println(" >>> variance of the estimator = ", T <: MG ? var(acc) / length(acc) : ve)
        println("*****************")
        println("*****************")
        println("*****************")
    end
end
