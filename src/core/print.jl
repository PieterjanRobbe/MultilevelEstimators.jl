# print.jl : collection of print statements for MultilevelEstimators.jl

# print spaces
spaces(n) = repeat(" ",n)

# short formatting
short(num) = @sprintf("%6.3f",num)

# print status of the estimator
function print_status(estimator::Estimator)
    print_index_set(estimator)
    n = 14
    nb = 2
    border = spaces(nb)
    level_name = ndims(estimator.method) == 1 ? "level" : "index"
    header = string(border,level_name,spaces(n-nb-length(level_name)))
    for name in ["E" "V" "N" "W"]
        header = string(header,name,spaces(n))
    end
    println(repeat("-",80))
    println(header)
    println(repeat("-",80))
    for index in keys(estimator)
        index_str = "$(index)"
        str = string(border,index_str,spaces(n-nb-length(index_str)-1))
        str = string(str,@sprintf("%12.5e",mean(estimator,index)),spaces(3))
        str = string(str,@sprintf("%12.5e",var(estimator,index)),spaces(4))
        nshifts = nb_of_shifts(estimator.number_generator)
        nsamples_str = nshifts == 1 ? "" : "$(nshifts) x "
        nsamples_str = string(nsamples_str,"$(estimator.nsamples[index])")
        str = string(str,@sprintf("%s",nsamples_str),spaces(n-length(nsamples_str)))
        str = string(str,@sprintf("%12.5e",cost(estimator,index)),spaces(4))
        println(str)
    end
    println(repeat("-",80))
end

# print convergence of the estimator
function print_convergence(estimator::Estimator,converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈",@sprintf("%12.5e",rmse(estimator)),"."))
    print_footer() 
end

# print warning about max level
function warn_max_level(estimator)
    warn("maximum level L = $(estimator.max_level) reached, no convergence")
end

# print warning about spill index
function warn_spill_index(index)
    warn("a maximum index $(index) was reached, ignoring its active neighbours")
end

print_largest_profit(index) = println("The index with largest profit is $(index), checking forward neighbours...")

# print rates
function print_rates(estimator::MultiLevelTypeEstimator)
	println(string("  ==> Rates: α ≈",short(α(estimator)),
                   ", β ≈",short(β(estimator)),
                   ", γ ≈",short(γ(estimator)),"."))
end

function print_rates(estimator::MultiIndexTypeEstimator)
    template = string("(",["%6.3f" for i in 1:ndims(estimator.method)]...,")")
	println(string("  ==> Rates: α ≈ (", join(short.(α.(estimator,1:ndims(estimator))),","), ")",
                   ", β ≈ (", join(short.(β.(estimator,1:ndims(estimator))),","), ")",
				   ", γ ≈ (", join(short.(γ.(estimator,1:ndims(estimator))),","), ")."))
end

# print MSE analysis
function print_mse_analysis(estimator::Estimator,ϵ::T where {T<:Real},θ::T where {T<:Real})
    print_status(estimator)
    println(string("Checking convergence..."))
    print_rates(estimator)
    println(string("  ==> Variance of the estimator ≈",@sprintf("%12.5e",varest(estimator)),"."))
    println(string("  ==> Bias of the estimator ≈",@sprintf("%12.5e",bias(estimator)),"."))
    θ != 1/2 && println(string("  ==> Non-trivial MSE splitting parameter ≈",@sprintf("%5.2f",θ),"."))
    if rmse(estimator) > ϵ
        println(string("No convergence yet. RMSE ≈",@sprintf("%12.5e",rmse(estimator))," > ",@sprintf("%9.3e",ϵ),"."))
        println(string("Adding an extra level..."))
    end
end

# print QMC convergence
function print_qmc_convergence(estimator::Estimator,ϵ::T where {T<:Real},θ::T where {T<:Real})
    println(string("Checking if variance of estimator is larger than θ*ϵ^2..."))
    converged = !(varest(estimator) > θ*ϵ^2)
    println(string("  ==> ",@sprintf("%11.5e",varest(estimator)),"$(converged ? " <" : " >")",@sprintf("%12.5e",θ*ϵ^2),": $(converged ? "yes!" : "no" )"))
end

# print header
function print_header(estimator::Estimator,ϵ::Float64)
    println(repeat("-",80))
    println("*** MultilevelEstimators.jl @$(now())")
    println("*** This is a $(print_name(estimator))")
    println("*** Simulating $(estimator.sample_function)")
    println(@sprintf("*** Tolerance on RMSE ϵ = %7.3e",ϵ))
    println(repeat("-",80))
end

# print footer
function print_footer()
    println(repeat("-",80))
    println("*** MultilevelEstimators.jl @$(now())")
    println("*** Successfull termination")
    println(repeat("-",80))
end

# print (optimal) number of samples
function print_number_of_samples(estimator,samples)
    println("Samples will be updated according to")
    n = 14
    nb = 2
    border = spaces(nb)
    level_name = ndims(estimator.method) == 1 ? "level" : "index"
    header = string(border,level_name,spaces(n-nb-length(level_name)))
    header = string(header,"N",spaces(n))
    println(repeat("-",29))
    println(header)
    println(repeat("-",29))
    for index in sort(collect(keys(samples)))
        index_str = "$(index)"
        str = string(border,index_str,spaces(n-nb-length(index_str)))
        str = string(str,@sprintf("%s",samples[index]))
        println(str)
    end
    println(repeat("-",29))
end

# print index set
print_index_set(estimator::Estimator) = nothing

print_index_set(estimator::MultiIndexTypeEstimator) = print_index_set(collect(estimator.current_index_set),Index{ndims(estimator)}[])

print_index_set(estimator::AdaptiveMultiIndexTypeEstimator) = print_index_set( collect(setdiff(estimator.old_index_set,estimator.spill_index_set)), collect(union(estimator.spill_index_set,active_set(estimator))) )

function print_index_set(old_set::Vector{Index{d}}, active_set::Vector{Index{d}}=Index{d}()) where {d}
    if d == 2
        unicode_old = "\u25FC"
        unicode_active = "\u25A3"
        n = isempty(old_set) ? 0 : maximum(maximum.(old_set)+1)
        n = isempty(active_set) ? n : max(maximum(maximum.(active_set)+1),n)
        A = zeros(Int64,n,n)
        for index in old_set
            A[index[1]+1,index[2]+1] = 1
        end
        B = zeros(Int64,n,n)
        for index in active_set
            B[index[1]+1,index[2]+1] = 1
        end
        str = "Shape of the index set$(isempty(active_set) ? "" : " ($(unicode_active) = active set)"):\n"
        for j = n:-1:1
            str = string(str,"  ")
            for i = 1:n
                if A[i,j] > 0
                    str = string(str,"$(unicode_old) ")
                end
                if B[i,j] > 0
                    str = string(str,"$(unicode_active) ")
                end
            end
            str = string(str,"\n")
        end
        print(str)
    end
end
