# print.jl : collection of print statements for MultilevelEstimators.jl

# print spaces
spaces(n) = repeat(" ",n)

# print status of the estimator
function print_status(estimator::Estimator)
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
        nsamples_str = "$(estimator.nsamples[index])"
        str = string(str,@sprintf("%s",nsamples_str),spaces(n-length(nsamples_str)))
        str = string(str,@sprintf("%12.5e",cost(estimator,index)),spaces(4))
        println(str)
    end
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

# print MSE analysis
function print_mse_analysis(estimator::Estimator,ϵ::T where {T<:Real},θ::T where {T<:Real})
    print_status(estimator)
    println(string("Checking convergence..."))
    println(string("  ==> Rates: α ≈",@sprintf("%6.3f",α(estimator)),
                              ", β ≈",@sprintf("%6.3f",β(estimator)),
                              ", γ ≈",@sprintf("%6.3f",γ(estimator)),"."))
    println(string("  ==> Variance of the estimator ≈",@sprintf("%12.5e",varest(estimator)),"."))
    println(string("  ==> Bias of the estimator ≈",@sprintf("%12.5e",bias(estimator)),"."))
    θ != 1/2 && println(string("  ==> Non-trivial MSE splitting parameter ≈",@sprintf("%5.2f",θ),"."))
    if rmse(estimator) > ϵ
        println(string("No convergence yet. RMSE ≈",@sprintf("%12.5e",rmse(estimator))," > ",@sprintf("%9.3e",ϵ),"."))
        println(string("Adding an extra level..."))
    end
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
end