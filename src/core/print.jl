## print.jl : useful print statements
#
# A collection of useful print statements for Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## convenience functions ##
l() = 81
l_short() = 29
n() = 15
spaces(n) = repeat(" ", n)
short(num) = @sprintf("%5.3f", num)
shorte(num) = @sprintf("%5.3e", num)
long(num) = @sprintf("%12.5e", num)
end_table_row(str) = println(string(str, spaces(l()-length(str)-1), "|")) 
hline() = println(string("+", repeat("-", l()-2), "+"))
table_hline(m) = println(string("+", repeat(string(repeat("-", n()), "+"), m)))

## print header ##
function print_header(estimator::Estimator, ϵ::Real)
    hline()
    str = string("| *** MultilevelEstimators.jl @", now())
    end_table_row(str)
    str = string("| *** This is a ", estimator)
    end_table_row(str)
    str = string("| *** Simulating ", shortname(estimator))
    end_table_row(str)
    str = string("| *** Tolerance on RMSE ϵ = ", shorte(ϵ))
    end_table_row(str)
    hline()
end

## print footer ##
function print_footer()
    hline()
    str = string("| *** MultilevelEstimators.jl @", now())
    end_table_row(str)
    str = "| *** Successfull termination"
    end_table_row(str)
    hline()
end

## print status ##
function print_status(estimator::Estimator)
    table_hline(5)
    header = "| "
    for name in [print_elname(estimator) "E" "V" "N" "W"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(5)
    for index in keys(estimator)
        index_str = string(index)
        str = "| "
        str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
        str = string(str, long(mean(estimator, index)), spaces(n()-12-1), " |")
        str = string(str, long(var(estimator, index)), spaces(n()-12-1), " |")
        samples_str = print_nb_of_samples(estimator, index)
        str = string(str, " ", samples_str, spaces(n()-length(samples_str)-2), " |")
        str = string(str, long(cost(estimator, index)), spaces(n()-12-1), " |")
        println(str)
    end
    table_hline(5)
end

## index or level ##
print_elname(::Estimator{<:SL}) = "level"
print_elname(::Estimator{<:AbstractML}) = "level"
print_elname(::Estimator{<:AbstractMI}) = "index"

## print_nb_of_samples ##
print_nb_of_samples(estimator::Estimator, index::Index) = print_nb_of_samples(estimator, nb_of_samples(estimator, index))
print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, n::Integer) = string(n)

## print_optimal_nb_of_samples ##
function print_optimal_nb_of_samples(estimator::Estimator, samples)
    println("Samples will be updated to")
    table_hline(2)
    header = "| "
    for name in [print_elname(estimator) "N"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(2)
    _print_optimal_nb_of_samples(estimator, samples)
    table_hline(2)
end

function _print_optimal_nb_of_samples(estimator::Estimator, samples::Dict)
    for index in sort(collect(keys(samples)))
        index_str = string(index)
        str = "| "
        str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
        samples_str = print_nb_of_samples(estimator, samples[index])
        str = string(str, " ", samples_str, spaces(n()-length(samples_str)-2), " |")
        println(str)
    end
end

## print_convergence ##
function print_convergence(estimator::Estimator, converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈", long(rmse(estimator)), "."))
    print_footer() 
end

## print_sample!_header ##
print_sample!_header(estimator::Estimator, index::Index, n::Integer, warm_up::Bool) = println(string("Taking ", print_nb_of_samples(estimator, n), print_warm_up(Val(warm_up)), " sample", print_with_s(n), " at ", print_elname(estimator), " ", index, "..."))

print_warm_up(::Val{true}) = " warm-up"
print_warm_up(::Val{false}) = " additional"

print_with_s(n::Integer) = n > 1 ? "s" : ""

## print_mse_analysis
function print_mse_analysis(estimator::Estimator, ϵ::Real, θ::Real)
    println("Checking convergence...")
    print_rates(estimator)
    println(string("  ==> Variance of the estimator ≈", long(varest(estimator)), "."))
    println(string("  ==> Bias of the estimator ≈", long(bias(estimator)), "."))
    println(string("  ==> MSE splitting parameter ≈ ", short(θ), "."))
    if !converged(estimator, ϵ, θ) > ϵ
        println(string("No convergence yet. RMSE ≈", long(rmse(estimator)), " > ", shorte(ϵ), "."))
        println("Adding an extra level...")
    end
end

function print_rates(estimator::Estimator{<:AbstractML})
    str = string("  ==> Rates: α ≈ ", short(α(estimator)))
    str = string(str, ", β ≈ ", short(β(estimator)))
    str = string(str, ", γ ≈ ", short(γ(estimator)), ".")
    println(str)
end

## warning when max level is reached ##
warn_max_level(estimator::Estimator{<:AbstractML}) = @warn string("Maximum level L = ", max_index_set_param(estimator), " reached, no convergence yet.")
