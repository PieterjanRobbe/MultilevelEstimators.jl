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
function print_header(estimator::Estimator, x::T) where T<:Real
    hline()
    str = string("| *** MultilevelEstimators.jl @", now())
    end_table_row(str)
    str = string("| *** This is a ", estimator)
    end_table_row(str)
    str = string("| *** Simulating ", shortname(estimator))
    end_table_row(str)
    str = string("| *** Tolerance on RMSE ϵ = ", shorte(x))
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
print_nb_of_samples(estimator::Estimator, index::Index) = print_nb_of_samples(estimator, length(first(samples(estimator))[index]))
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

## print_weights ##
function print_weights(estimator::Estimator{<:MG})
    println("Weight factors:")
    table_hline(2)
    header = "| "
    for name in [print_elname(estimator) "weight"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(2)
    for index in sort(collect(keys(estimator)))
        index_str = string(index)
        str = "| "
        str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
		weight_str = short(weight(estimator, index))
        str = string(str, " ", weight_str, spaces(n()-length(weight_str)-2), " |")
        println(str)
    end
    table_hline(2)
end

## print_convergence ##
function print_convergence(estimator::Estimator, converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈", long(rmse(estimator)), "."))
    print_footer() 
end

## print_sample!_header ##
print_sample!_header(estimator::Estimator, index::Index, n::Integer, warm_up::Bool) = print(string("Taking ", print_nb_of_samples(estimator, n), print_warm_up(Val(warm_up)), " sample", print_with_s(n), " at ", print_elname(estimator), " ", index, "..."))
print_sample!_footer() = println("done")

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

function print_rates(estimator::Estimator)
    str = string("  ==> Rates: α ≈ ", print_rate(estimator, α))
    str = string(str, ", β ≈ ", print_rate(estimator, β))
    str = string(str, ", γ ≈ ", print_rate(estimator, γ), ".")
    println(str)
end

print_rate(estimator::Estimator{<:AbstractML}, f::Function) = short(f(estimator))
print_rate(estimator::Estimator{<:AbstractMI}, f::Function) = string("(", join(short.(f(estimator)), ", "), ")")

print_rate_r(estimator::Estimator{<:MG}) = println(string("Using exponential rate r ≈ ", short(r(estimator)), "."))

## warning when max level is reached ##
warn_max_level(estimator::Estimator) = @warn string("Maximum ", _warn_max_level_name(estimator), " L = ", max_index_set_param(estimator), " reached, no convergence yet.")
_warn_max_level_name(estimator::Estimator{<:AbstractML}) = "level"
_warn_max_level_name(estimator::Estimator{<:AbstractMI}) = "index set parameter"

## print level ##
print_level(estimator::Estimator{<:Union{SL, AbstractML}}, level::Integer) = println(string("Currently running on level ", level, "."))
print_level(estimator::Estimator{<:AbstractMI}, L::Integer) = println(string("Currently running with L = ", L, "."))

## print_index_set ##
print_index_set(estimator::Estimator, index_set) = nothing
print_index_set(estimator::Estimator{<:AbstractMI}, index_set) = ndims(estimator) == 2 ? _print_index_set(union(keys(estimator), index_set)) : nothing

function _print_index_set(index_set)
    char = "\u25FC"
    n = maximum(maximum.(index_set))
    R = CartesianIndices(UnitRange.((0, 0), (n, n)))
    A = map(i -> Index(i.I...) ∈ index_set, R)
    str = "Shape of the index set:\n"
    for j in n+1:-1:1
        str = string(str, "  ")
        for i in 1:n+1
            str = A[i,j] ? string(str, char, " ") : str
        end
        str = string(str,"\n")
    end
    print(str)
end
