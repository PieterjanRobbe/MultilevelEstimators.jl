## print.jl : useful print statements
#
# A collection of useful print statements for Estimators.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

#
# convenience functions
#
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

#
# header and footer
#
function print_header(estimator::Estimator, ϵ::Real)
    hline()
    str = string("| *** MultilevelEstimators.jl @", now())
    end_table_row(str)
    str = string("| *** This is a ", estimator)
    end_table_row(str)
    str = string("| *** Simulating ", first(split(estimator[:name], ".")))
    end_table_row(str)
    str = string("| *** Tolerance on RMSE ϵ = ", shorte(ϵ))
    end_table_row(str)
    hline()
end

function print_footer()
    hline()
    str = string("| *** MultilevelEstimators.jl @", now())
    end_table_row(str)
    str = "| *** Successfull termination."
    end_table_row(str)
    hline()
end

#
# status
#
function print_status(estimator::Estimator)
    table_hline(5)
    header = "| "
    for name in [print_elname(estimator) "E" "V" "N" "W"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(5)
    for index in keys(estimator)
		if nb_of_samples(estimator, index) > 0
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
    end
    table_hline(5)
end

#
# index or level
#
print_elname(::Estimator{<:SL}) = "level"

print_elname(::Estimator{<:AbstractML}) = "level"

print_elname(::Estimator{<:AbstractMI}) = "index"

#
# optimal number of samples
#
function print_optimal_nb_of_samples(estimator::Estimator, samples)
    println("Samples will be updated ", estimator isa U ? "with" : "to")
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
		if samples[index] > 0
			index_str = string(index)
			str = "| "
			str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
			samples_str = print_nb_of_samples(estimator, index, samples[index])
			str = string(str, " ", samples_str, spaces(n()-length(samples_str)-2), " |")
			println(str)
		end
    end
end

print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, n::Integer) = string(n)

print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index, n::Integer) = string(estimator[:nb_of_shifts](index), " x ", n)

print_nb_of_samples(estimator::Estimator, index::Index) = print_nb_of_samples(estimator, index, nb_of_samples(estimator, index))

#
# profits
#
function print_largest_profit(estimator, max_index, max_profit, indices, profits)
    println("Profit indicators:")
    table_hline(2)
    header = "| "
    for name in [print_elname(estimator) "profit"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(2)
    for i in 1:length(indices)
        index_str = string(indices[i])
        str = "| "
        str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
        profit_str = shorte(profits[i])
        str = string(str, " ", profit_str, spaces(n()-length(profit_str)-2), " |")
        println(str)
    end
    table_hline(2)
    println("Index with max profit is ", max_index, " (value = ", long(max_profit), ").")
end

#
# probabilities
#
function print_pmf(estimator::Estimator)
	P = pmf(estimator)
    println("Using approximate probability mass function:")
    table_hline(2)
    header = "| "
	for name in [print_elname(estimator) "P"]
        header = string(header, name, spaces(n()-length(name)-1), "| ")
    end
    println(header)
    table_hline(2)
	for idx in sort(collect(keys(P)))
        index_str = string(idx)
        str = "| "
        str = string(str, index_str, spaces(n()-length(index_str)-2), " |")
        pmf_str = shorte(P[idx])
        str = string(str, " ", pmf_str, spaces(n()-length(pmf_str)-2), " |")
        println(str)
    end
    table_hline(2)
end

#
# convergence
#
function print_convergence(estimator::Estimator, converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈", long(rmse(estimator)), "."))
    print_footer() 
end

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

function print_qmc_convergence(estimator::Estimator{<:AbstractIndexSet, <:QMC}, ϵ::Real, θ::Real)
    println(string("Checking if variance of estimator is larger than θ*ϵ^2..."))
    converged = !(varest(estimator) > θ*ϵ^2)
    println(string("  ==> ", converged ? "no!" : "yes", " (", long(varest(estimator))[2:end], converged ? " < " : " > ", short(θ), " *", long(ϵ^2), ")"))
end

function print_unbiased_convergence(estimator::Estimator{<:U}, ϵ::Real)
    println(string("Checking if variance of estimator is larger than ϵ^2..."))
    converged = !(varest(estimator) > ϵ^2)
    println(string("  ==> ", converged ? "no!" : "yes", " (", long(varest(estimator))[2:end], converged ? " < " : " > ", long(ϵ^2)[2:end], ")"))
    print_rates(estimator)
end

#
# taking ...
#
print_sample!_header(estimator::Estimator, index::Index, n::Integer, warm_up::Bool) = print(string("Taking ", print_nb_of_samples(estimator, index, n), print_warm_up(Val(warm_up)), " sample", print_with_s(n), " at ", print_elname(estimator), " ", index, "..."))

print_warm_up(::Val{true}) = " warm-up"

print_warm_up(::Val{false}) = " additional"

print_with_s(n::Integer) = n > 1 ? "s" : ""

print_sample!_footer() = println(" done")

#
# rates
# 
function print_rates(estimator::Estimator)
    str = string("  ==> Rates: α ≈ ", print_rate(estimator, α))
    str = string(str, ", β ≈ ", print_rate(estimator, β))
    str = string(str, ", γ ≈ ", print_rate(estimator, γ), ".")
    println(str)
end

print_rate(estimator::Estimator{<:AbstractML}, f::Function) = short(f(estimator)[1])

print_rate(estimator::Estimator{<:AbstractMI}, f::Function) = string("(", join(short.(f(estimator)), ", "), ")")

#
# max level/index
#
warn_max_level(estimator::Estimator) = @warn string("Maximum ", _warn_max_level_name(estimator), " L = ", estimator[:max_index_set_param], " reached, no convergence yet.")

_warn_max_level_name(estimator::Estimator{<:AbstractML}) = "level"

_warn_max_level_name(estimator::Estimator{<:AbstractMI}) = "index set parameter"

warn_max_index(estimator::Estimator{<:AD}, index) = @warn string("Trying to add index ", index, " outside search space to index set, ignoring...")

#
# level parameter
#
print_level(estimator::Estimator{<:SL}, L) = println(string("Currently running on finest level."))

print_level(estimator::Estimator{<:ML}, L) = println(string("Currently running on level ", L, "."))

print_level(estimator::Estimator, L) = println(string("Currently running with L = ", L, "."))

#
# index set
#
print_index_set(estimator::Estimator, index_set) = nothing

print_index_set(estimator::Estimator{<:AbstractMI}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set:")
        print(union(keys(estimator), index_set))
    end
end

print_index_set(estimator::Estimator{<:U}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set:")
		print(filter(i -> has_samples_at_index(estimator, i), keys(estimator)))
    end
end

print_index_set(estimator::Estimator{<:AD}, index_set) = begin
    if ndims(estimator) == 2
        println("Shape of the index set (\u25A3 = active set):")
        print(collect(old_set(estimator)), collect(active_set(estimator)))
    end
end
