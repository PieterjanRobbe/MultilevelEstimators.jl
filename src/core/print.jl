# print.jl : collection of print statements for MultilevelEstimators.jl

## convenience functions ##
l() = 80
l_short() = 29
n() = 14
nb() = 2
spaces(n) = repeat(" ", n)
border() = spaces(nb())
short(num) = @sprintf("%6.3f", num)
long(num) = @sprintf("%12.5e", num)

## print header ##
function print_header(estimator::Estimator, ϵ::Real)
    println(repeat("-", l()))
    println(string("*** MultilevelEstimators.jl @", now()))
    println(string("*** This is a ", estimator))
    println(string("*** Simulating ", name(estimator)))
    println(string("*** Tolerance on RMSE ϵ = ", short(ϵ)))
    println(repeat("-", l()))
end

## print footer ##
function print_footer()
    println(repeat("-", l()))
    println(string("*** MultilevelEstimators.jl @", now()))
    println("*** Successfull termination")
    println(repeat("-", l()))
end

## print status ##
function print_status(estimator::Estimator)
    print_index_set(estimator)
    header = print_pre_header(estimator)
    for name in ["E" "V" "N" "W"]
        header = string(header, name, spaces(n()))
    end
    print_sub_header(header, l())
    for index in keys(estimator)
        index_str = string(index)
        str = string(border(), index_str, spaces(n()-nb()-length(index_str)-1))
        str = string(str, long(mean(estimator, index)), spaces(3))
        str = string(str, long(var(estimator, index)), spaces(4))
        samples_str = print_nb_of_samples(estimator)
        str = string(str, samples_str, spaces(n-length(samples_str)))
        str = string(str, long(cost(estimator, index)), spaces(4))
        println(str)
    end
    print_sub_footer(l())
end

print_sub_header(header, length) = println(repeat("-", length)), println(header), println(repeat("-", length))
print_sub_footer(length) = println(repeat("-", length))

print_pre_header(estimator) = string(border(), print_elname(estimator), spaces(n()-nb()-length(print_elname(estimator))))

## index or level ##
print_elname(::Estimator{<:AbstractML}) = "level"
print_elname(::Estimator{<:AbstractMI}) = "index"

## print_nb_of_samples ##
print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index, b::Integer) = string(b)
print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:MC}, index::Index) = print_nb_of_samples(estimator, index, nb_of_samples_at_index(estimator, index))
print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index, n::Integer) = string(nb_of_shifts_at_index(estimator), " x ", n)
print_nb_of_samples(estimator::Estimator{<:AbstractIndexSet, <:QMC}, index::Index) = string(estimator, index, nb_of_samples_at_index(estimator, index))

## print_optimal_nb_of_samples ##
function print_optimal_nb_of_samples(estimator::Estimator, samples)
    println("Samples will be updated according to")
    header = print_pre_header(estimator)
    header = string(header, "N", spaces(n))
    print_sub_header(header, l_short()) 
    print_optimal_nb_of_samples(samples)
    print_sub_footer(l_short())
end

function print_optimal_nb_of_samples(samples::Dict)
    for index in sort(collect(keys(samples)))
        index_str = string(index)
        str = string(border,index_str(), spaces(n-()nb()-length(index_str)))
        str = string(str, samples[index])
        println(str)
    end
end

print_optimal_nb_of_samples(samples::Integer) = println(string(border(), "0", spaces(n()-nb()-1), samples))

## print_convergence ##
function print_convergence(estimator::Estimator, converged::Bool)
    print_status(estimator)
    converged && println(string("Convergence reached. RMSE ≈",@sprintf("%12.5e",rmse(estimator)),"."))
    print_footer() 
end

## print sample! header ##
print_sample!_header(estimator::Estimator, index::Index, n::Integer, warm_up::Bool) = println(string("Taking ", print_nb_of_samples(estimator, index, n), print_warm_up(Val(warm_up)), " sample", print_with_s(n), " at ", print_elname(estimator), " ", index, "..."))

print_warm_up(::Type{Val{true}}) = " warm-up "
print_warm_up(::Type{Val{false}}) = " additional "

print_with_s(n::integer) = n > 1 ? "s" : ""















#=
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
function print_rates(estimator::Estimator{<:AbstractML})
	println(string("  ==> Rates: α ≈",short(α(estimator)),
                   ", β ≈",short(β(estimator)),
                   ", γ ≈",short(γ(estimator)),"."))
end

function print_rates(estimator::Estimator{<:AbstractMI})
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
    println(string("  ==> ",@sprintf("%11.5e",varest(estimator)),"$(converged ? " <" : " >")",@sprintf("%12.5e",θ*ϵ^2),": $(converged ? "no!" : "yes" )"))
end


# print (optimal) number of samples
function print_number_of_samples(estimator::Estimator, samples::Dict)
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
        index_str = string(index)
        str = string(border,index_str,spaces(n-nb-length(index_str)))
        str = string(str,@sprintf("%s",samples[index]))
        println(str)
    end
    println(repeat("-",29))
end

# print index set
print_index_set(estimator::Estimator) = nothing

print_index_set(estimator::Estimator{<:AbstractMI}) = print_index_set(collect(estimator.current_index_set),Index{ndims(estimator)}[])

print_index_set(estimator::Estimator{<:AbstractAD}) = print_index_set( collect(setdiff(estimator.old_index_set,estimator.spill_index_set)), collect(union(estimator.spill_index_set,active_set(estimator))) )

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
=#
