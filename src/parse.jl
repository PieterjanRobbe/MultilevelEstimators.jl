## parse.jl : functions to parse input arguments of the Estimator
#
# Parse inputs passed to the Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

## parse! ##
# Make a function parse! for the key key_name.
# This macro has three arguments:
#  - an expression with the key name;
#  - an expression with the default value;
#  - an expression with the checks to perform on this input argument.
macro parse!(key_name, default_value, checks_to_perform)
    @eval begin
        function parse!(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, options::Dict{Symbol, T} where T, symbol::Val{$key_name})
            key = eltype(symbol)
            if !haskey(options, key)
                options[key] = $default_value
            else
                val = options[key]
                $checks_to_perform
            end
        end
    end
end

parse!(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, options::Dict{Symbol, T} where T, symbol::Symbol) = parse!(index_set, sample_method, options, Val(symbol))

to_string(key, val) = Estimator, val, string("optional key ", key)

eltype(::Type{<:Val{T}}) where {T} = T

## nb_of_warm_up_samples ##
@parse!(:nb_of_warm_up_samples,
		sample_method isa QMC ? (index_set isa U ? 2 : 1) : 20,
        begin
            check_type(to_string(key, val)..., Signed)
            check_larger_than(to_string(key, val)..., 1)
        end
       )

## nb_of_qoi ##
@parse!(:nb_of_qoi,
        1,
        begin
            check_type(to_string(key, val)..., Signed)
            check_larger_than(to_string(key, val)..., 0)
        end
       )

## continuate ##
@parse!(:continuate,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## nb_of_tols ##
@parse!(:nb_of_tols,
        10,
        begin
            check_type(to_string(key, val)..., Integer)
            check_larger_than(to_string(key, val)..., 0)
        end
       )

## continuation_mul_factor ##
@parse!(:continuation_mul_factor,
        1.3,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
            check_larger_than(to_string(key, val)..., 1)
        end
       )

## min_splitting ##
@parse!(:min_splitting,
        0.5,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
            check_larger_than(to_string(key, val)..., 0)
            check_smaller_than_or_equal_to(to_string(key, val)..., 1)
            haskey(options, :max_splitting) || parse!(index_set, sample_method, options, :max_splitting)
            check_ordered_or_equal(Estimator, val, options[:max_splitting], "optional key min_splitting", "optional key max_splitting")

        end
       )

## max_splitting ##
@parse!(:max_splitting,
        0.99,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
            check_larger_than(to_string(key, val)..., 0)
            check_smaller_than_or_equal_to(to_string(key, val)..., 1)
            haskey(options, :min_splitting) || parse!(index_set, sample_method, options, :min_splitting)
            check_ordered_or_equal(Estimator, options[:min_splitting], val, "optional key min_splitting", "optional key max_splitting")
        end
       )

## folder ##
@parse!(:folder,
        pwd(),
        begin
            check_type(to_string(key, val)..., String)
            haskey(options, :save) || parse!(index_set, sample_method, options, :save)
            if options[:save]
                isdir(val) || throw(ArgumentError(string(val, "is not a directory!")))
                ispath(val) || makepath(val)
            end
        end
       )

## name ##
@parse!(:name,
        get_valid_filename(index_set, sample_method, options),
        begin
            check_type(to_string(key, val)..., String)
            parse!(index_set, sample_method, options, Val(:folder))
            !endswith(val, ".jld2") && occursin(".", val) && throw(ArgumentError("in Estimator, optional key name must not contain a ."))
            val = endswith(val, ".jld2") ? val : string(val, ".jld2")
            options[key] = val
            isfile(joinpath(options[:folder], val)) && @warn string("filename ", val, " exists, will be overwritten!")
        end
       )

function get_valid_filename(index_set, sample_method, options)
    parse!(index_set, sample_method, options, Val(:folder))
    filename = "UntitledEstimator"
    cntr = 0
    if isfile(joinpath(options[:folder], string(filename, ".jld2")))
        cntr += 1
        while isfile(joinpath(options[:folder], string(filename, cntr, ".jld2")))
            cntr = cntr+1
        end
    end
    string(filename, cntr == 0 ? "" : cntr, ".jld2")
end

## save_samples ##
@parse!(:save_samples,
        false,
        check_type(to_string(key, val)..., Bool)
       )

## save ##
@parse!(:save,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## verbose ##
@parse!(:verbose,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## cost_model ##
@parse!(:cost_model,
        EmptyFunction(),
        check_type(to_string(key, val)..., Function)
       )

struct EmptyFunction <: Function end

## robustify_bias_estimate ##
@parse!(:robustify_bias_estimate,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## do_mse_splitting ##
@parse!(:do_mse_splitting,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## max_index_set_param ##
@parse!(:max_index_set_param,
        10 * ndims(index_set),
        begin
            check_type(to_string(key, val)..., Signed)
            check_larger_than(to_string(key, val)..., 0)
            haskey(options, :min_index_set_param) || parse!(index_set, sample_method, options, :min_index_set_param)
            check_ordered_or_equal(Estimator, options[:min_index_set_param], val, "optional key min_index_set_param", "optional key max_index_set_param")
        end
       )

## min_index_set_param ##
@parse!(:min_index_set_param,
        2,
        begin
            check_type(to_string(key, val)..., Signed)
            check_larger_than(to_string(key, val)..., 0)
            haskey(options, :max_index_set_param) || parse!(index_set, sample_method, options, :max_index_set_param)
            check_ordered_or_equal(Estimator, val, options[:max_index_set_param], "optional key min_index_set_param", "optional key max_index_set_param")
        end
       )

## sample_mul_factor ##
@parse!(:sample_mul_factor,
        1.2,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
        end
       )

## nb_of_workers ##
@parse!(:nb_of_workers,
        i -> nworkers(),
        begin
            check_type(to_string(key, val)..., Union{Integer, Function})
            if val isa Integer
                check_larger_than(to_string(key, val)..., 0)
                delete!(options, val)
                options[key] = i -> val
            end
        end
       )

## nb_of_shifts ##
@parse!(:nb_of_shifts,
        i -> 10,
        begin
            check_type(to_string(key, val)..., Union{Integer, Function})
            val isa Function && index_set isa U && throw(ArgumentError("in Estimator, optional key nb_of_shifts cannot be a Function when using unbiased estimation"))
            if val isa Integer
                check_larger_than(to_string(key, val)..., 1)
                delete!(options, val)
                options[key] = i -> val
            end
        end
       )

## point_generator ##
@parse!(:point_generator,
        LatticeRule32(length(options[:distributions])),
        check_type(to_string(key, val)..., AbstractRNG)
       )

## do_regression ##
@parse!(:do_regression,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## max_search_space ##
@parse!(:max_search_space,
        TD(ndims(index_set)),
        begin
            check_type(to_string(key, val)..., AbstractIndexSet)
            check_equal_to(Estimator, ndims(val), "dimensions of max_search_space", ndims(index_set))
        end
       )

## nb_of_uncertainties ##
@parse!(:nb_of_uncertainties,
        i -> length(options[:distributions]),
        check_type(to_string(key, val)..., Function)
       )

## penalization ##
@parse!(:penalization,
        0.5,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
            check_larger_than_or_equal_to(to_string(key, val)..., 0)
            check_smaller_than_or_equal_to(to_string(key, val)..., 1)
        end
       )

## acceptance_rate ##
@parse!(:acceptance_rate,
        1.0,
        begin
            check_type(to_string(key, val)..., Real)
            check_finite(to_string(key, val)...)
            check_larger_than_or_equal_to(to_string(key, val)..., 0)
            check_smaller_than_or_equal_to(to_string(key, val)..., 1)
        end
       )

## qoi_with_max_var
@parse!(:qoi_with_max_var,
        0,
        begin
            check_larger_than(to_string(key, val)..., 0)
            check_smaller_than_or_equal_to(to_string(key, val)..., length(options[:nb_of_qoi]))
        end
       )
