## parsers.jl : functions to parse input arguments of the Estimator
#
# Parse inputs passed to the Estimator type.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## parse! ##
# make a function parse! for the key key_name
# This macro has three arguments:
#  - an expression with the key name;
#  - an expression with the default value;
#  - an expression with the checks to perform on this input argument.
macro parse!(key_name, default_value, checks_to_perform)
    eval(
         quote
             function parse!(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, settings::Dict{Symbol, T} where T, symbol::Val{$key_name})
                 key = eltype(symbol)
                 if !haskey(settings, key)
                     settings[key] = $default_value
                 else
                     val = settings[key]
                     $checks_to_perform
                 end
             end
         end)
end

parse!(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, settings::Dict{Symbol, T} where T, symbol::Symbol) = parse!(index_set, sample_method, settings, Val(symbol))

## nb_of_warm_up_samples ##
@parse!(:nb_of_warm_up_samples,
        sample_method isa QMC ? 1 : 20,
        (check_type(to_string(key, val)..., Signed);
         check_larger_than(to_string(key, val)..., 1))
       )


## nb_of_qoi ##
@parse!(:nb_of_qoi,
        1,
        (check_type(to_string(key, val)..., Signed);
         check_larger_than(to_string(key, val)..., 0))
       )

## continuate ##
# TODO must be true for multigrid estimator?
@parse!(:continuate,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## nb_of_tols ##
@parse!(:nb_of_tols,
        10,
        (check_type(to_string(key, val)..., Integer);
         check_larger_than(to_string(key, val)..., 0))
       )

## continuation_mul_factor ##
@parse!(:continuation_mul_factor,
        1.2,
        (check_type(to_string(key, val)..., Real);
         check_finite(to_string(key, val)...);
         check_larger_than(to_string(key, val)..., 1))
       )

## splitting ##
@parse!(:splitting,
        1/2,
        (check_type(to_string(key, val)..., Real);
         check_finite(to_string(key, val)...);
         check_larger_than(to_string(key, val)..., 0);
         check_smaller_than_or_equal_to(to_string(key, val)..., 1))
)

## folder ##
@parse!(:folder,
        pwd(),
        (check_type(to_string(key, val)..., String);
         isdir(val) || throw(ArgumentError(string(val, "is not a directory!")));
         ispath(val) || makepath(val))
       )

## name ##
@parse!(:name,
        get_valid_filename(index_set, sample_method, settings),
        (check_type(to_string(key, val)..., String);
         parse!(index_set, sample_method, settings, Val(:folder));
         val = endswith(val, ".jld") ? val : string(val, ".jld");
         settings[key] = val;
         isfile(joinpath(settings[:folder], val)) && @warn string("filename ", val, " exists, will be overwritten!"))
       )

function get_valid_filename(index_set::AbstractIndexSet, sample_method::AbstractSampleMethod, settings::Dict{Symbol, T} where T)
    parse!(index_set, sample_method, settings, Val(:folder))
    filename = "UntitledEstimator";
    cntr = 0;
    if isfile(joinpath(settings[:folder], string(filename, ".jld")))
        cntr += 1
        while isfile(joinpath(settings[:folder], string(filename, cntr, ".jld")))
            cntr = cntr+1;
        end
    end
    string(filename, cntr == 0 ? "" : cntr, ".jld")
end

## save_samples ##
@parse!(:save_samples,
        false,
        check_type(to_string(key, val)..., Bool)
       )

## verbose ##
@parse!(:verbose,
        true,
        check_type(to_string(key, val)..., Bool)
       )

## cost_model ##
struct EmptyFunction <: Function end
@parse!(:cost_model,
        EmptyFunction(),
        check_type(to_string(key, val)..., Function)
       )

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
        100,
        (check_type(to_string(key, val)..., Signed);
         check_larger_than(to_string(key, val)..., 0))
       )

## sample_mul_factor ##
@parse!(:sample_mul_factor,
        1.2,
        (check_type(to_string(key, val)..., Real);
         check_finite(to_string(key, val)...))
       )

## nb_of_workers ##
@parse!(:nb_of_workers,
        i -> nworkers(),
        (check_type(to_string(key, val)..., Union{Integer, Function});
         if val isa Integer;
             check_larger_than(to_string(key, val)..., 0);
             delete!(settings, val);
             settings[key] = i -> val;
         end)
       )

## nb_of_shifts ##
@parse!(:nb_of_shifts,
        i -> 10,
        (check_type(to_string(key, val)..., Union{Integer, Function});
         if val isa Integer;
             check_larger_than(to_string(key, val)..., 0);
             delete!(settings, val);
             settings[key] = i -> val;
         end)
       )

## point_generator ##
@parse!(:point_generator,
        LatticeRule32(length(settings[:distributions])),
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
        (check_type(to_string(key, val)..., AbstractIndexSet);
         check_ndims(to_string(key, val)..., ndims(index_set), "dimensions of search space and index set do not agree"))
       )

to_string(key, val) = val, "Estimator", string("optional key ", key)
eltype(::Type{<:Val{T}}) where {T} = T
