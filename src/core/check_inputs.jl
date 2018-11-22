## check_inputs.jl : common functions to perform input argument checking
#
# Common functions for input argument checking.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

## checks ##
function check_finite(x, type_name, parameter_name)
    isfinite(x) || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be finite")))
end

function check_ordered(a, b, type_name, parameter_name_1, parameter_name_2)
    a < b || throw(ArgumentError(string("in ", type_name, ", ", parameter_name_1, " must be smaller than ", parameter_name_2)))
end

for (descr, sym) in zip(("larger than", "larger than or equal to", "smaller than", "smaller than or equal to"), (:>, :≥, :<, :≤))
    eval(
         quote
             $(Symbol("check_", join(split(descr), "_")))(x, type_name, parameter_name, y) = $sym(x, y) || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be ", $descr, " ", y)))
         end)
end

check_positive(x, type_name, parameter_name) = check_larger_than(x, type_name, parameter_name, 0)

function check_type(x, type_name, parameter_name, requested_type)
    x isa requested_type || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be of type ", requested_type, ", got ", typeof(x))))
end

function check_ndims(x, type_name, parameter_name, m, msg)
    ndims(x) == m || throw(ArgumentError(string("in ", type_name, ", ", msg, ", got ", ndims(x), ", expected ", m)))
end
