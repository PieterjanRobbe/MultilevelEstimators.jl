## check_input.jl : common functions to perform input argument checking
#
# Common functions for input argument checking.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods

## checks ##
for (descr, sym) in zip(("equal to", "larger than", "larger than or equal to", "smaller than", "smaller than or equal to"), (:(==),:>, :≥, :<, :≤))
    @eval begin
        $(Symbol("check_", join(split(descr), "_")))(type_name, x, parameter_name, y) = all($(sym).(x, y)) || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be ", $(descr), " ", y, ".")))
    end
end

# function check_ndims(x, type_name, parameter_name, m, msg)
#     ndims(x) == m || throw(ArgumentError(string("in ", type_name, ", ", msg, ", got ", ndims(x), ", expected ", m)))
# end

check_finite(type_name, x, parameter_name) = isfinite(x) || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be finite.")))

check_ordered(type_name, a, b, parameter_name_1, parameter_name_2) = a < b || throw(ArgumentError(string("in ", type_name, ", ", parameter_name_1, " must be smaller than ", parameter_name_2)))

check_ordered_or_equal(type_name, a, b, parameter_name_1, parameter_name_2) = a ≤ b || throw(ArgumentError(string("in ", type_name, ", ", parameter_name_1, " must be smaller than or equal to ", parameter_name_2)))

check_type(type_name, x, parameter_name, requested_type) = x isa requested_type || throw(ArgumentError(string("in ", type_name, ", ", parameter_name, " must be of type ", requested_type, ", got ", typeof(x))))
