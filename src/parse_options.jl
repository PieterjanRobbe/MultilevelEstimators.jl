## parse_args.jl : parse extra options passed with julia

# parse extra input arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--verbose"
        help = "print debug info"
        action = :store_true
        "-v"
        help = "print debug info"
        action = :store_true
    end

    return parse_args(s)
end

# parse command line and set DEBUG variable
function main()
    parsed_args = parse_commandline()
    global DEBUG = get(parsed_args,"debug",false)
end

# definition of @debug macro
macro debug(expression)
    if DEBUG
        return esc(:($expression))
    end
end

# check if we are running in debug mode and toggle @debug macro
main()
