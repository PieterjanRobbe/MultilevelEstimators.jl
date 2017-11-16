using ArgParse

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--debug"
			help = "print debug info"
			action = :store_true
	end

	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	global DEBUG = parsed_args["debug"]
end

macro debug(expression)
	if DEBUG
		return esc(:($expression))
	end
end

main()
