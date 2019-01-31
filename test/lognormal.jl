## lognormal.jl : testing for LognormalDiffusionProblems.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods (c) Pieterjan Robbe, 2019

# tmp folder where history files will be stored
my_temp_folder = tempdir()

# handling printing of output for tests
short_type_name(T) = split(string(typeof(T)), "{")[1]
name(index_set, sample_method) =  begin
	str = join(["(lognormal)", short_type_name(index_set), short_type_name(sample_method)], " ")
	string(str, repeat(" ", 29 - length(str))) 
end

# handling of additional keys
extras(index_set, sample_method) = begin
	if index_set isa SL
		options = (max_index_set_param = 3,)
	else
		options = ()
	end
	pairs(options)
end

for index_set in [SL(), ML(), TD(2), AD(2), U(1), U(2)]
	for sample_method in [MC(), QMC()]
		@testset "$(name(index_set, sample_method))" begin
			#
			# single qoi
			#
			estimator = init_lognormal(index_set, sample_method,
									   grf_generator = KarhunenLoeve(100),
									   covariance_function = Matern(1, 2),
									   folder = my_temp_folder,
									   verbose = false;
									   extras(index_set, sample_method)...
									   )
			run(estimator, 1e-2)
			#
			# multiple qoi
			#
			estimator = init_lognormal(index_set, sample_method,
									   qoi = Qoi3(),
									   nb_of_qoi = 16,
									   nb_of_coarse_dofs = 8,
									   cost_model = (level) -> 16*2^level[1],
									   folder = my_temp_folder,
									   verbose = false;
									   extras(index_set, sample_method)...
									   )
			run(estimator, 1e-1)
		end
	end
end
