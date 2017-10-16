## sampler.jl : defintion of a Sampler

mutable struct Sampler
    index_set
    number_generator
    sample_function

    maximum_level
    cost_model
    nb_of_qoi
    nb_of_warm_up_samples
    gaussian_random_field_sampler
    nb_of_tol
    store_samples
    user_type
    maximum_index_set

    current_index_set
    # TODO make generator states obsolete by using a new number generator at each index
    samples
    # TODO make dicts as functions: E(sampler,index), V(sampler, index); E(sampler), V(sampler)
    ml_sample_fun
end
