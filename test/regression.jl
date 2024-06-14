## regression.jl : regression tests for MultilevelEstimators.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for
# Multilevel Monte Carlo Methods

@testset "Regression                   " begin

    # number of qoi
    n = 13

    # define test function
    function problem_18(ell::Level, xi::Vector{<:Real})
        ℓ = ell[1]
        ξ = xi[1]
        x = range(0, 6, n)
        f = ones(n)
        ξ3 = ξ * ξ * ξ
        for j in 1:n
            f[j] = x[j] ≤ 3 ? (x[j] - 2)^2 : 2*log(x[j] - 2) + 1
            if ℓ == 3
                f[j] += ξ3
            elseif ℓ == 2
                f[j] += 1.1 * ξ3
            elseif ℓ == 1
                f[j] += (x[j]/60 + 1.2) * ξ3
            else
                f[j] += 3/2 * ξ3
            end
        end
        f
    end

    # interface
    function test_function(ℓ, ξ)
        Qf = problem_18(ℓ, ξ)
        dQ = copy(Qf)
        if ℓ != Level(0)
            dQ .-= problem_18(Level(ℓ[1] - 1), ξ)
        end
        dQ, Qf
    end

    # run MLQMC
    for sample_method in (MC(), QMC())
        @suppress_err @capture_out h = run(Estimator(ML(), sample_method, test_function, Uniform(-0.5, 0.5); nb_of_qoi=n, max_index_set_param=3, save_samples=true, folder=tempdir()), 0.001)
    end

end