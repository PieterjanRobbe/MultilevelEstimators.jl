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
        for (k, v) in diff(ℓ)
            dQ .+= v * problem_18(k, ξ)
        end
        dQ, Qf
    end

    # run MLQMC
    for sample_method in (MC(), QMC())
        @suppress_err @capture_out h = run(Estimator(ML(), sample_method, test_function, Uniform(-0.5, 0.5); nb_of_qoi=n, max_index_set_param=3, save_samples=true, folder=tempdir()), 0.001)
    end

    # define test function
    function problem_18(ell::Index, xi::Vector{<:Real})
        ℓ1 = ell[1]
        ξ1 = xi[1]
        ℓ2 = ell[2]
        ξ2 = xi[2]
        x = range(0, 6, n)
        f = ones(n)
        ξ13 = ξ1 * ξ1 * ξ1
        ξ23 = ξ2 * ξ2 * ξ2
        for j in 1:n
            f[j] = x[j] ≤ 3 ? (x[j] - 2)^2 : 2*log(x[j] - 2) + 1
            if ℓ1 == 3
                f[j] += ξ13
            elseif ℓ1 == 2
                f[j] += 1.1 * ξ13
            elseif ℓ1 == 1
                f[j] += (x[j]/60 + 1.2) * ξ13
            else
                f[j] += 3/2 * ξ13
            end
            if ℓ2 == 3
                f[j] += ξ23
            elseif ℓ2 == 2
                f[j] += 1.1 * ξ23
            elseif ℓ2 == 1
                f[j] += (x[j]/60 + 1.2) * ξ23
            else
                f[j] += 3/2 * ξ23
            end
            if ℓ1 > 0 && ℓ2 > 0
                f[j] += 1/32 * ξ13 * ξ23
            end
        end
        f
    end

    for index_method in (TD, FT, HC, ZC, AD)
        for sample_method in (MC(), QMC())
            @suppress_err @capture_out h = run(Estimator(TD(2), sample_method, test_function, [Uniform(-0.5, 0.5), Uniform(-0.5, 0.5)]; nb_of_qoi=n, max_index_set_param=3, save_samples=true, folder=tempdir()), 0.001)
        end
    end

end