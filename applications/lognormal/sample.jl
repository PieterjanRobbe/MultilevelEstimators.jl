## sample.jl : sample methods for lognormal diffusion problem
#
# Sample emthods for the 2d lognormal diffusion problem.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

function sample_lognormal(index::Index, x::Vector{<:AbstractFloat}, grf::GaussianRandomField, damping::Real, qoi::AbstractQoi)

    # sample grf
    # TODO for QMC, reorden inputs!!!!
    Z = sample(grf, xi=view(x, 1:randdim(grf)))
    k = exp.(Z)
    sz = size(k).-1

    # direct-discretization function
    g(n, m) = begin
        step = div.(size(k), (n, m))
        range = StepRange.(1, step, size(k))
        elliptic2d(view(k, range...))
    end

    # solve
    xf = FMG_solve(g, sz, damping)
    Qf = apply_qoi(reshape(xf, sz.-1), qoi)

    # compute difference
    dQ = Qf
    for (key, val) in diff(index)
        szc = div.(sz, max.(1, (index.-key).*2))
        xc = FMG_solve(g, szc, damping)
        Qc = apply_qoi(reshape(xc, szc.-1), qoi)
        dQ += val*Qc
    end

    (dQ, Qf)
end

function apply_qoi(x, ::Qoi1)
    sz = size(x) .+ 1
    x[div.(sz, 2)...]
end

function apply_qoi(x, ::Qoi2)
    sz = size(x) .+ 1
    i_end = div.(sz, 2)
    i_start = div.(i_end, 2)
    16*trapz(trapz(view(x, UnitRange.(i_start, i_end)...), 1), 2)[1]
end

function apply_qoi(x, ::Qoi3)
    xp = PaddedView(0, x, size(x).+2, (2,2))
    itp = interpolate(xp, BSpline(Linear()))
    sz = size(x) .+ 1
    itp(div(sz[1], 2), range(1, stop=size(xp,2), length=16))
end

function trapz(A, dim)
    sz = size(A)
    Rpre = CartesianIndices(sz[1:dim-1])
    Rpost = CartesianIndices(sz[dim+1:end])
    szs = [sz...]
    n = szs[dim]
    szs[dim] = 1
    B = Array{eltype(A)}(undef, szs...)
    trapz!(B, A, Rpre, Rpost, n)
end

@noinline function trapz!(B, A, Rpre, Rpost, n)
    fill!(B, zero(eltype(B)))
    for Ipost in Rpost
        for Ipre in Rpre
            for i = 2:n
                B[Ipre, 1, Ipost] += A[Ipre, i, Ipost] + A[Ipre, i-1, Ipost]
            end
        end
    end
    B./(2(n-1))
end
