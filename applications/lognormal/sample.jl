## sample.jl : sample methods for lognormal diffusion problem
#
# Sample methods for the 2d lognormal diffusion problem.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

function sample_lognormal(index::Index, x::Vector{<:AbstractFloat}, grf::GaussianRandomField, qoi::AbstractQoi, solver::AbstractSolver, reuse::R) where R<:AbstractReuse

    # wrap the sample code in a try-catch  
    @repeat 3 try

        # sample grf
        # TODO for QMC, reorden inputs!!!!
        Z = my_grf_sample(grf, view(x, 1:randdim(grf)))
        k = exp.(Z)
        sz = size(k).-1

        # direct-discretization function
        f(n, m) = begin
            step = div.(size(k), (n, m))
            range = StepRange.(1, step, size(k))
            view(k, range...)
        end
        g(n, m) = elliptic2d(f(n, m))

        # solve
        xfs, szs = FMG_solve(g, sz, solver, reuse)
        Qf = apply_qoi(xfs, f, szs, index, reuse, qoi)

        # compute difference
        dQ = copy(Qf)
        if R <: NoReuse
            for (key, val) in diff(index)
                szc = div.(sz, max.(1, (index.-key).*2))
                xcs, szcs = FMG_solve(g, szc, solver, reuse)
                Qc = apply_qoi(xcs, f, szcs, index, reuse, qoi)
                dQ += val*Qc
            end
        else
            for i in CartesianIndices(dQ)
                index_ = Index(i-one(i))
                for (key, val) in diff(index_)
                    dQ[i] += val*Qf[(key.+1)...]
                end
            end
        end

        if all(mean.(dQ) .< 100)
            return dQ, Qf
        else
            throw(ErrorException("Something went wrong computing this sample, rethrowing error after 3 tries :("))
        end
    catch e
        @retry if true
            randn!(x)
        end
    end
end

## custom GRF sampling ##
# NOTE: FFT plans cannot deal with pmap (unique C pointer cannot be serialized; afaik)
my_grf_sample(grf::GaussianRandomField, x::AbstractVector) = sample(grf, xi=x)
function my_grf_sample(grf::GaussianRandomField{CirculantEmbedding}, x::AbstractVector)
    v = grf.data[1]

    # compute multiplication with square root of circulant embedding via FFT
    y = v .* reshape(x, size(v))
    w = fft!(complex(y)) # this is slower than using the plan, but works in parallel

    # extract realization of random field
    z = Array{eltype(grf.cov)}(undef, length.(grf.pts))
    @inbounds for i in CartesianIndices(z)
        wi = w[i]
        z[i] = real(wi) + imag(wi)
    end
    z
end

## apply QOI ##
apply_qoi(xfs, f, szs, index, ::NoReuse, qoi) = apply_qoi(reshape(xfs, szs.-1), f(szs...), qoi)

function apply_qoi(xfs, f, szs, index, ::Reuse, qoi)
    R = CartesianIndices(index.+one(index))
    xfs_view = view(xfs, R)
    szs_view = view(szs, R)
    map(i->apply_qoi(reshape(xfs_view[i], szs_view[i].-1), f(szs_view[i]...), qoi), Base.Iterators.reverse(eachindex(xfs_view)))
end

function apply_qoi(x, k, ::Qoi1)
    sz = size(x) .+ 1
    x[div.(sz, 2)...]
end

function apply_qoi(x, k, ::Qoi2)
    sz = size(x) .+ 1
    i_end = div.(sz, 2)
    i_start = div.(i_end, 2)
    16*trapz(trapz(view(x, UnitRange.(i_start, i_end)...), 1), 2)[1]
end

function apply_qoi(x, k, ::Qoi3)
    xp = PaddedView(0, x, size(x).+2, (2,2))
    itp = interpolate(xp, BSpline(Linear()))
    sz = size(x) .+ 1
    itp(div(sz[1], 2), range(1, stop=size(xp,2), length=16))
end

function apply_qoi(x, k, ::Qoi4)
    px = PaddedView(zero(eltype(x)), x, size(x).+2, (2,2))
    n, m = size(px)
    # TODO k * \partial p / \partial x (effective permeability ) in FLOW CELL geometry
    trapz((m-1)*view(px, :, m-1), 1)
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

# function used to analyze performance of Multigrid method
function analyze_lognormal(index::Index, x::Vector{<:AbstractFloat}, grf::GaussianRandomField, qoi::AbstractQoi, solver::AbstractSolver)

    # sample grf
    Z = my_grf_sample(grf, view(x, 1:randdim(grf)))
    k = exp.(Z)
    sz = size(k).-1

    # direct-discretization function
    g(n, m) = begin
        step = div.(size(k), (n, m))
        range = StepRange.(1, step, size(k))
        elliptic2d(view(k, range...))
    end

    # solve
    V_cycle_solve(g, sz, solver)
end
