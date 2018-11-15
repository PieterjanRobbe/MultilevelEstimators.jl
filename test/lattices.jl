## lattices.jl : unit testing for lattices.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

@testset "LatticeRule32                " begin

lattice_rule = LatticeRule32([0x00000001,0x00000022], 2, 0x00000037)
point = MultilevelEstimators.get_point(lattice_rule, 2)
@test lattice_rule isa LatticeRule32
@test length(point) == 2
@test all(0 .≤ point .< 1)

z_file_1 = joinpath(dirname(pathof(MultilevelEstimators)), "core", "generating_vectors", "K_3600_32.txt")
z_file_2 = joinpath(dirname(pathof(MultilevelEstimators)), "core", "generating_vectors", "CKN_250_20.txt")
lattice_rule = LatticeRule32(z_file_1, 101)
point = MultilevelEstimators.get_point(lattice_rule, 100_000)
@test lattice_rule isa LatticeRule32
@test length(point) == 101
@test all(0 .≤ point .< 1)

function approx_π(n::Integer)
    lattice_rule = LatticeRule32(z_file_2, 2)
    count = 0
    for i in 0:n-1
        x = MultilevelEstimators.get_point(lattice_rule, i)
        if x[1]*x[1] + x[2]*x[2] < 1
            count += 1
        end
    end
    return 4*count/n
end

for (n, ϵ) in zip(2 .^(0:20), exp10.([0 0 0 0 0 0 1 1 2 2 2 3 3 3 3 3 3 4 4 5]))
    @test approx_π(n) ≈ π atol=ϵ
end

@test LatticeRule32(3600) isa LatticeRule32
@test_throws ArgumentError LatticeRule32(3601)
@test LatticeRule32(z_file_2, 250) isa LatticeRule32
@test_throws ArgumentError LatticeRule32(z_file_2, 251)
@test_throws ArgumentError LatticeRule32(-2)
@test LatticeRule32(z_file_2, 1) isa LatticeRule32
@test_throws ArgumentError LatticeRule32(0)
@test_throws ArgumentError LatticeRule32(z_file_1, 24, 4294967297)

end
