## index_set.jl : unit testing for index_set.jl
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

@testset "IndexSet                     " begin

for (T, n, i) in Base.Iterators.zip(["SL" "ML"], [1 6], [5, 0])
    @eval begin
        t = $(Symbol(T))()
        @test isa(t, $(Symbol(T)))
        @test isa(t, $(Symbol(T)){1})
        idxset = collect(get_index_set(t, 5))
        @test length(idxset) == $n
        @test idxset[1] == Level($i)
        @test idxset[end] == Level(5)
        @test_throws MethodError $(Symbol(T))(2)
        @sprintf "%s" t
    end
end

for (T, n3, n2) in Base.Iterators.zip(["FT" "TD" "HC" "ZC"], [512 120 38 98], [28 17 11 15])
    @eval begin
        #
        # 3d index set
        #
        t3 = $(Symbol(T))(3)
        @test isa(t3, $(Symbol(T)){3})
        idxset3 = collect(get_index_set(t3, 7))
        @test length(idxset3) == $n3
        @test idxset3[1] == Index(0, 0, 0)
        @test Index(1, 1, 1) ∈ idxset3
        @test Index(8, 0, 0) ∉ idxset3
        @test Index(0, 8, 0) ∉ idxset3
        @test Index(0, 0, 8) ∉ idxset3
        #
        # weighted d2 index set
        #
        t2 = $(Symbol(T))(1/1, 3/5)
        @test isa(t2, $(Symbol(T)){2})
        idxset2 = collect(get_index_set(t2, 6))
        @test length(idxset2) == $n2
        @test idxset2[1] == Index(0, 0)
        @test Index(1, 1) ∈ idxset2
        @test Index(0, 4) ∉ idxset2
        @test Index(7, 0) ∉ idxset2
        @test_throws ArgumentError $(Symbol(T))(1)
        @test_throws ArgumentError $(Symbol(T))(1, 0)
        @sprintf "%s" t2
    end
end

ad = AD(2)
@test ad isa AD
@test ad isa AD{2}
@test_throws ArgumentError get_index_set(ad, 4)
@test_throws ArgumentError AD(1)
@test_throws MethodError AD()
@sprintf "%s" ad

u = U(1)
@test u isa U
@test u isa U{1}
@test_throws ArgumentError get_index_set(u, 1)
@test_throws MethodError U()
@test_throws MethodError U(1, 2)
@sprintf "%s" u

m = MG(U(2))
@test m isa MG
@test_throws ArgumentError get_index_set(m, 4)
@sprintf "%s" m

@test length(collect(get_index_set(TD(2), 0))) == 1
@test_throws ArgumentError get_index_set(TD(2), -1)

end
