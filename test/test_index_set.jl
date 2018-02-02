# test_index_set.jl : test IndexSet

@testset "IndexSet                     " begin

sl = SL()
@test isa(sl,SL)
@test isa(sl,SL{1})
idxset = get_indexset(sl,5)
@test length(idxset) == 1
@test idxset[1] == (5,)
@test_throws MethodError SL(2)
str = @sprintf("%s",sl)

ml = ML()
@test isa(ml,ML)
@test isa(ml,ML{1})
idxset = get_indexset(ml,4)
@test length(idxset) == 5
@test sort(idxset)[3] == (2,)
@test_throws MethodError ML(2)
str = @sprintf("%s",ml)

ft = FT(3)
@test isa(ft,FT{3})
idxset = get_indexset(ft,4)
@test length(idxset) == 125
@test sort(idxset)[1] == (0,0,0)
@test sort(idxset)[end] == (4,4,4)
@test in((1,2,3),idxset)
@test !in((5,0,0),idxset)
ft2 = FT(2,δ=[2,1])
@test isa(ft2,FT{2})
idxset2 = get_indexset(ft2,3)
@test length(idxset2) == 8
@test sort(idxset2)[end] == (1,3)
@test_throws BoundsError FT(1)
@test_throws ArgumentError FT(2,δ=[1,1,1])
@test_throws ArgumentError FT(2,δ=[1,0])
@test_throws MethodError FT()
str = @sprintf("%s",ft)

td = TD(3)
@test isa(td,TD{3})
idxset = get_indexset(td,4)
@test length(idxset) == 35
@test sort(idxset)[1] == (0,0,0)
@test sort(idxset)[end] == (4,0,0)
@test in((1,2,0),idxset)
@test !in((4,0,1),idxset)
td2 = TD(2,δ=[1,2])
@test isa(td2,TD{2})
idxset2 = get_indexset(td2,4)
@test length(idxset2) == 9
@test sort(idxset2)[end] == (4,0)
@test_throws BoundsError TD(1)
@test_throws ArgumentError TD(2,δ=[1,1,1])
@test_throws ArgumentError TD(2,δ=[1,0])
@test_throws MethodError TD()
str = @sprintf("%s",td)

hc = HC(3)
@test isa(hc,HC{3})
idxset = get_indexset(hc,4)
@test length(idxset) == 13
@test sort(idxset)[1] == (0,0,0)
@test sort(idxset)[end] == (3,0,0)
@test in((0,1,1),idxset)
@test !in((4,0,0),idxset)
hc2 = HC(2,δ=[1,2])
@test isa(hc2,HC{2})
idxset2 = get_indexset(hc2,12)
@test length(idxset2) == 21
@test sort(idxset2)[end] == (11,0)
@test_throws BoundsError HC(1)
@test_throws ArgumentError HC(2,δ=[1,1,1])
@test_throws ArgumentError HC(2,δ=[1,0])
@test_throws MethodError HC()
str = @sprintf("%s",hc)

ad = AD(4)
@test isa(ad,AD)
@test isa(ad,AD{4})
@test_throws MethodError get_indexset(ad,4)
@test_throws BoundsError AD(1)
@test_throws MethodError AD()
str = @sprintf("%s",ad)

end
