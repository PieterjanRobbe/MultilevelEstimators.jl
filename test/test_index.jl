# test_index.jl : test Index and Level

level_0 = Level(0)
level_2 = Level(2)
@test typeof(level_0) <: Level
@test typeof(level_2) <: Level
@test level_0 .+ level_2 == level_2
@test_throws ArgumentError Level(-1)
@test_throws MethodError Level(0.)

index_0_0 = Index(0,0)
index_2_3 = Index(2,3)
@test typeof(index_0_0) <: Index{2}
@test typeof(index_2_3) <: Index{2}
@test index_0_0 .+ index_2_3 == index_2_3
@test_throws ArgumentError Index(-1,0)
@test_throws MethodError Index(0.,2)
