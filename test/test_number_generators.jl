## test_number_generators.jl : tests for number_generators.jl

## uniform MC generator ##
verbose && print("testing uniform MC generator...")

umc = UniformMCgenerator(200)
@show umc
println(umc)
umc = UniformMCgenerator(20,-ones(20),ones(20))

i = Index(1,2,3)
@test typeof(i) <: Index
@test i == Index([1,2,3])
@test i == Index("[1,2,3]")
@test Index(0) == Index("[0]")

@test_throws ArgumentError Index("1,2,3")
@test_throws ArgumentError Index("1,b,3")
@test_throws ArgumentError Index("")
@test_throws MethodError Index(0.)

verbose && println("done")

