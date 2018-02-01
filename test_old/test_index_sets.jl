## test_index_sets.jl : tests for index_sets.jl

## Indexset ##
verbose && print("testing index set creation...")

td = TD(3)
sl = SL()
@test ndims(td) == 3
@test ndims(sl) == 1
@test_throws MethodError SL(2)
@test_throws MethodError ML(2)
@test_throws BoundsError TD(0)
@test_throws BoundsError FT(0)
@test_throws BoundsError AD(0)
@test_throws BoundsError TD(-1)
@test_throws BoundsError FT(-1)
@test_throws BoundsError AD(-1)
@test_throws ArgumentError TD(Float64[])
@test_throws ArgumentError TD([1.,1.,0.])
@test_throws MethodError TD([1,1,0])
@test_throws ArgumentError TD([1.,1.,-1.])

verbose && println("done")

verbose && print("testing index set methods...")

@test length(get_index_set(td,2)) == 10
@test length(get_boundary(get_index_set(td,2))) == 6
indexset = get_index_set(td,3)
@test is_admissable(indexset,Index(4,0,0))
@test_throws MethodError is_admissable(indexset,Index(0,0))
@test !is_admissable(indexset,Index(0,0,0))
@test !is_admissable(indexset,Index(0,0,6))
indexset_sorted = sort(indexset)
@test indexset_sorted[1] == Index(0,0,0)
@test indexset_sorted[end] == Index(3,0,0)

verbose && println("done")

verbose && println("testing print method...")

str = pretty_print(get_index_set(TD(2),2))
verbose && print(str)

verbose && println("done")
