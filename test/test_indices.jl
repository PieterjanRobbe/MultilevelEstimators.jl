## test_indices.jl : tests for indices.jl

verbose && print("testing index creation...")

i = Index(1,2,3)
@test typeof(i) <: Index
@test i == Index([1,2,3])
@test i == Index("[1,2,3]")
@test Index(0) == Index("[0]")

@test_throws ArgumentError Index("1,2,3")
@test_throws ArgumentError Index("1,b,3")
@test_throws ArgumentError Index("")
@test_throws MethodError Index(0.)
# note: not testing for "negative" indices, since they are needed in the drop algorithm

verbose && println("done")

verbose && print("testing index operators...")

@test zero(i) == Index(0,0,0)
@test one(i) == Index(1,1,1)
@test Index(2,4,6) == 2*i
@test Index(2,4,6) == i*2
@test Index(2,4,6) == Index(2,2,2)*i
@test Index(0,0,0) < i
@test i > Index(0,2,3)
@test i > Index(1,2,2)
@test i > Index(0,2,2)
@test Index(2,4,6) == i+i
@test Index(2,4,6) == Index(1,3,5)+1
@test Index(2,3,4)-1 == i
@test i[1] == 1 && i[2] == 2 && i[3] == 3
@test maximum(i) == 3
@test indmax(i) == 3
@test sum(i) == 6
@test prod(i) == 6
i[1] = 2
@test i[1] == 2 && i[2] == 2 && i[3] == 3
@test diff(i,Index(1,2,3)) == 1
@test i == copy(i)

A = eye(2)
@test A[Index(0,0)] == 1
@test A[Index(0,1)] == 0
@test A[Index(1,1)] == 1
@test A[Index(1,0)] == 0

verbose && println("done")
