 module HomogeneousNormalMatrixGen

using Distributions
#abstract type HomogeneousNormalMatrixGen end

struct HomogeneousNormalMatrix{}

    vx
    vy
    pts
end

#function HomogeneousNormalMatrix(Dist::Distribution,vx::StepRangeLen{Float64},vy::StepRangeLen{Float64})
#    println("test")
#end


function sample(grf::HomogeneousNormalMatrix,ξ::Float64)
    mat=fill(ξ,size(grf.vx)[1],size(grf.vy)[1])

    return mat
end
end
