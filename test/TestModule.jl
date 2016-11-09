module TestModule

	using Interpolations

	using Reexport

    @reexport using MultilevelEstimators

	export parametrizedPDEPointEvaluation, parametrizedPDEEffectiveConductivity

	include("../2d/epde2.jl")

	include("../2d/darcy2.jl")
end