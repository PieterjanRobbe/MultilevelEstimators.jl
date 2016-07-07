module TestModule

	using Interpolations

	using Reexport

    @reexport using MultilevelEstimators

	export parametrizedPDEpointEvaluation, parametrizedPDEEffectiveConductivity

	include("../2d/epde2.jl")

	include("../2d/darcy2.jl")
end