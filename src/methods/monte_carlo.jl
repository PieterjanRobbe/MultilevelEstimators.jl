## monte_carlo.jl : standard Monte Carlo method
#
# Implementation of the standard Monte Carlo (MC) method.
#
# This file is part of MultilevelEstimators.jl - A Julia toolbox for Multilevel Monte
# Carlo Methods (c) Pieterjan Robbe, 2018

# TODO how to make comparison with QMC method?
#  - implement as qmc, with MersenneTwister and deticated var functions?
#  - dedicated method?

## convergence ##
converged(estimator::Estimator{<:SL}) = true

## rates ##
for f in (:α, :β, :γ)
    eval(
         quote 
             $f(estimator::Estimator{<:SL}) = nothing
         end
        )
end

## bias ##
bias(estimator::Estimator{<:SL}) = 0.0

## splitting ##
splitting(estimator::Estimator{<:SL}) = 1.0
do_mse_splitting(estimator::Estimator{<:SL}) = false

## regression ##
do_regression(estimator::Estimator{<:SL}) = false
