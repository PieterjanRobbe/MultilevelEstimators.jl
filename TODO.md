# TODO

This is a list of some TODO's for the package overhauling!

## Applications

We want to make some example applications of usage in the `applciations/...` folder.

* SDE example (from the Giles paper)
* SPDE example (using `GaussianRandomFields.jl`)

Two tutorial files for these applications will be added.

## Sampler

Make an argument parser that sets the correct options for the sampler.

Interesting options: `lowmem`, `continuate`

Should we use a `Dict` to store mean and var etc., or should we make a `LevelContainer`, such that we can call `mean(sampler[1,0])`? In our algorithms, we only need to call `mean(sampler)` and `var(sampler)`, so it should shield implementation. So it is easier to store them as `sampler.E`.

## Algorithms

We would like to have a `sample` function that covers all algorithms. This function must contain the general structure of all multilevel algorithms: MLMC, MLQMC, MIMC, MIQMC, MGMLMC etc. The appropriate functions will be called based on the type of the sampler (hurray for multiple dispatch!).

