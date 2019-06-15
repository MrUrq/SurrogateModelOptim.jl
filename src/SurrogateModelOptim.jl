module SurrogateModelOptim

export  smoptimize,
        model_infill


using LatinHypercubeSampling
import LatinHypercubeSampling.Continuous, LatinHypercubeSampling.Categorical 
using ScatteredInterpolation
using BlackBoxOptim
using Parameters
using Distances
using Statistics
using StatsBase
using LinearAlgebra
using Distributed
using StaticArrays
import NearestNeighbors: KDTree, knn
using Printf

include("types.jl")
include("LHC_sampling_plan.jl")
include("rbf_hypers_opt_utilities.jl")
include("rbf_hypers_opt.jl")
include("sample_infill_utilities.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")
include("search_range.jl")

#TODO

#####Fancy results printing
#####Call and continue optimisation with the results type?
#####Call the method with a test function + options shortcut
#####Documentation
#####Doc strings
#####Examples
#####Function headers
#####Assert that inputs have correct dimensions

end # module
