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
include("rbf_hypers_opt.jl")
include("rbf_hypers_opt_utilities.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")
include("search_range.jl")

#TODO

#####Fancy results printing
#####Call the method with the results type already?
#####Call the method with a test function + options shortcut
#####PlotlyJS parallel coordinates
#####Plot accuracy as a function of "time" predicted and actual
#####Integer test case
#####Documentation
#####Examples
#####Function headers

end # module
