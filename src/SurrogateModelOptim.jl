module SurrogateModelOptim

export  smoptimize,
        surrogate_model,
        model_infill,
        scaled_LHC_sampling_plan


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

end # module
