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
include("model_infill_utilities.jl")
include("model_infill.jl")
include("surrogate_model_utilities.jl")
include("surrogate_model.jl")
include("smoptimize.jl")

end # module