module SurrogateModelOptim

export  smoptimize,
        rosenbrock2d,
        minimum,
        maximum,
        std,
        median,
        mean,
        model_infill

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

using LatinHypercubeSampling: LHCoptim
using ScatteredInterpolation
using BlackBoxOptim: bboptimize
using Parameters: reconstruct, type2dict, @with_kw, @unpack
using Distances
using Statistics
using LinearAlgebra
using Distributed
using StaticArrays
using NearestNeighbors

#Extended methods
import Base.minimum
import Base.maximum
import Statistics.std
import Statistics.median
import Statistics.mean

include("types.jl")
include("interface.jl")
include("rbf_hypers_opt.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")


#TODO
    #Assert that the infill points are not the same as the samples
    #Handle several infill methods
    #Pareto front type of infill? 
    #Tracing
    #How to define the distance between samples of a pareto front?
    
    #Result type containing
        # Original samples
        # Infill points
        # Smallest evaluated function value per iteration
        # Options?
        # Num iterations
        # ?



end # module
