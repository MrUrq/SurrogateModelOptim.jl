module SurrogateModelOptim

export  smoptimize,
        rosenbrock2d,
        rotatedHyperElipsoid,
        styblinskiTang,
        minimum,
        maximum,
        std,
        median,
        mean,
        model_infill

function rosenbrock2d(x)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
function rotatedHyperElipsoid(x)
    out = 0.0
    xy = [x[1],x[2]]
    for i = 1:2
        inner = 0.0
        for j = 1:i
            inner = inner + xy[j]^2
        end
        out += inner
    end
    return out
end
function styblinskiTang(x)
    out = 0.0    
    xy = [x[1],x[2]]
    for i = 1:2
        out += xy[i]^4 - 16*xy[i]^2 + 5*xy[i]
    end
    out *= 0.5
    return out
end



using LatinHypercubeSampling
using ScatteredInterpolation
using BlackBoxOptim
using Parameters
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
include("rbf_hypers_opt_utilities.jl")
include("sample_infill.jl")
include("scaling.jl")
include("smoptimize.jl")
include("search_range.jl")


#TODO
#####Tracing
#####Smoothing - Optimize the smoothness level 
#####Result type containing
        # Original samples
        # Infill points
        # Smallest evaluated function value per iteration
        # Options?
        # Num iterations
        # ?

#####Infill
    #Verify that the infill works as expected with small test function
    #Assert that the infill points are not the same as the samples
    #Handle several infill methods
    #Pareto front type of infill? 
    #How to define the distance between samples of a pareto front?
    #New designs.
    #1. Pure Exploration of design space. 
#Done   #1. Infill distance, find largest empty space - GA to find this, KD-tree to increase speed
        #2. Think of dimension scaling. Not possible with several interpolants?
#Done   #3. Find maximum STD
        #4. Large gradient
    #2. Pure Exploitation
#Done   #1. Take overall min prediction
#Done   #2. Take min median prediction
    #3. Both
#Done   #1. min(Median prediction - 1STD)
#Done   #2. min(Median prediction - 2STD)
        #3. min(Median prediction) >= c Distance from any other point
        #4. min(Median prediction - xSTD) where X is the smallest value acheiving distance
        #   >= c from any other point

            
        

end # module
