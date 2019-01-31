module SurrogateModelOptim

export  smoptimize,
        model_infill,
        TestFunction


using LatinHypercubeSampling
using ScatteredInterpolation
using BlackBoxOptim
using Distributions
using Parameters
using Distances
using Statistics
using StatsBase
using LinearAlgebra
using Distributed
using StaticArrays
using NearestNeighbors
using Printf

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



#####Fancy results printing
#####Call the method with the results type already?

#####Change distance to negative instead of 1/x? same with std?

#####Does std samples, possibly constrain the deviation in other locations? Fit an elephant, put in paper.

#####Call the method with a test function + options shortcut

#####Create plot of optim method vs optim method, shaded area thing. 

#####Find which parameters to use when optimizing, which infills? how many interpolants? rippa? 
#####Use std and min
#####Find which of the interpolants to use, A,B,C or D? Dependent on problem? Rosen2 hart6? Initial vs final convergence? How does this work with integer test case?


#####PlotlyJS parallel coordinates
#####Plot accuracy as a function of "time" predicted and actual

#####Performance of several infill points at once vs one at the time?



#####Infill
    #Verify that the infill works as expected with small test function. 
    #Plot it and show the where the points are added after each iteration
    #Pareto front infill option possibly?
    
    #1. Pure Exploration of design space. 
#Done   #1. Infill distance, find largest empty space - GA to find this, KD-tree to increase speed
        #2. Think of dimension scaling. Not possible with several interpolants?   - Outside of project scope
#Done   #3. Find maximum STD
    #2. Pure Exploitation
#Done   #1. Take overall min prediction
#Done   #2. Take min median prediction
    #3. Both
#Done   #1. min(Median prediction - 1STD)
#Done   #2. min(Median prediction - 2STD)
        #3. min(Median prediction) >= c Distance from any other point
        #4. min(Median prediction - xSTD) where X is the smallest value achieving distance
        #   >= c from any other point

#####Integer test case
#####Documentation
#####Examples
#####Function headers

end # module
