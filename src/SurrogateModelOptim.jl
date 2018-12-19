module SurrogateModelOptim

export  smoptimize,
        sm_hypers_opt

using LatinHypercubeSampling: LHCoptim
using ScatteredInterpolation: interpolate
using BlackBoxOptim: bboptimize
using Parameters: reconstruct, type2dict, @with_kw, @unpack
using Distances
using Statistics

include("default_parameters.jl")
include("rbf_hypers_opt.jl")
include("sample_infill.jl")
include("scaling.jl")
include("setup.jl")
include("smoptimize.jl")

end # module
