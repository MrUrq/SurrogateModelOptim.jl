module SurrogateModelOptim

export  smoptimize

import LatinHypercubeSampling
import ScatteredInterpolation
import BlackBoxOptim
import Parameters

include("smoptimize.jl")
include("setup.jl")
include("sample_infill.jl")
include("default_parameters.jl")
include("scaling.jl")

end # module
