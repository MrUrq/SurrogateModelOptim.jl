using SurrogateModelOptim
using Test

include("samples_infill.jl")
include("examples.jl")




# @testset "Interface" begin
#     @testset "Dimensional scaling" begin
#         @test SurrogateModelOptim._update_options([(-1.0,1.0),(1.0,2.0)]).variable_dim_scaling == true
#         @test SurrogateModelOptim._update_options([(-1.0,1.0)]).variable_dim_scaling == false
#     end

#     @testset "Rule of thumb initial sample number" begin
#         n_samples = 3
#         rule_of_thumb = 10
#         @test SurrogateModelOptim._update_options(repeat([(-1.0,1.0)],n_samples)).num_start_samples == n_samples*rule_of_thumb
#     end
# end

# @testset "Sampling plan" begin
#     search_range = [(-2.0,1.0),(-3.0,0.5)]
#     num_start_samples = 5
#     sampling_plan_opt_gens = 100
#     @test extrema(SurrogateModelOptim._LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)[1,:]) == search_range[1]
#     @test extrema(SurrogateModelOptim._LHC_sampling_plan(search_range,num_start_samples,sampling_plan_opt_gens)[2,:]) == search_range[2]
# end