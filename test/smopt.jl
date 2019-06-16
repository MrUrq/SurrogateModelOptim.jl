@testset "smoptimize" begin
   
    search_range = [(-10.0,10.0), (0.0,2.0), (25.0,40.0)]
    n_samples = 6
    plan = SurrogateModelOptim.scaled_LHC_sampling_plan(search_range,n_samples,1000;trace=false)
    samples = mapslices(sum,plan,dims=1)

    #Make sure correct input format is used
    Test.@test_throws ErrorException surrogate_model(permutedims(plan),samples)
end