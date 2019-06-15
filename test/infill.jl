@testset "infill" begin
    using Statistics

    #KD-tree for the plan as well as the samples from the plan
    search_range = [(-10.0,10.0), (0.0,2.0), (25.0,40.0)]
    n_samples = 6
    plan = SurrogateModelOptim.scaled_LHC_sampling_plan(search_range,n_samples,1000;trace=false)
    samples = mapslices(sum,plan,dims=1)

    combined_points_samples = [plan;samples]
    kdtree = SurrogateModelOptim.KDTree(combined_points_samples)  

    @test SurrogateModelOptim._nearest_point(kdtree,x->sum(x),plan[:,1]) == 0
    @test SurrogateModelOptim.distance_infill(plan,samples,x->sum(x))(plan[:,1]) == 0
    @test SurrogateModelOptim.minimum_infill(x->(x, 2x, 4x))(1) == 1
    @test SurrogateModelOptim.median_infill(x->(x, 2x, 4x))(1) == 2
    @test SurrogateModelOptim.mean_infill(x->(x, 2x, 4x))(1) == mean([1,2,4])
    @test SurrogateModelOptim.std_infill(x->(x, 2x, 4x))(1) == -1*std([1 2 4])

    #test for catch block showing correct exception
end