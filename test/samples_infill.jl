@testset "infill" begin
    
    #KD-tree for the plan as well as the samples from the plan
    search_range = [(-10.0,10.0), (0.0,2.0), (25.0,40.0)]
    n_samples = 6
    plan = SurrogateModelOptim._LHC_sampling_plan(search_range,n_samples,1000;trace=true)
    samples = mapslices(sum,plan,dims=1)

    combined_points_samples = [plan;samples]
    kdtree = SurrogateModelOptim.KDTree(combined_points_samples)  

    @test SurrogateModelOptim._nearest_point(kdtree,x->sum(x),plan[:,1]) == 0
    
end