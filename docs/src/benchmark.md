# Benchmark
Benchmark results from the paper are available in the folder
```julia
julia> joinpath(dirname(pathof(SurrogateModelOptim)), "../benchmark")
```
There is also a script for convenience to run the benchmark yourself. The solver
settings used for the benchmark can be found with
```julia
julia> SurrogateModelOptim.paper_bench_opts()
```
The benchmark options are similar to the default options. To update the benchmark
options, for example to turn off parallel creation of the surrogate model, simply 
update the default benchmark options with the new field as 

```julia
julia> SurrogateModelOptim.Options(SurrogateModelOptim.paper_bench_opts(); parallel_surrogate=false)
```

The file 
```julia
julia> joinpath(dirname(pathof(SurrogateModelOptim)), "../examples/test_functions.jl"))
```
contains several additional benchmark functions which were not tested in the paper.
No guarantee is given that the implementation of the functions is correct.

If this package is used in work leading to publication, please consider citing the paper.