using Documenter, SurrogateModelOptim

makedocs(;
        modules = [SurrogateModelOptim],
        format = Documenter.HTML(   assets = ["assets/favicon.ico"],
                                prettyurls = get(ENV, "CI", nothing) == "true"),
        sitename = "SurrogateModelOptim.jl",
        strict = true,
        clean = true,
        checkdocs = :none,
        pages = [
                "Home" => "index.md",
                "Manual" => [
                        "smoptimize.md",
                        "surrogate_model.md",
                        "model_infill.md",
                        "options.md",
                        "examples.md",
                        "benchmark.md",
                ]
        ]
)

deploydocs(
    repo = "github.com/MrUrq/SurrogateModelOptim.jl.git"
)