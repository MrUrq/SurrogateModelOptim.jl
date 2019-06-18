using Documenter, SurrogateModelOptim

makedocs(;
    modules=[SurrogateModelOptim],
    format = Documenter.HTML(   assets = ["assets/favicon.ico"],
                                prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
                "Home" => "index.md",
                "Manual" => [
                        "smoptimize.md",
                        "surrogate_model.md",
                        "options.md"
                ],
        ],
    repo="https://github.com/MrUrq/SurrogateModelOptim.jl/blob/{commit}{path}#L{line}",
    sitename="SurrogateModelOptim.jl",
    #authors="Magnus Urquhart",
    #strict = true,
    #clean = true,
    #checkdocs = :none,
)

deploydocs(;
    repo="github.com/MrUrq/SurrogateModelOptim.jl",
    target="build",
    deps=nothing,
    make=nothing,
)
