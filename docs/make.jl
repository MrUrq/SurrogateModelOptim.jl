using Documenter, SurrogateModelOptim

makedocs(;
    modules=[SurrogateModelOptim],
    format = Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MrUrq/SurrogateModelOptim.jl/blob/{commit}{path}#L{line}",
    sitename="SurrogateModelOptim.jl",
    authors="Magnus Urquhart",
    assets = ["assets/favicon.ico"],
)

deploydocs(;
    repo="github.com/MrUrq/SurrogateModelOptim.jl",
    target="build",
    deps=nothing,
    make=nothing,
)
