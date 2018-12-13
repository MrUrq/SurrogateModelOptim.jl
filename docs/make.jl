using Documenter, SurrogateModelOptim

makedocs(;
    modules=[SurrogateModelOptim],
    format=:html,
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MrUrq/SurrogateModelOptim.jl/blob/{commit}{path}#L{line}",
    sitename="SurrogateModelOptim.jl",
    authors="Magnus Urquhart",
    assets=[],
)

deploydocs(;
    repo="github.com/MrUrq/SurrogateModelOptim.jl",
    target="build",
    julia="1.0",
    deps=nothing,
    make=nothing,
)
