using Documenter
using NonconvexUno

makedocs(;
    modules=[NonconvexUno],
    sitename="NonconvexUno.jl",
    authors="Mohamed Tarek and contributors",
    repo="github.com/JuliaNonconvex/NonconvexUno.jl.git",
    checkdocs=:all,
    doctest=true,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true", edit_link=nothing
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Options" => "options.md",
        "JuMP integration" => "jump.md",
        "Workspaces and results" => "workspaces.md",
        "API" => "api.md",
        "Development" => "development.md",
    ],
)
