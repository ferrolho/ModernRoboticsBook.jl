using Documenter, ModernRobotics

makedocs(;
    modules=[ModernRobotics],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/ferrolho/ModernRobotics.jl/blob/{commit}{path}#L{line}",
    sitename="ModernRobotics.jl",
    authors="Henrique Ferrolho",
    assets=[],
)

deploydocs(;
    repo="github.com/ferrolho/ModernRobotics.jl",
)
