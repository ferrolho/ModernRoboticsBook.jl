using Documenter, ModernRobotics

makedocs(;
    modules=[ModernRobotics],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "man/examples.md",
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
        ],
    ],
    repo="https://github.com/ferrolho/ModernRobotics.jl/blob/{commit}{path}#L{line}",
    sitename="ModernRobotics.jl",
    authors="Henrique Ferrolho",
    analytics = "UA-72743607-4",
    assets=[],
)

deploydocs(;
    repo="github.com/ferrolho/ModernRobotics.jl",
)
