using Documenter, ModernRoboticsBook

makedocs(
    modules = [ModernRoboticsBook],
    format = Documenter.HTML(
        analytics = "UA-72743607-4",
        assets = [],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/examples.md",
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
        ],
    ],
    repo = "https://github.com/ferrolho/ModernRoboticsBook.jl/blob/{commit}{path}#L{line}",
    sitename = "ModernRoboticsBook.jl",
    authors = "Henrique Ferrolho",
)

deploydocs(
    repo = "github.com/ferrolho/ModernRoboticsBook.jl",
)
