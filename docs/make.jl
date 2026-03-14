using Documenter, ModernRoboticsBook
import LinearAlgebra as LA

DocMeta.setdocmeta!(
    ModernRoboticsBook,
    :DocTestSetup,
    quote
        using ModernRoboticsBook
        import LinearAlgebra as LA
        # 3-link robot model for dynamics and control examples
        M01 = [1 0 0 0; 0 1 0 0; 0 0 1 0.089159; 0 0 0 1]
        M12 = [0 0 1 0.28; 0 1 0 0.13585; -1 0 0 0; 0 0 0 1]
        M23 = [1 0 0 0; 0 1 0 -0.1197; 0 0 1 0.395; 0 0 0 1]
        M34 = [1 0 0 0; 0 1 0 0; 0 0 1 0.14225; 0 0 0 1]
        G1 = LA.Diagonal([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7])
        G2 = LA.Diagonal([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393])
        G3 = LA.Diagonal([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275])
        link_frames = [M01, M12, M23, M34]
        spatial_inertias = [G1, G2, G3]
        screw_axes = [1 0 1 0 1 0; 0 1 0 -0.089 0 0; 0 1 0 -0.089 0 0.425]'
    end;
    recursive = true,
)

makedocs(
    modules = [ModernRoboticsBook],
    format = Documenter.HTML(analytics = "UA-72743607-4", assets = []),
    pages = [
        "Home" => "index.md",
        "Examples" => "man/examples.md",
        "API Reference" => [
            "lib/ch3.md",
            "lib/ch4.md",
            "lib/ch5.md",
            "lib/ch6.md",
            "lib/ch8.md",
            "lib/ch9.md",
            "lib/ch11.md",
        ],
    ],
    repo = "https://github.com/ferrolho/ModernRoboticsBook.jl/blob/{commit}{path}#L{line}",
    sitename = "ModernRoboticsBook.jl",
    authors = "Henrique Ferrolho",
)

deploydocs(repo = "github.com/ferrolho/ModernRoboticsBook.jl")
