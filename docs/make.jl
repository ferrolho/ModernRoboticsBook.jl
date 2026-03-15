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

# Truncate floats after 10 decimal digits to handle cross-platform precision differences
makedocs(
    modules = [ModernRoboticsBook],
    doctestfilters = [r"(\d\.\d{10})\d+" => s"\1"],
    format = Documenter.HTML(
        canonical = "https://ferrolho.github.io/ModernRoboticsBook.jl/stable",
        description = "Julia implementation of the Modern Robotics textbook algorithms for kinematics, dynamics, trajectory generation, and robot control",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "man/examples.md",
        "API Reference" => [
            "lib/robot-model.md",
            "lib/rigid-body-motions.md",
            "lib/forward-kinematics.md",
            "lib/velocity-kinematics.md",
            "lib/inverse-kinematics.md",
            "lib/dynamics.md",
            "lib/trajectory-generation.md",
            "lib/robot-control.md",
        ],
    ],
    repo = "https://github.com/ferrolho/ModernRoboticsBook.jl/blob/{commit}{path}#L{line}",
    sitename = "ModernRoboticsBook.jl",
    authors = "Henrique Ferrolho",
)

deploydocs(repo = "github.com/ferrolho/ModernRoboticsBook.jl")
