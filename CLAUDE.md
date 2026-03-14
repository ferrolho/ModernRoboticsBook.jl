# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Julia implementation of algorithms from the *Modern Robotics: Mechanics, Planning, and Control* textbook by Kevin Lynch and Frank Park. The package provides functions for rigid-body motions, forward/inverse kinematics, dynamics, trajectory generation, and robot control.

## Commands

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run doctests only
julia --project docs/make.jl

# Start Julia REPL with the package loaded
julia --project -e 'using ModernRoboticsBook'

# Format code (default style)
julia -e 'using JuliaFormatter; format(".")'
```

JuliaFormatter is installed globally (not a project dependency). To install it: `julia -e 'using Pkg; Pkg.add("JuliaFormatter")'`. See https://domluna.github.io/JuliaFormatter.jl/stable/ for configuration options.

## Architecture

- **Single-module, single-file package**: All code lives in `src/ModernRoboticsBook.jl` — one module with ~50 exported functions organized by textbook chapter (3, 4, 5, 6, 8, 9, 11).
- **Only stdlib dependency**: `LinearAlgebra` (aliased as `LA`).
- **Tests**: `test/runtests.jl` includes [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality checks alongside chapter-organized test sets.
- **Docs**: Built with Documenter.jl via `docs/make.jl`.

## Conventions

- **Function naming**: PascalCase (e.g., `MatrixExp3`, `FKinBody`, `TransInv`) to match the textbook's notation. The only lowercase export is `ad`.
- **Type signatures**: Functions accept `AbstractVector`/`AbstractMatrix` — no custom types are defined.
- **Docstrings**: Include `jldoctest` examples with `setup = :(using ModernRoboticsBook)`. These are tested via Documenter.
