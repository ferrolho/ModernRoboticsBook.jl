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

## Pre-commit Checks

Before committing, always run these three checks and fix any failures:

```bash
# 1. Formatter — must return true (no changes needed)
julia -e 'using JuliaFormatter; @assert format(".") "Formatter made changes — stage them and retry"'

# 2. Doctests and doc build — must exit cleanly (also catches missing docstrings in docs pages)
julia --project docs/make.jl

# 3. Full test suite
julia --project -e 'using Pkg; Pkg.test()'
```

Common issues these catch:
- **Formatter**: unformatted code after edits.
- **Doc build**: new exported symbols missing from `docs/src/lib/*.md` pages (the `Pages` filter in `@autodocs` blocks must match the source file where the docstring lives, e.g., `robot.jl` not `ModernRoboticsBook.jl`).
- **Tests**: regressions from renames or refactors.

## Architecture

- **Single-module package**: Core algorithms in `src/ModernRoboticsBook.jl`, robot model loading in `src/robot.jl` (included via `include("robot.jl")`).
- **Dependencies**: `LinearAlgebra` (aliased as `LA`) and `JSON`.
- **Tests**: `test/runtests.jl` includes [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality checks alongside chapter-organized test sets.
- **Docs**: Built with Documenter.jl via `docs/make.jl`. Each `docs/src/lib/*.md` page uses `@autodocs` with a `Pages` filter — when moving code between source files, update the corresponding doc page filter.

## Conventions

- **Function naming**: snake_case (e.g., `matrix_exp3`, `forward_kinematics_body`, `transform_inv`) following Julia conventions.
- **Type signatures**: Functions accept `AbstractVector`/`AbstractMatrix` — the `Robot` struct is the only custom type.
- **Docstrings**: Include `jldoctest` examples with `setup = :(using ModernRoboticsBook)`. These are tested via Documenter.
- **Terminology**: "configuration" refers to joint positions, "pose" or "transform" refers to SE(3) task-space quantities.
