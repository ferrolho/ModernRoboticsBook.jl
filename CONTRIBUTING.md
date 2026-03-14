# Contributing to ModernRoboticsBook.jl

## Prerequisites

- Julia 1.10 or later
- [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) (installed globally):
  ```bash
  julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
  ```

## Running tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Building docs locally

The docs project intentionally does not list `ModernRoboticsBook` in its `Project.toml` — this ensures you always build against the local checkout rather than the registered version. You need to add it manually with `Pkg.develop`:

```bash
julia --project=docs -e '
  using Pkg
  Pkg.develop(PackageSpec(path=pwd()))
  Pkg.instantiate()'
```

Then build (this also runs doctests):

```bash
julia --project=docs docs/make.jl
```

The generated site is in `docs/build/` — open `docs/build/index.html` to preview.

## Running doctests only

Doctests are run as part of the docs build above. There is no separate doctest-only command.

## Formatting

All code must be formatted with JuliaFormatter (default style). CI will reject unformatted code.

```bash
julia -e 'using JuliaFormatter; format(".")'
```

## Releasing a new version

1. Update the `version` field in `Project.toml`
2. Commit and push to `master`
3. Invoke [JuliaRegistrator](https://github.com/JuliaRegistries/Registrator.jl) by commenting `@JuliaRegistrator register` on the commit (with release notes)
4. The registry PR auto-merges after a waiting period
5. [TagBot](https://github.com/JuliaRegistries/TagBot) automatically creates the git tag and GitHub release

**Important:** Do not create git tags manually — TagBot must create them so that it also creates the GitHub release.
