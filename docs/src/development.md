# Development

## Tests and coverage

From the package directory:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test(coverage=true)'
```

The suite exercises first-order and exact modes, bounds, scalar/vector and mixed
constraints, custom Hessians, changing derivative nonzero values, structured
variables, JuMP derivative routing, objective scaling, presets, option types,
termination, validation errors, and workspace reuse.

CI tests Julia 1.10 and current Julia on Linux, plus current Julia on macOS and
Windows. Coverage is produced as an artifact. A separate job builds the docs
and executes their examples. Successful `main` docs builds deploy to GitHub Pages;
pull requests produce downloadable documentation artifacts.

## Build documentation

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

If using an unregistered local NonconvexCore revision, develop its checkout in
the docs environment at the same time. Open `docs/build/index.html` after the
build. `docs/fragment.toml` supplies page metadata for the sibling packages'
DocumenterFragments-based documentation assembly.

## Formatting

The package uses JuliaFormatter's Blue style, configured in
`.JuliaFormatter.toml`. With JuliaFormatter installed in a tooling environment:

```julia
using JuliaFormatter
format("path/to/NonconvexUno.jl")
```

## Releases

Use JuliaRegistrator and TagBot: bump `Project.toml` and `NEWS.md`, merge the
change, and trigger registration on that commit. Once the registry update
merges, TagBot creates the version tag and release. The included workflow
supports that flow after the package has a GitHub repository and registration.
