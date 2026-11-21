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
Windows. Coverage is uploaded to Codecov. A separate job builds the docs
and executes their examples. Successful builds deploy to GitHub Pages under
`dev` on `main`, `stable` on tags, and `preview/PRnn` for pull requests.

The remaining workflows mirror the other JuliaNonconvex packages:

- `Aqua` quality assurance and a blocking `JET.test_package` run in the test suite.
- A non-blocking JET report posted as a pull request comment (`test/jet_report.jl`).
- A JuliaFormatter Blue-style check that posts a comment on pull requests, and a
  nightly run that opens an automatic formatting pull request.
- CompatHelper bumping compat bounds and TagBot creating release tags after
  JuliaRegistrator registrations.

## Build documentation

The documentation is built as a [DocumenterFragments](https://github.com/PumasAI/DocumenterFragments.jl)
fragment: `docs/fragment.toml` declares the pages, and `docs/make.jl` calls
`build_fragment`, so the pages assemble into the org-wide Nonconvex.jl
documentation site. Cross-references to sibling packages only resolve in that
assembled site; the standalone build tolerates them as broken links. The
repository's GitHub Pages deploy serves this standalone preview under
`dev`/`stable`/`preview/PRnn`, and doctests run in every build.

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` after the build.

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
