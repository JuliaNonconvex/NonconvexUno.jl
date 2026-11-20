# Contributing

## Development environment

```sh
git clone https://github.com/JuliaNonconvex/NonconvexUno.jl
cd NonconvexUno.jl
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test(coverage=true)'
```

Tests use JuMP as a test-only dependency. The package itself depends on
NonconvexCore and UnoSolver, and defaults to first-order L-BFGS. Keep that
default when adding solver features.

## Derivative regression tests

For ordinary function models, check derivative coordinates and signs against
analytic values, including entries that are zero at the initial point. For
JuMP models, verify that all derivative callbacks use the shared evaluator;
the tests prohibit fallback AD during complete converted-model solves.

## Documentation

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build checks all package docstrings and runs the examples. The generated
site is in `docs/build/`. The documentation workflow uploads it as an artifact
and deploys successful `main` builds to GitHub Pages. Pull requests only build
the documentation and upload an artifact.

## Formatting and generated files

Use JuliaFormatter with `.JuliaFormatter.toml` (Blue style). Do not commit
manifests, coverage output, documentation builds, or local dependency paths.

## Releases

1. Check the latest registered version and choose the next valid SemVer version.
2. Update `Project.toml` and `NEWS.md`, run tests and docs, and merge to `main`.
3. Comment `@JuliaRegistrator register` on the commit, including a `Release notes:` block.
4. Address the General registry checks and wait for registration to merge.
5. Let TagBot create the tag and release; do not manually tag while TagBot is active.

The initial release also requires registering the package itself in General.
