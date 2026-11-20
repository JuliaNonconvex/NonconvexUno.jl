# NonconvexUno.jl

[![CI](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions/workflows/Documentation.yml)

A Nonconvex wrapper for [Uno](https://github.com/cvanaret/Uno), using the
registered `UnoSolver.jl` package and its precompiled solver binaries.

See the [documentation](https://JuliaNonconvex.github.io/NonconvexUno.jl/) for
usage, options, JuMP integration, and API details. The [documentation sources](docs/src/index.md)
and downloadable workflow artifacts are also available.

## Install from GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaNonconvex/NonconvexUno.jl", rev="main")
```

## Local installation

From this workspace's root, in your chosen Julia environment:

```julia
using Pkg
Pkg.develop([
    PackageSpec(path="repos/NonconvexCore.jl"),
    PackageSpec(path="repos/NonconvexUno.jl"),
])
```

## Usage

```julia
using NonconvexUno

model = Model(x -> sqrt(x[2]))
addvar!(model, [1e-4, 1e-4], [10.0, 10.0])
add_ineq_constraint!(model, x -> (2x[1])^3 - x[2])
add_ineq_constraint!(model, x -> (1-x[1])^3 - x[2])

result = optimize(model, UnoAlg(), [1.234, 2.345];
    options=UnoOptions(logger="SILENT"))
result.minimizer  # approximately [1/3, 8/27]
result.minimum   # approximately sqrt(8/27)
result.status    # :first_order on convergence
```

Like NonconvexIpopt, this package reexports NonconvexCore and provides algorithm,
options, workspace, and result types. Use `using NonconvexUno` directly.

### Options and derivatives

- `first_order=true` (default): Uno's `LBFGS` Hessian approximation.
- `first_order=false`: exact Lagrangian Hessians through NonconvexCore's AD.
- `preset="ipopt"` (default): Uno's interior-point preset. For SQP, use
  `preset="filtersqp", first_order=false` or `preset="funnelsqp", first_order=false`.
- `hessian_model`: explicitly override the above Hessian choice using Uno's name.
  For SLP, use `preset="filterslp", hessian_model="zero"`.
- `linear_constraints=true`: compute the constraint Jacobian once.
- Other keywords are forwarded to Uno, e.g. `max_iterations=1000`,
  `primal_tolerance=1e-8`, `dual_tolerance=1e-8`, `time_limit=60.0`.
  These are **Uno option names**, not Ipopt option names.

All native Uno options are available through keyword arguments; there is no
wrapper-side name whitelist. Uno validates option names, types, and values.
Use strings for named strategies, booleans for switches, integers for counts,
and floating-point values for tolerances and time limits. See the
[Uno documentation](https://unosolver.readthedocs.io) for the full option reference.

Ordinary function models use dense fixed Jacobian/lower-triangular Hessian
patterns, including entries that vanish at the initial point. It supports scalar
and vector constraints, equality/inequality right-hand sides, and infinite bounds.
Sparse/symbolic option switches and iteration callbacks are not implemented;
passing a callback raises an error. Integer and semidefinite constraints are not
supported by this solver interface.

### JuMP models

```julia
import JuMP
jm = JuMP.Model()
JuMP.@variable(jm, x[1:2], start=0.5)
JuMP.@objective(jm, Min, (x[1]-2)^2 + (x[2]-1)^2)
JuMP.@constraint(jm, x[1] + x[2] == 2)
model = NonconvexCore.DictModel(jm)
result = optimize(model, UnoAlg(); options=UnoOptions(logger="SILENT"))
```

Converted JuMP models use the shared JuMP/MOI evaluator directly for objective
gradients, constraint values, sparse Jacobians, and (when requested) sparse
Lagrangian Hessians. No fallback AD is used for these derivatives. The test suite
disables all fallback derivative entry points during JuMP solves to verify this.

For an analytic Hessian, pass `lagrangian_hessian=(x, sigma, lambda) -> H` to
`optimize` together with `first_order=false`. Return the Hessian of
`sigma*f(x) + dot(lambda, c(x))`, with inequalities first and equalities second.

`result.result` contains Uno's full execution statistics, including multipliers,
feasibility residuals, iteration count, and native status codes. `result.fcalls`
is Uno's objective evaluation count. `result.problem` is the native model.

### Reusing a workspace

```julia
vm, x0, _ = NonconvexCore.tovecmodel(model, [1.234, 2.345])
using NonconvexCore: optimize!, reset!
w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(logger="SILENT"))
r = optimize!(w)
reset!(w, r.minimizer)
r = optimize!(w)
```

Each solve starts at `w.x0`; use `reset!` to change it. Model structure and
derivative settings are fixed when constructing the workspace.

## Tests

```julia
using Pkg
Pkg.test("NonconvexUno")
# Instrument the package for line coverage:
Pkg.test("NonconvexUno"; coverage=true)
```

## Build the documentation

From the package directory:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html`. See [CONTRIBUTING.md](CONTRIBUTING.md) for development
and release details.
