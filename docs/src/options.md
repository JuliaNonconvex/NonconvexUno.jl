# Options

[`UnoOptions`](@ref) accepts every native Uno option via keyword arguments;
option names are not restricted to a wrapper-side whitelist. Uno validates
names, types, and values. Consult the [Uno reference](https://unosolver.readthedocs.io)
for the complete native option descriptions.

## Wrapper settings

| Keyword | Default | Meaning |
|:--|:--|:--|
| `first_order` | `true` | `true` maps to `hessian_model="LBFGS"`; `false` maps to `"exact"` |
| `linear_constraints` | `false` | Cache the constraint Jacobian at the start |
| `preset` | `"ipopt"` | Uno strategy preset |
| `hessian_model` | `nothing` | Explicit override of the `first_order` mapping |

`first_order` is only a convenience mapping onto Uno's native `hessian_model`.
For other Hessian models, set `hessian_model` directly: `"LBFGS"`, `"LSR1"`,
`"identity"`, or `"zero"`. An explicit `hessian_model` takes precedence over
`first_order`. Presets are applied before individual options, so explicit native
settings override them. `linear_constraints=true` is appropriate only when
**all** constraints are linear.

Inconsistent combinations are rejected when the options are constructed:
Uno's L-SR1 quasi-Newton model is not supported with the interior-point method
(the `"ipopt"` or `"auto"` preset, or an explicit
`inequality_handling_method="interior_point"`) because the barrier line search
fails. Pair `"LSR1"` with the SQP presets (`"filtersqp"`, `"funnelsqp"`) or an
explicit `inequality_handling_method="inequality_constrained"` instead.

```julia
options = UnoOptions(
    max_iterations=1000,
    time_limit=60.0,
    primal_tolerance=1e-8,
    dual_tolerance=1e-8,
    quasi_newton_memory_size=10,
    linear_solver="MUMPS",
    print_solution=false,
    logger="SILENT",
)
```

Use strings for strategy names, booleans for switches, integers for counts,
and floating-point numbers for tolerances and time limits. Floating-point
types such as `Float32` are normalized to `Float64`.

## Presets

```julia
UnoOptions()                                       # L-BFGS interior point
UnoOptions(first_order=false)                      # exact-Hessian interior point
UnoOptions(preset="filtersqp", first_order=false)  # filter SQP
UnoOptions(preset="funnelsqp", first_order=false)  # funnel SQP
UnoOptions(preset="filtersqp", hessian_model="LSR1") # filter SQP with L-SR1
UnoOptions(preset="filterslp", hessian_model="zero") # filter SLP
```

Native option names differ from Ipopt's: for example, use `max_iterations`,
not `max_iter`. Available linear/QP solvers depend on the Uno binary build.

## Custom Hessian

Pass `lagrangian_hessian=(x, sigma, lambda) -> H` to `optimize` or `Workspace`,
with `first_order=false`. Return the Hessian of
`sigma*f(x) + dot(lambda, c(x))`. Constraint rows are ordered inequalities first,
then equalities; residuals already subtract their right-hand sides. A custom
Hessian takes precedence over JuMP's evaluator Hessian too.

First-order mode never calls this Hessian callback; providing one while the
effective Hessian model is not `"exact"` logs a warning from `Workspace`.
