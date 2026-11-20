# Workspaces and results

## Reuse a workspace

```@example workspace
using NonconvexUno
using NonconvexCore: optimize!, reset!

model = Model(x -> sum(abs2, x .- 2))
addvar!(model, [-Inf, -Inf], [Inf, Inf])
vm, x0, _ = NonconvexCore.tovecmodel(model, [0.0, 0.0])
w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(logger="SILENT"))
result = optimize!(w)
reset!(w, result.minimizer)
result = optimize!(w)
@assert isapprox(result.minimizer, [2.0, 2.0]; atol=1e-5)
result.minimizer
```

Every solve starts at `w.x0`. Use `reset!` to change it; `reset!(w)` leaves the
initial point unchanged. Each result owns a copy of its minimizer. Rebuild the
workspace when changing model structure or derivative settings. Workspaces
operate on flattened vectors; `optimize` restores structured minimizers for you.

## Inspect results

[`NonconvexUno.UnoResult`](@ref) provides:

| Field | Meaning |
|:--|:--|
| `minimizer` | Final variables, restored to the original structure by `optimize` |
| `minimum` | Final objective value |
| `status` | UnoSolver's symbolic termination status |
| `fcalls` | Uno's objective evaluation count |
| `problem` | Native Uno model |
| `result` | Full `UnoExecutionStats`, retaining the native model and solver |

`status == :first_order` indicates convergence to a first-order stationary
point, including when exact Hessians were used. Other statuses include
`:max_iter`, `:max_time`, `:infeasible`, `:small_step`, `:unbounded`, and `:exception`.
Inspect `result.result.optimization_status` and `result.result.solution_status`
for native status codes.

The native statistics also expose multipliers, primal/dual feasibility,
complementarity, iteration count, elapsed time, and derivative evaluation counts.
