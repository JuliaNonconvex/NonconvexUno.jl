# JuMP integration

Load JuMP to activate NonconvexCore's conversion extension, and construct a
`DictModel` from the JuMP model:

```@example jump
using NonconvexUno
import JuMP

jm = JuMP.Model()
JuMP.@variable(jm, x[1:2], start=0.5)
JuMP.@objective(jm, Min, x[1]^2 + x[2]^2)
JuMP.@constraint(jm, x[1] + x[2] == 2)
model = NonconvexCore.DictModel(jm)
result = optimize(model, UnoAlg(); options=UnoOptions(logger="SILENT"))
@assert isapprox(result.minimum, 2.0; atol=1e-5)
result.minimizer
```

## Derivative routing

When the objective and constraints belong to the same complete converted
JuMP model, NonconvexUno calls NonconvexCore's shared JuMP/MOI evaluator for:

- Objective gradients (including objective-sense conversion and scaling).
- Constraint values and sparse Jacobian values/structure.
- Sparse Lagrangian Hessians when `first_order=false`.

It does not differentiate through JuMP using fallback AD, and it does not
allocate dense derivative structures for this path. The regression tests
replace fallback gradient, Jacobian, and Hessian entry points with throwing
functions, then solve JuMP models in both first- and second-order modes.
Results are compared with a direct JuMP–Uno solve.

If you add extra Nonconvex constraints after conversion, the complete-model
fast path is disabled to ensure those constraints are honored. First-order
mode can use the converted functions' derivative rules. Exact Hessians of
such mixed models require a suitable custom `lagrangian_hessian`, since nested
AD through the JuMP evaluator is unsupported.
