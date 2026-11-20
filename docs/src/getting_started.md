# Getting started

## Install locally

From the parent Nonconvex workspace:

```julia
using Pkg
Pkg.develop([
    PackageSpec(path="repos/NonconvexCore.jl"),
    PackageSpec(path="repos/NonconvexUno.jl"),
])
```

## Solve a constrained problem

Nonconvex inequality constraints use the convention `c(x) <= rhs`, with a
default right-hand side of zero. Equality constraints use `c(x) == rhs`.

```@example start
using NonconvexUno

model = Model(x -> sqrt(x[2]))
addvar!(model, [1e-4, 1e-4], [Inf, Inf])
add_ineq_constraint!(model, x -> (2x[1])^3 - x[2])
add_ineq_constraint!(model, x -> (1-x[1])^3 - x[2])
result = optimize(model, UnoAlg(), [1.234, 2.345];
    options=UnoOptions(logger="SILENT"))
@assert isapprox(result.minimizer, [1/3, 8/27]; atol=1e-5)
(minimizer=result.minimizer, minimum=result.minimum, status=result.status)
```

To use exact Hessians, change the options to
`UnoOptions(first_order=false, logger="SILENT")`.

## Vector constraints

```julia
add_ineq_constraint!(model,
    FunctionWrapper(x -> [x[1] - 5, x[2] - 5], 2))
```

## Structured variables

```@example structured
using NonconvexUno
model = DictModel(d -> sum(abs2, d[:x] .- [1.0, 2.0]))
addvar!(model, :x, [-Inf, -Inf], [Inf, Inf]; init=[0.0, 0.0])
result = optimize(model, UnoAlg(); options=UnoOptions(logger="SILENT"))
@assert isapprox(result.minimizer[:x], [1.0, 2.0]; atol=1e-5)
result.minimizer
```
