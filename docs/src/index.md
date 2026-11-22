# NonconvexUno

NonconvexUno exposes [Uno](https://github.com/cvanaret/Uno) through the Nonconvex
modeling interface. It uses the registered `UnoSolver.jl` package and its
precompiled binaries.

The default is **first-order L-BFGS**, using Uno's `ipopt` preset. Exact Hessians
are opt-in. Every native Uno option can be supplied through [`UnoOptions`](@ref).

Supported models have continuous variables, finite or infinite bounds, and
nonlinear equality and inequality constraints. Scalar and vector constraints
and Nonconvex's structured variables are supported. Converted JuMP models use
their evaluator's analytic/sparse derivatives directly.

Start with [Getting started](@ref), then see [Options](@ref) and
[JuMP integration](@ref).

## Scope

- Ordinary Julia function models currently use dense derivative structures.
- JuMP-converted models use the evaluator's sparse structures.
- Sparse/symbolic switches for ordinary function models are not implemented.
- Iteration callbacks are not exposed; supplying a callback raises an error.
- Integer/binary variables and semidefinite constraints are unsupported.

The package is loaded directly with `using NonconvexUno`.
