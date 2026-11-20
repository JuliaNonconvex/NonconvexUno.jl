# v0.1.0

- Initial Nonconvex interface to Uno through UnoSolver.jl.
- First-order L-BFGS by default; optional exact or user-supplied Lagrangian Hessians.
- Nonlinear equality/inequality constraints, vector constraints, and variable bounds.
- Native Uno options and execution statistics, reusable workspaces, and solver capability declarations.
- Direct JuMP/MOI evaluator gradients, sparse Jacobians, and sparse Lagrangian Hessians.
- Executable documentation, derivative-routing regression tests, coverage, and CI workflows.
