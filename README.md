<p align="center">
  <img src="assets/logo.png" width="400" alt="NonconvexUno.jl logo">
</p>

[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://julianonconvex.github.io/NonconvexUno.jl/stable)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://julianonconvex.github.io/NonconvexUno.jl/dev)
[![Build Status](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaNonconvex/NonconvexUno.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaNonconvex/NonconvexUno.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNonconvex/NonconvexUno.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-Blue-4473c9.svg)](https://github.com/invenia/BlueStyle)
[![JET](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

# NonconvexUno.jl

Unified nonlinear optimization with Uno

First-order (L-BFGS) and exact-Hessian interior point and SQP optimization for nonlinearly constrained problems, with direct JuMP evaluator derivatives.

NonconvexUno.jl is part of the [JuliaNonconvex](https://github.com/JuliaNonconvex) ecosystem built on [NonconvexCore.jl](https://github.com/JuliaNonconvex/NonconvexCore.jl) under the [Nonconvex.jl](https://github.com/JuliaNonconvex/Nonconvex.jl) umbrella.

## Quick start

```julia
using NonconvexUno

alg = UnoAlg()
options = UnoOptions()
result = optimize(model, alg, x0, options = options)
```

## Documentation

The full documentation is available at [https://julianonconvex.github.io/NonconvexUno.jl/stable](https://julianonconvex.github.io/NonconvexUno.jl/stable) (stable) and [https://julianonconvex.github.io/NonconvexUno.jl/dev](https://julianonconvex.github.io/NonconvexUno.jl/dev) (development).

## Installation

The package is not registered in the General registry yet. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaNonconvex/NonconvexUno.jl", rev = "main")
```

Once registered, install it with:

```julia
using Pkg
Pkg.add("NonconvexUno")
```
