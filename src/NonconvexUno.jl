module NonconvexUno

using Reexport
@reexport using NonconvexCore
using LinearAlgebra: dot
using UnoSolver: UnoSolver
using NonconvexCore:
    VecModel,
    AbstractOptimizer,
    AbstractResult,
    NoCallback,
    getinit,
    getnineqconstraints,
    getneqconstraints,
    value_gradient,
    value_jacobian,
    value_hessian
import NonconvexCore: Workspace, optimize!, reset!

export UnoAlg, UnoOptions

"""
    UnoAlg()

Uno nonlinear optimizer, accessed through UnoSolver.jl. Supports continuous
variables, finite/infinite bounds, and nonlinear equality/inequality constraints.
Configure the algorithm with [`UnoOptions`](@ref).
"""
struct UnoAlg <: AbstractOptimizer end

function NonconvexCore.solver_capabilities(::UnoAlg)
    return ConstraintCapabilities(;
        linear_eq=true,
        nonlinear_eq=true,
        linear_ineq=true,
        nonlinear_ineq=true,
        linear_sd=false,
        nonlinear_sd=false,
        box=true,
        continuous=true,
        binary=false,
        integer=false,
        infinite_bounds=true,
    )
end

"""
    UnoOptions(; first_order=true, linear_constraints=false,
                 preset="ipopt", hessian_model=nothing, kwargs...)

`first_order=true` (the default) maps to `hessian_model="LBFGS"`: a
first-order, limited-memory quasi-Newton method. `first_order=false` maps to
`hessian_model="exact"`. An explicit `hessian_model` (`"exact"`, `"LBFGS"`,
`"LSR1"`, `"identity"`, or `"zero"`) overrides `first_order`.

`linear_constraints=true` caches the constraint Jacobian. Remaining keywords
are Uno options (for example `max_iterations`, `primal_tolerance`,
`dual_tolerance`, and `logger="SILENT"`). Other presets include `"filtersqp"`,
`"funnelsqp"`, and `"filterslp"`.

Inconsistent combinations are rejected: Uno's L-SR1 Hessian model is not
supported with the interior-point method (the `"ipopt"` or `"auto"` preset,
or an explicit `inequality_handling_method="interior_point"`); use an SQP
preset such as `preset="filtersqp"` with LSR1 instead. Unknown option names
are also rejected here, in Julia, against the option table of the pinned Uno
binary (Uno 2.9.3 via UnoSolver.jl 0.3.10), before any C call.

Derivatives use NonconvexCore's AD backend and dense, fixed coordinate patterns.
"""
struct UnoOptions{T<:NamedTuple}
    nt::T
    linear_constraints::Bool
end

function UnoOptions(;
    first_order::Bool=true,
    linear_constraints::Bool=false,
    preset="ipopt",
    hessian_model=nothing,
    kwargs...,
)
    hessian = isnothing(hessian_model) ? (first_order ? "LBFGS" : "exact") : hessian_model
    # Uno's presets choose the inequality handling method themselves, and some
    # explicit overrides of a preset are known to fail at solve time: the L-SR1
    # quasi-Newton model breaks the barrier line search of the interior-point
    # method (Uno's own CI exercises LSR1 only with the SQP `filtersqp` preset).
    ih = get(kwargs, :inequality_handling_method, nothing)
    interior_point = if ih !== nothing
        ih == "interior_point"
    else
        # `ipopt` uses the barrier method; `auto` may pick it at solve time.
        preset == "ipopt" || preset == "auto"
    end
    if hessian == "LSR1" && interior_point
        throw(
            ArgumentError(
                "Uno's L-SR1 Hessian model is inconsistent with the interior-point " *
                "method (the `ipopt` or `auto` preset, or " *
                "`inequality_handling_method=\"interior_point\"`): the barrier line " *
                "search fails. Use an SQP preset (`preset=\"filtersqp\"` or " *
                "`preset=\"funnelsqp\"`) with `hessian_model=\"LSR1\"`, or the exact " *
                "Hessian with the interior-point method.",
            ),
        )
    end
    nt = (; preset, hessian_model=hessian, kwargs...)
    for (key, value) in pairs(nt)
        value isa Union{AbstractString,Real} || throw(
            ArgumentError("Uno option $key must be a string, boolean, or real number")
        )
        # Validate the option names in Julia: Uno 2.9.3's C API lets
        # std::out_of_range escape its option setter on macOS arm64, aborting
        # the whole process for an unknown option. The list is the complete
        # option table of the pinned Uno binary (Uno_jll 2.9.3 through
        # UnoSolver.jl 0.3.10); update it when the compat bounds move.
        String(key) in _UNO_OPTION_NAMES || throw(
            ArgumentError(
                "`$key` is not an option of Uno v$(UnoSolver.uno_version()) " *
                "(UnoSolver.jl v$(pkgversion(UnoSolver))).",
            ),
        )
    end
    # UnoSolver's option dispatcher accepts Float64 and String, not all subtypes.
    normalized = map(nt) do value
        if value isa AbstractString
            String(value)
        elseif value isa Integer
            value
        else
            Float64(value)
        end
    end
    return UnoOptions(normalized, linear_constraints)
end

"""
    Workspace(model::NonconvexCore.VecModel, UnoAlg(), x0; options=UnoOptions(),
              lagrangian_hessian=nothing)

Construct a reusable `UnoWorkspace`. Convert other Nonconvex models with
`NonconvexCore.tovecmodel` first. `NonconvexCore.optimize!(workspace)` solves from
`workspace.x0`; `NonconvexCore.reset!(workspace, x0)` replaces that initial point.
Iteration callbacks are currently unsupported. A custom `lagrangian_hessian`
receives `(x, objective_multiplier, constraint_multipliers)` and returns a matrix.
"""
mutable struct UnoWorkspace{M,P,O,H} <: Workspace
    model::M
    problem::P
    x0::Vector{Float64}
    options::O
    lagrangian_hessian::H
end

"""
    UnoResult

Result of a Uno solve. Fields: `minimizer`, `minimum`, native `problem`, full Uno
execution statistics in `result`, symbolic `status`, and objective evaluation
count `fcalls`. On convergence `status == :first_order`; this denotes the
optimality condition, regardless of the Hessian mode used.
"""
struct UnoResult{X,P,R} <: AbstractResult
    minimizer::X
    minimum::Float64
    problem::P
    result::R
    status::Symbol
    fcalls::Int
end

function Workspace(
    model::VecModel,
    ::UnoAlg,
    x0::AbstractVector=getinit(model);
    options::UnoOptions=UnoOptions(),
    callback=NoCallback(),
    lagrangian_hessian=nothing,
)
    callback isa NoCallback || throw(
        ArgumentError(
            "UnoSolver.jl does not expose iteration callbacks; callback is unsupported"
        ),
    )
    if lagrangian_hessian !== nothing && options.nt.hessian_model != "exact"
        @warn "A custom `lagrangian_hessian` was provided but the Hessian model is " *
            "`$(options.nt.hessian_model)`, so it will never be called. Use " *
            "`first_order=false` (or `hessian_model=\"exact\"`) to enable it." maxlog = 1
    end
    x = Float64.(x0)
    length(x) == length(getmin(model)) || throw(DimensionMismatch("Invalid initial point"))
    problem = uno_problem(model, x, options, lagrangian_hessian)
    return UnoWorkspace(model, problem, copy(x), options, lagrangian_hessian)
end

function optimize!(workspace::UnoWorkspace)
    (; problem, x0, options) = workspace
    UnoSolver.uno_set_initial_primal_iterate(problem, x0)
    solver = uno_solver(options)
    GC.@preserve solver problem UnoSolver.uno_optimize(solver, problem)
    stats = UnoSolver.uno_statistics(problem, solver)
    return UnoResult(
        copy(stats.solution),
        stats.objective,
        problem,
        stats,
        stats.status,
        stats.number_objective_evaluations,
    )
end

# Use the typed native setters so invalid presets are checked as well, and the
# solver already has its finalizer if setting any option fails.
function uno_solver(options::UnoOptions)
    solver = UnoSolver.uno_solver()
    for (key, value) in pairs(options.nt)
        success = if key == :preset
            UnoSolver.uno_set_solver_preset(solver, value)
        elseif value isa String
            UnoSolver.uno_set_solver_string_option(solver, string(key), value)
        elseif value isa Bool
            UnoSolver.uno_set_solver_bool_option(solver, string(key), value)
        elseif value isa Integer
            UnoSolver.uno_set_solver_integer_option(solver, string(key), Cint(value))
        else
            UnoSolver.uno_set_solver_double_option(solver, string(key), value)
        end
        success || throw(ArgumentError("Uno rejected option $key=$(repr(value))"))
    end
    return solver
end

function reset!(workspace::UnoWorkspace, x0=nothing)
    if x0 !== nothing
        length(x0) == length(workspace.x0) ||
            throw(DimensionMismatch("Invalid initial point"))
        workspace.x0 .= x0
    end
    return workspace
end

# Differentiate each constraint separately, preserving custom derivative rules
# and avoiding AD through VectorOfFunctions' concatenation implementation.
function constraint_dot(constraints, x, multipliers)
    value = zero(eltype(x))
    offset = 0
    for c in constraints.fs
        dim = NonconvexCore.getdim(c)
        value += dot(c(x), view(multipliers, (offset + 1):(offset + dim)))
        offset += dim
    end
    return value
end

# The complete option table of the pinned Uno binary (Uno 2.9.3 through
# UnoSolver.jl 0.3.10), used to reject unknown option names in Julia before
# any C call: the Uno 2.9.3 C API lets std::out_of_range escape its option
# setter on macOS arm64 and abort the process.
const _UNO_OPTION_NAMES = Set((
    "primal_tolerance",
    "dual_tolerance",
    "tight_unscaled_stationarity_tolerance",
    "tight_unscaled_primal_infeasibility_tolerance",
    "tight_unscaled_complementarity_tolerance",
    "acceptable_tolerance",
    "acceptable_unscaled_stationarity_tolerance",
    "acceptable_unscaled_primal_infeasibility_tolerance",
    "acceptable_unscaled_complementarity_tolerance",
    "acceptable_objective_change_tolerance",
    "acceptable_iteration_threshold",
    "max_iterations",
    "time_limit",
    "print_solution",
    "diverging_iterate_threshold",
    "unbounded_objective_threshold",
    "logger",
    "constraint_relaxation_strategy",
    "inequality_handling_method",
    "globalization_mechanism",
    "globalization_strategy",
    "hessian_model",
    "inertia_correction_strategy",
    "progress_norm",
    "residual_norm",
    "residual_scaling_threshold",
    "protect_actual_reduction_against_roundoff",
    "protected_actual_reduction_macheps_coefficient",
    "print_subproblem",
    "use_function_scaling",
    "function_scaling_threshold",
    "function_scaling_lower_bound",
    "print_minor_iterations",
    "print_extended_statistics",
    "armijo_decrease_fraction",
    "armijo_tolerance",
    "switching_delta",
    "switching_merit_exponent",
    "switching_infeasibility_exponent",
    "sufficient_infeasibility_decrease_ratio",
    "filter_beta",
    "filter_gamma",
    "filter_ubd",
    "filter_fact",
    "filter_capacity",
    "filter_sufficient_infeasibility_decrease_factor",
    "filter_reset_iteration_threshold",
    "max_number_filter_resets",
    "filter_ubd_raise_trigger",
    "filter_ubd_raise_factor",
    "filter_ubd_max_raises",
    "filter_merit_max_increase",
    "filter_margin_uses_entry_infeasibility",
    "funnel_kappa",
    "funnel_beta",
    "funnel_gamma",
    "funnel_ubd",
    "funnel_fact",
    "funnel_update_strategy",
    "funnel_require_acceptance_wrt_current_iterate",
    "LS_backtracking_ratio",
    "LS_scale_duals_with_step_length",
    "gamma_alpha",
    "LS_watchdog_shortened_iterations_trigger",
    "LS_watchdog_max_trial_iterations",
    "quasi_newton_memory_size",
    "quasi_newton_delta_lower_bound",
    "quasi_newton_delta_upper_bound",
    "LBFGS_max_skips_before_reset",
    "LSR1_pivot_max_magnitude",
    "regularization_failure_threshold",
    "regularization_increase_factor",
    "primal_regularization_initial_factor",
    "dual_regularization_fraction",
    "primal_regularization_lb",
    "primal_regularization_decrease_factor",
    "primal_regularization_fast_increase_factor",
    "primal_regularization_slow_increase_factor",
    "threshold_unsuccessful_attempts",
    "regularize_all_variables",
    "iterative_refinement_min_steps",
    "iterative_refinement_max_steps",
    "iterative_refinement_residual_ratio_max",
    "iterative_refinement_residual_ratio_singular",
    "iterative_refinement_improvement_factor",
    "TR_radius",
    "TR_increase_factor",
    "TR_decrease_factor",
    "TR_aggressive_decrease_factor",
    "TR_activity_tolerance",
    "TR_min_radius",
    "TR_radius_reset_threshold",
    "switch_to_optimality_requires_linearized_feasibility",
    "l1_constraint_violation_coefficient",
    "use_proximal_term",
    "barrier_function",
    "barrier_initial_parameter",
    "barrier_default_multiplier",
    "barrier_tau_min",
    "barrier_k_sigma",
    "barrier_k_mu",
    "barrier_theta_mu",
    "barrier_k_epsilon",
    "barrier_update_fraction",
    "barrier_regularization_exponent",
    "barrier_small_direction_factor",
    "barrier_push_variable_to_interior_k1",
    "barrier_push_variable_to_interior_k2",
    "barrier_damping_factor",
    "barrier_small_infeasibility_factor",
    "least_square_multiplier_max_norm",
    "bound_multiplier_max_norm",
    "constraint_violation_tolerance",
    "bound_relaxation_factor",
    "BQPD_kmax_heuristic",
    "QP_solver",
    "LP_solver",
    "linear_solver",
    "preset",
    "option_file",
    "write_solution_to_file",
    "SOC_max_iterations",
    "SOC_infeasibility_fraction",
    "libhsl_path",
    "MA57_use_scaling",
    "dump_kkt_path",
    "dump_kkt_frequency",
))

# Keep the AD entry points separate so regression tests can forbid fallback AD
# when a converted JuMP model supplies all derivatives itself.
_value_gradient = (args...) -> value_gradient(args...)
_value_jacobian = (args...) -> value_jacobian(args...)
_value_hessian = (args...) -> value_hessian(args...)

function raw_function(f)
    f isa NonconvexCore.FlattenedFunction &&
        return raw_function(NonconvexCore.getrawfunction(f))
    f isa FunctionWrapper && return raw_function(NonconvexCore.getfunction(f))
    return f
end

function jump_evaluator(model)
    raw = raw_function(NonconvexCore.getfunction(getobjective(model)))
    raw isa NonconvexCore.JuMPModelObjective || return nothing
    evaluator = raw.m
    for (constraints, T, count) in (
        (model.ineq_constraints, NonconvexCore.JuMPModelIneq, evaluator.nineq),
        (model.eq_constraints, NonconvexCore.JuMPModelEq, evaluator.neq),
    )
        if count == 0
            isempty(constraints.fs) || return nothing
        else
            length(constraints.fs) == 1 || return nothing
            c = only(constraints.fs)
            f = raw_function(c.f)
            f isa T && f.m === evaluator && NonconvexCore.getdim(c) == count ||
                return nothing
            all(iszero, c.rhs) || return nothing
        end
    end
    return evaluator
end

function uno_problem(model, x0, options, lagrangian_hessian)
    objective = getobjective(model)
    objective === nothing && throw(ArgumentError("The model has no objective."))
    ineq, eq = model.ineq_constraints, model.eq_constraints
    ni, ne = getnineqconstraints(model), getneqconstraints(model)
    n, m = length(x0), ni + ne
    jump = jump_evaluator(model)
    # Include structural zeros: a derivative zero at x0 may become nonzero later.
    if jump === nothing
        jrows = Cint[i for j in 1:n for i in 1:m]
        jcols = Cint[j for j in 1:n for i in 1:m]
    else
        jr, jc = NonconvexCore.jump_jacobian_structure(jump)
        jrows, jcols = Cint.(jr), Cint.(jc)
    end
    exact = options.nt.hessian_model == "exact"
    if !exact
        hrows, hcols = Cint[], Cint[]
    elseif jump !== nothing && lagrangian_hessian === nothing
        hr, hc = NonconvexCore.jump_hessian_structure(jump)
        hrows, hcols = Cint.(hr), Cint.(hc)
    else
        hrows = Cint[i for j in 1:n for i in j:n]
        hcols = Cint[j for j in 1:n for i in j:n]
    end

    eval_objective = (_, x) -> Float64(objective(x))
    function eval_gradient(_, g, x)
        if jump !== nothing
            NonconvexCore.jump_objective_gradient!(jump, g, x)
            g .*= objective.multiple[]
            return nothing
        end
        grad = _value_gradient(objective, x)[2]
        grad === nothing ? fill!(g, 0.0) : copyto!(g, grad)
        return nothing
    end
    function eval_constraints(_, c, x)
        if jump !== nothing
            NonconvexCore.jump_constraint_values!(jump, c, x)
            return nothing
        end
        ni == 0 || (c[1:ni] .= ineq(x))
        ne == 0 || (c[(ni + 1):m] .= eq(x))
        return nothing
    end
    function jacobian(x)
        Ji = ni == 0 ? zeros(0, n) : _value_jacobian(ineq, x)[2]
        Je = ne == 0 ? zeros(0, n) : _value_jacobian(eq, x)[2]
        return vcat(Ji, Je)
    end
    J0 = if !options.linear_constraints
        nothing
    elseif jump === nothing
        vec(jacobian(x0))
    else
        values = zeros(length(jrows))
        NonconvexCore.jump_jacobian!(jump, values, x0)
        values
    end
    function eval_jacobian(_, values, x)
        if J0 !== nothing
            copyto!(values, J0)
        elseif jump !== nothing
            NonconvexCore.jump_jacobian!(jump, values, x)
        else
            copyto!(values, vec(jacobian(x)))
        end
        return nothing
    end
    function eval_hessian(_, values, x, multipliers, objective_multiplier)
        if jump !== nothing && lagrangian_hessian === nothing
            NonconvexCore.jump_hessian_lagrangian!(
                jump, values, x, objective_multiplier * objective.multiple[], multipliers
            )
            return nothing
        end
        H = if lagrangian_hessian !== nothing
            lagrangian_hessian(x, objective_multiplier, multipliers)
        else
            L =
                z ->
                    objective_multiplier * objective(z) +
                    constraint_dot(ineq, z, view(multipliers, 1:ni)) +
                    constraint_dot(eq, z, view(multipliers, (ni + 1):m))
            _value_hessian(L, x)
        end
        for k in eachindex(values)
            values[k] = H[hrows[k], hcols[k]]
        end
        return nothing
    end
    return UnoSolver.uno_model(
        "NLP",
        true,
        n,
        m,
        Float64.(getmin(model)),
        Float64.(getmax(model)),
        vcat(fill(-Inf, ni), zeros(ne)),
        zeros(m),
        jrows,
        jcols,
        length(jrows),
        hrows,
        hcols,
        length(hrows),
        eval_objective,
        eval_constraints,
        eval_gradient,
        eval_jacobian,
        exact ? eval_hessian : nothing,
        nothing,
        nothing,
        nothing,
        model,
        'L',
        1,
        copy(x0),
        zeros(m),
    )
end

end
