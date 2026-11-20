using NonconvexUno
using NonconvexCore: optimize!, reset!
import UnoSolver
import JuMP

# The do-block passes the function as the first argument: `run("label") do ... end`
# calls `run(thunk, "label")`.
function run(thunk::Function, label::String)
    println("=== ", label)
    flush(stdout)
    try
        r = thunk()
        println("=== ", label, " -> ", r isa NonconvexUno.UnoResult ? r.status : "ok")
    catch e
        println("=== ", label, " -> Julia exception: ", sprint(showerror, e))
    end
    flush(stdout)
end

# Replicate the first suite testset exactly: both first-order modes, four
# constraint layouts, and the same option tolerances.
f(x::AbstractVector) = sqrt(x[2])
g(x::AbstractVector, a, b) = (a * x[1] + b)^3 - x[2]
for first_order in (true, false)
    run("suite 1: simple constraints, first_order=$first_order") do
        m = Model(f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        optimize(m, UnoAlg(), [1.234, 2.345]; options=UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8,
        ))
    end
    run("suite 1: equality constraints, first_order=$first_order") do
        m = Model(f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        add_eq_constraint!(m, x -> sum(x) - 1 / 3 - 8 / 27)
        optimize(m, UnoAlg(), [1.234, 2.345]; options=UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8,
        ))
    end
    run("suite 1: block constraints, first_order=$first_order") do
        m = Model(f)
        addvar!(m, [1e-4, 1e-4], [10.0, 10.0])
        add_ineq_constraint!(m, FunctionWrapper(x -> [g(x, 2, 0), g(x, -1, 1)], 2))
        optimize(m, UnoAlg(), [1.234, 2.345]; options=UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8,
        ))
    end
    run("suite 1: infinite upper bound, first_order=$first_order") do
        m = Model(f)
        addvar!(m, [1e-4, 1e-4], [Inf, Inf])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        optimize(m, UnoAlg(), [1.234, 2.345]; options=UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8,
        ))
    end
    run("suite 1: infinite lower bound, first_order=$first_order") do
        m = Model(f)
        addvar!(m, [-Inf, 1e-4], [10, 10])
        add_ineq_constraint!(m, x -> g(x, 2, 0))
        add_ineq_constraint!(m, x -> g(x, -1, 1))
        optimize(m, UnoAlg(), [1.234, 2.345]; options=UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8,
        ))
    end
end

run("nonlinear inequalities first_order=true") do
    model = Model(x -> sqrt(x[2]))
    addvar!(model, [1e-4, 1e-4], [Inf, Inf])
    add_ineq_constraint!(model, FunctionWrapper(x -> [(2x[1])^3 - x[2], (1-x[1])^3 - x[2]], 2))
    optimize(model, UnoAlg(), [1.234, 2.345]; options=UnoOptions(logger="SILENT"))
end
run("nonlinear inequalities first_order=false") do
    model = Model(x -> sqrt(x[2]))
    addvar!(model, [1e-4, 1e-4], [Inf, Inf])
    add_ineq_constraint!(model, FunctionWrapper(x -> [(2x[1])^3 - x[2], (1-x[1])^3 - x[2]], 2))
    optimize(model, UnoAlg(), [1.234, 2.345]; options=UnoOptions(first_order=false, logger="SILENT"))
end
run("equality constraints") do
    f(x) = sqrt(x[2])
    g(x, a, b) = (a * x[1] + b)^3 - x[2]
    model = Model(f)
    addvar!(model, [1e-4, 1e-4], [Inf, Inf])
    add_ineq_constraint!(model, x -> g(x, 2, 0))
    add_ineq_constraint!(model, x -> g(x, -1, 1))
    add_eq_constraint!(model, x -> sum(x) - 1 / 3 - 8 / 27)
    optimize(model, UnoAlg(), [1.234, 2.345]; options=UnoOptions(first_order=false, logger="SILENT"))
end
run("mixed ipopt exact linear_constraints") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf, -Inf], [Inf, Inf])
    add_eq_constraint!(model, x -> sum(x), 3.0)
    add_ineq_constraint!(model, x -> 2 - x[1])
    optimize(model, UnoAlg(), [2.5, 0.5]; options=UnoOptions(;
        first_order=false, linear_constraints=true, preset="ipopt", logger="SILENT",
    ))
end
run("mixed filtersqp exact linear_constraints") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf, -Inf], [Inf, Inf])
    add_eq_constraint!(model, x -> sum(x), 3.0)
    add_ineq_constraint!(model, x -> 2 - x[1])
    optimize(model, UnoAlg(), [2.5, 0.5]; options=UnoOptions(;
        first_order=false, linear_constraints=true, preset="filtersqp", logger="SILENT",
    ))
end
run("workspace reuse and reset") do
    model = Model(x -> sum(abs2, x .- 2))
    addvar!(model, [-Inf, -Inf], [Inf, Inf])
    vm, x0, _ = NonconvexCore.tovecmodel(model, [0.5, 0.5])
    w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(logger="SILENT"))
    optimize!(w)
    reset!(w, [0.25, 0.75])
    optimize!(w)
end
run("nonlinear equality exact and custom hessian") do
    model = Model(x -> -x[1])
    addvar!(model, [-Inf, -Inf], [Inf, Inf])
    add_eq_constraint!(model, x -> sum(abs2, x), 1.0)
    hessian = (x, sigma, lambda) -> 2lambda[1] * Matrix{Float64}(I, 2, 2)
    optimize(model, UnoAlg(), [0.8, 0.6];
        options=UnoOptions(first_order=false, logger="SILENT"), lagrangian_hessian=hessian)
end
run("max_iterations=0") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf], [Inf])
    optimize(model, UnoAlg(), [2.0]; options=UnoOptions(max_iterations=0, logger="SILENT"))
end
run("filterslp zero") do
    model = Model(x -> -sum(x))
    addvar!(model, zeros(2), ones(2))
    optimize(model, UnoAlg(), [0.5, 0.5]; options=UnoOptions(; preset="filterslp", hessian_model="zero", logger="SILENT"))
end
run("funnelsqp exact") do
    model = Model(x -> -sum(x))
    addvar!(model, zeros(2), ones(2))
    optimize(model, UnoAlg(), [0.5, 0.5]; options=UnoOptions(; preset="funnelsqp", hessian_model="exact", logger="SILENT"))
end
run("filtersqp LSR1") do
    model = Model(x -> -sum(x))
    addvar!(model, zeros(2), ones(2))
    optimize(model, UnoAlg(), [0.5, 0.5]; options=UnoOptions(; preset="filtersqp", hessian_model="LSR1", logger="SILENT"))
end
run("JuMP evaluator first_order=true") do
    jm = JuMP.Model()
    JuMP.@variable(jm, x[1:2], start=0.5)
    JuMP.@objective(jm, Min, (x[1]-2)^2 + (x[2]-1)^2 + 0.1 * sin(x[1])^2)
    JuMP.@constraint(jm, x[1] + x[2] == 2)
    JuMP.@constraint(jm, x[1]^2 + x[2]^2 <= 4)
    JuMP.@constraint(jm, exp(x[2]) >= 1.1)
    model = NonconvexCore.DictModel(jm)
    optimize(model, UnoAlg(); options=UnoOptions(logger="SILENT"))
end
run("JuMP evaluator first_order=false") do
    jm = JuMP.Model()
    JuMP.@variable(jm, x[1:2], start=0.5)
    JuMP.@objective(jm, Min, (x[1]-2)^2 + (x[2]-1)^2 + 0.1 * sin(x[1])^2)
    JuMP.@constraint(jm, x[1] + x[2] == 2)
    JuMP.@constraint(jm, x[1]^2 + x[2]^2 <= 4)
    JuMP.@constraint(jm, exp(x[2]) >= 1.1)
    model = NonconvexCore.DictModel(jm)
    optimize(model, UnoAlg(); options=UnoOptions(first_order=false, logger="SILENT"))
end
run("JuMP objective-only maximization and scaling") do
    jm = JuMP.Model()
    JuMP.@variable(jm, x, start=0.5)
    JuMP.@objective(jm, Max, -(x-2)^2)
    dm = NonconvexCore.DictModel(jm)
    vm, x0, _ = NonconvexCore.tovecmodel(dm, NonconvexCore.getinit(dm))
    NonconvexCore.getobjective(vm).multiple[] = 3.0
    w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(first_order=false, logger="SILENT"))
    optimize!(w)
end
println("=== ALL SCENARIOS DONE")
flush(stdout)
