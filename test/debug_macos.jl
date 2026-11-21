using NonconvexUno
using NonconvexCore: optimize!, reset!
import UnoSolver
import JuMP

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

run("options: solver creation with explicit options") do
    opts = UnoOptions(
        logger="SILENT",
        print_solution=false,
        max_iterations=42,
        primal_tolerance=1e-7,
        quasi_newton_memory_size=8,
        linear_solver="MUMPS",
    )
    NonconvexUno.uno_solver(opts)
end

run("options: getter round-trip") do
    opts = UnoOptions(
        logger="SILENT",
        print_solution=false,
        max_iterations=42,
        primal_tolerance=1e-7,
        quasi_newton_memory_size=8,
        linear_solver="MUMPS",
    )
    solver = NonconvexUno.uno_solver(opts)
    for (key, value) in pairs(opts.nt)
        key == :preset && continue
        getter = if value isa String
            UnoSolver.uno_get_solver_string_option
        elseif value isa Bool
            UnoSolver.uno_get_solver_bool_option
        elseif value isa Integer
            UnoSolver.uno_get_solver_integer_option
        else
            UnoSolver.uno_get_solver_double_option
        end
        getter(solver, string(key))
    end
end

run("options: constructor-only checks") do
    UnoOptions(logger="SILENT")
    UnoOptions()
    UnoOptions(first_order=false)
    UnoOptions(hessian_model="identity")
    UnoOptions(primal_tolerance=1.0f-6)
end

run("options: constructor throws") do
    for f in (
        () -> UnoOptions(logger=[:SILENT]),
        () -> UnoOptions(hessian_model="LSR1"),
        () -> UnoOptions(preset="auto", hessian_model="LSR1"),
        () -> UnoOptions(preset="filtersqp", hessian_model="LSR1", inequality_handling_method="interior_point"),
    )
        try
            f()
            error("expected an ArgumentError")
        catch e
            e isa ArgumentError || rethrow()
        end
    end
end

run("options: callback=identity throws") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf], [Inf])
    try
        optimize(model, UnoAlg(), [2.0]; callback=identity)
        error("expected an ArgumentError")
    catch e
        e isa ArgumentError || rethrow()
    end
end

run("options: not_an_uno_option throws") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf], [Inf])
    try
        optimize(
            model,
            UnoAlg(),
            [2.0];
            options=UnoOptions(not_an_uno_option=true, logger="SILENT"),
        )
        error("expected an ArgumentError")
    catch e
        e isa ArgumentError || rethrow()
    end
end

run("options: max_iterations=0 solve") do
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf], [Inf])
    optimize(model, UnoAlg(), [2.0]; options=UnoOptions(max_iterations=0, logger="SILENT"))
end

println("=== ALL OPTION STEPS DONE")
flush(stdout)
