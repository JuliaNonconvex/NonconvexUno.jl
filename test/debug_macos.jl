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
    GC.gc(true)
    flush(stdout)
end

function plain_model()
    model = Model(x -> sum(abs2, x))
    addvar!(model, [-Inf], [Inf])
    return model
end

run("a: unknown option name rejected in Julia (no C call)") do
    try
        UnoOptions(not_an_uno_option=true, logger="SILENT")
        error("expected an ArgumentError")
    catch e
        e isa ArgumentError || rethrow()
    end
end

run("b: invalid string VALUE (valid name) on the C++ side") do
    s = UnoSolver.uno_solver()
    UnoSolver.uno_set_solver_string_option(s, "logger", "NOT_A_LEVEL")
end

run("c: invalid integer VALUE on the C++ side") do
    s = UnoSolver.uno_solver()
    UnoSolver.uno_set_solver_integer_option(s, "max_iterations", -1)
end

run("d: a normal solve still works") do
    optimize(
        plain_model(),
        UnoAlg(),
        [2.0];
        options=UnoOptions(max_iterations=200, logger="SILENT"),
    )
end

run("e: full suite options flow") do
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

println("=== ALL FINE-GRAINED STEPS DONE")
flush(stdout)
