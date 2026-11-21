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

run("a: create the Uno model only") do
    vm, x0, _ = NonconvexCore.tovecmodel(plain_model(), [2.0])
    NonconvexUno.uno_problem(vm, x0, UnoOptions(logger="SILENT"), nothing)
end

run("b: create the solver only") do
    UnoSolver.uno_solver()
end

run("c: solver + preset") do
    s = UnoSolver.uno_solver()
    UnoSolver.uno_set_solver_preset(s, "ipopt")
end

run("d: solver + preset + hessian_model") do
    s = UnoSolver.uno_solver()
    UnoSolver.uno_set_solver_preset(s, "ipopt")
    UnoSolver.uno_set_solver_string_option(s, "hessian_model", "LBFGS")
end

run("e: solver + preset + hessian_model + bad bool setter (caught, no Julia throw)") do
    s = UnoSolver.uno_solver()
    UnoSolver.uno_set_solver_preset(s, "ipopt")
    UnoSolver.uno_set_solver_string_option(s, "hessian_model", "LBFGS")
    UnoSolver.uno_set_solver_bool_option(s, "not_an_uno_option", true)
end

run("f: full failing path, Julia throw caught") do
    try
        optimize(
            plain_model(),
            UnoAlg(),
            [2.0];
            options=UnoOptions(not_an_uno_option=true, logger="SILENT"),
        )
    catch e
        e isa ArgumentError || rethrow()
    end
end

run("g: full failing path with an unknown string option") do
    try
        optimize(
            plain_model(),
            UnoAlg(),
            [2.0];
            options=UnoOptions(not_an_uno_option="x", logger="SILENT"),
        )
    catch e
        e isa ArgumentError || rethrow()
    end
end

run("h: full failing path with an unknown double option") do
    try
        optimize(
            plain_model(),
            UnoAlg(),
            [2.0];
            options=UnoOptions(not_an_uno_option=1.5, logger="SILENT"),
        )
    catch e
        e isa ArgumentError || rethrow()
    end
end

println("=== ALL FINE-GRAINED STEPS DONE")
flush(stdout)
