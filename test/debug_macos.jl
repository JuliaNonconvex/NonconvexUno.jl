using NonconvexUno
import UnoSolver
import JuMP

println("uno version: ", UnoSolver.uno_version())
flush(stdout)
if isdefined(UnoSolver, :HSL_jll)
    nothing
end
try
    import HSL_jll
    println("HSL functional: ", @ccall(HSL_jll.libhsl.LIBHSL_isfunctional()::Bool))
catch e
    println("HSL check failed: ", e)
end
flush(stdout)

solver = UnoSolver.uno_solver()
for name in ("linear_solver", "QP_solver", "LP_solver")
    try
        println("$name = ", UnoSolver.uno_get_solver_string_option(solver, name))
    catch e
        println("$name getter failed: ", e)
    end
    flush(stdout)
end

function sqrt_problem()
    model = Model(x -> sqrt(x[2]))
    addvar!(model, [1e-4, 1e-4], [Inf, Inf])
    add_ineq_constraint!(model, FunctionWrapper(x -> [(2x[1])^3 - x[2], (1-x[1])^3 - x[2]], 2))
    return model
end

for (label, options) in (
    ("default first_order=true", UnoOptions(logger="INFO")),
    ("first_order=false", UnoOptions(first_order=false, logger="INFO")),
    ("MUMPS", UnoOptions(linear_solver="MUMPS", logger="INFO")),
    ("SSIDS", UnoOptions(linear_solver="SSIDS", logger="INFO")),
)
    println("=== solving: $label")
    flush(stdout)
    try
        r = optimize(sqrt_problem(), UnoAlg(), [1.234, 2.345]; options)
        println("=== status: ", r.status, " min: ", r.minimum)
    catch e
        println("=== Julia exception: ", sprint(showerror, e))
    end
    flush(stdout)
end

println("=== direct JuMP + UnoSolver (upstream path)")
flush(stdout)
jm = JuMP.Model(() -> UnoSolver.Optimizer(preset="ipopt", logger="INFO"))
JuMP.@variable(jm, x[1:2] >= 1e-4, start=1.234)
JuMP.set_start_value(x[2], 2.345)
JuMP.@objective(jm, Min, sqrt(x[2]))
JuMP.@constraint(jm, (2x[1])^3 - x[2] <= 0)
JuMP.@constraint(jm, (1-x[1])^3 - x[2] <= 0)
try
    JuMP.optimize!(jm)
    println("=== JuMP termination: ", JuMP.termination_status(jm))
catch e
    println("=== JuMP Julia exception: ", sprint(showerror, e))
end
flush(stdout)
