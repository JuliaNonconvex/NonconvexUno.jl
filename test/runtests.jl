using NonconvexUno
using NonconvexCore: optimize!, reset!
using Test
using LinearAlgebra
using UnoSolver: UnoSolver
using JuMP: JuMP

@testset "NonconvexUno" begin
    @testset "Nonlinear inequalities, first_order=$first_order" for first_order in
                                                                    (true, false)
        model = Model(x -> sqrt(x[2]))
        addvar!(model, [1e-4, 1e-4], [Inf, Inf])
        add_ineq_constraint!(
            model, FunctionWrapper(x -> [(2x[1])^3 - x[2], (1 - x[1])^3 - x[2]], 2)
        )
        options = UnoOptions(;
            first_order, logger="SILENT", primal_tolerance=1e-8, dual_tolerance=1e-8
        )
        r = optimize(model, UnoAlg(), [1.234, 2.345]; options)
        @test r.minimizer ≈ [1 / 3, 8 / 27] atol = 1e-5
        @test r.minimum ≈ sqrt(8 / 27) atol = 1e-6
        @test r.status == :first_order
        @test r.fcalls > 0
        @test if first_order
            r.result.number_hessian_evaluations == 0
        else
            r.result.number_hessian_evaluations > 0
        end
    end

    @testset "Mixed constraints and presets" for preset in ("ipopt", "filtersqp")
        model = Model(x -> sum(abs2, x))
        addvar!(model, [-Inf, -Inf], [Inf, Inf])
        add_eq_constraint!(model, x -> sum(x), 3.0)
        add_ineq_constraint!(model, x -> 2 - x[1])
        r = optimize(
            model,
            UnoAlg(),
            [2.5, 0.5];
            options=UnoOptions(;
                first_order=false, linear_constraints=true, preset, logger="SILENT"
            ),
        )
        @test r.minimizer ≈ [2.0, 1.0] atol = 1e-5
        @test r.minimum ≈ 5.0 atol = 1e-5
    end

    @testset "Unconstrained, bounds, and workspace reuse" begin
        for bounds in (([-Inf, -Inf], [Inf, Inf]), ([0.0, 0.0], [1.0, 1.0]))
            model = Model(x -> sum(abs2, x .- 2))
            addvar!(model, bounds...)
            vm, x0, _ = NonconvexCore.tovecmodel(model, [0.5, 0.5])
            w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(logger="SILENT"))
            r = optimize!(w)
            target = isinf(bounds[2][1]) ? [2.0, 2.0] : [1.0, 1.0]
            @test r.minimizer ≈ target atol = 1e-5
            saved = copy(r.minimizer)
            reset!(w, [0.25, 0.75])
            @test optimize!(w).minimizer ≈ target atol = 1e-5
            @test r.minimizer == saved
            @test_throws DimensionMismatch reset!(w, [0.0])
        end
    end

    @testset "Nonlinear equality sign, custom=$custom" for custom in (false, true)
        model = Model(x -> -x[1])
        addvar!(model, [-Inf, -Inf], [Inf, Inf])
        add_eq_constraint!(model, x -> sum(abs2, x), 1.0)
        calls = Ref(0)
        hessian = (x, sigma, lambda) -> begin
            calls[] += 1
            2lambda[1] * Matrix{Float64}(I, 2, 2)
        end
        r = optimize(
            model,
            UnoAlg(),
            [0.8, 0.6];
            options=UnoOptions(first_order=false, logger="SILENT"),
            lagrangian_hessian=custom ? hessian : nothing,
        )
        @test r.minimizer ≈ [1.0, 0.0] atol = 1e-5
        @test custom ? calls[] > 0 : calls[] == 0
    end

    @testset "Options and termination" begin
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
            @test getter(solver, string(key)) == value
        end
        @test UnoOptions().nt.hessian_model == "LBFGS"
        @test UnoOptions(first_order=false).nt.hessian_model == "exact"
        @test UnoOptions(hessian_model="identity").nt.hessian_model == "identity"
        @test UnoOptions(primal_tolerance=1.0f-6).nt.primal_tolerance isa Float64
        @test_throws ArgumentError UnoOptions(logger=[:SILENT])
        # Inconsistent combinations: L-SR1 breaks the interior-point barrier.
        @test_throws ArgumentError UnoOptions(hessian_model="LSR1")
        @test_throws ArgumentError UnoOptions(preset="auto", hessian_model="LSR1")
        @test_throws ArgumentError UnoOptions(
            preset="filtersqp",
            hessian_model="LSR1",
            inequality_handling_method="interior_point",
        )
        # LSR1 is fine with SQP presets or an explicit SQP method override.
        @test UnoOptions(preset="filtersqp", hessian_model="LSR1").nt.hessian_model ==
            "LSR1"
        @test UnoOptions(preset="funnelsqp", hessian_model="LSR1").nt.hessian_model ==
            "LSR1"
        @test UnoOptions(
            preset="ipopt",
            hessian_model="LSR1",
            inequality_handling_method="inequality_constrained",
        ).nt.hessian_model == "LSR1"
        model = Model(x -> sum(abs2, x))
        addvar!(model, [-Inf], [Inf])
        @test_throws ArgumentError optimize(model, UnoAlg(), [2.0]; callback=identity)
        @test_throws ArgumentError optimize(
            model,
            UnoAlg(),
            [2.0];
            options=UnoOptions(not_an_uno_option=true, logger="SILENT"),
        )
        r = optimize(
            model, UnoAlg(), [2.0]; options=UnoOptions(max_iterations=0, logger="SILENT")
        )
        @test r.status == :max_iter
    end

    @testset "JuMP evaluator derivatives, first_order=$first_order" for first_order in
                                                                        (true, false)
        jm = JuMP.Model()
        JuMP.@variable(jm, x[1:2], start=0.5)
        JuMP.@objective(jm, Min, (x[1] - 2)^2 + (x[2] - 1)^2 + 0.1 * sin(x[1])^2)
        JuMP.@constraint(jm, x[1] + x[2] == 2)
        JuMP.@constraint(jm, x[1]^2 + x[2]^2 <= 4)
        JuMP.@constraint(jm, exp(x[2]) >= 1.1)
        model = NonconvexCore.DictModel(jm)
        vm, x0, _ = NonconvexCore.tovecmodel(model, NonconvexCore.getinit(model))
        @test NonconvexUno.jump_evaluator(vm) !== nothing
        original = (
            NonconvexUno._value_gradient,
            NonconvexUno._value_jacobian,
            NonconvexUno._value_hessian,
        )
        forbidden = (args...) -> error("JuMP model unexpectedly used fallback AD")
        NonconvexUno._value_gradient = forbidden
        NonconvexUno._value_jacobian = forbidden
        NonconvexUno._value_hessian = forbidden
        try
            r = optimize(
                model, UnoAlg(); options=UnoOptions(; first_order, logger="SILENT")
            )
            @test r.status == :first_order
            @test r.result.number_objective_gradient_evaluations > 0
            @test r.result.number_jacobian_evaluations > 0
            @test if first_order
                r.result.number_hessian_evaluations == 0
            else
                r.result.number_hessian_evaluations > 0
            end
            JuMP.set_optimizer(jm, UnoSolver.Optimizer)
            JuMP.set_attribute(jm, "preset", "ipopt")
            JuMP.set_attribute(jm, "logger", "SILENT")
            JuMP.optimize!(jm)
            @test JuMP.is_solved_and_feasible(jm)
            @test r.minimum ≈ JuMP.objective_value(jm) atol = 1e-5
        finally
            NonconvexUno._value_gradient,
            NonconvexUno._value_jacobian,
            NonconvexUno._value_hessian = original
        end
    end
end

include("edge_cases.jl")

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(NonconvexUno)
end

# The JET analysis re-evaluates the package definitions in the live session,
# which invalidates the keyword constructors of the package types for anything
# that runs afterwards: it must be the last testset.
@testset "JET" begin
    using JET
    JET.test_package(NonconvexUno; target_modules=(NonconvexUno,))
end
