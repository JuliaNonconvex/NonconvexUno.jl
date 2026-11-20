@testset "Derivative coordinates retain initially zero entries" begin
    model = Model(x -> x[1]^4 + x[2]^2)
    addvar!(model, [-Inf, -Inf], [Inf, Inf])
    add_ineq_constraint!(model, x -> x[1] * x[2], 2.0)
    vm, x0, _ = NonconvexCore.tovecmodel(model, zeros(2))
    w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(first_order=false))
    p = w.problem
    x = [2.0, 3.0]
    g, c, j, h = zeros(2), zeros(1), zeros(2), zeros(3)
    p.eval_gradient(p.user_model, g, x)
    p.eval_constraints(p.user_model, c, x)
    p.eval_jacobian(p.user_model, j, x)
    p.eval_hessian(p.user_model, h, x, [5.0], 2.0)
    @test g ≈ [32.0, 6.0]
    @test c ≈ [4.0]
    @test j ≈ [3.0, 2.0]
    @test h ≈ [96.0, 5.0, 4.0]
end

@testset "JuMP linear Jacobian caching" begin
    jm = JuMP.Model()
    JuMP.@variable(jm, x[1:2], start=0.5)
    JuMP.@objective(jm, Min, sum(x[i]^2 for i in 1:2) + 0.1 * sin(x[1]))
    JuMP.@constraint(jm, x[1] + x[2] == 2)
    dm = NonconvexCore.DictModel(jm)
    r = optimize(dm, UnoAlg(); options=UnoOptions(linear_constraints=true, logger="SILENT"))
    # min x1^2 + x2^2 + 0.1*sin(x1) s.t. x1 + x2 == 2: optimum near [0.9865, 1.0135]
    @test r.minimizer[:x] ≈ [0.9865, 1.0135] atol = 1e-2
    @test r.minimum ≈ 2.0837 atol = 1e-3
end

@testset "Constant objective and initial-point defaults" begin
    model = Model(x -> 3.0)
    addvar!(model, [-1.0], [1.0]; init=[0.0])
    vm, _, _ = NonconvexCore.tovecmodel(model)
    w = Workspace(vm, UnoAlg(); options=UnoOptions(logger="SILENT"))
    g = ones(1)
    w.problem.eval_gradient(w.problem.user_model, g, [0.0])
    @test g == [0.0]
    @test reset!(w) === w
    r = optimize!(w)
    @test r.minimum == 3.0
    @test r.status == :first_order
    @test_throws DimensionMismatch Workspace(vm, UnoAlg(), zeros(2))
    empty_vm = NonconvexCore.VecModel(
        nothing,
        vm.eq_constraints,
        vm.ineq_constraints,
        vm.sd_constraints,
        vm.box_min,
        vm.box_max,
        vm.init,
        vm.integer,
    )
    @test_throws ArgumentError Workspace(empty_vm, UnoAlg(), [0.0])
end

@testset "Structured variables" begin
    model = DictModel(d -> sum(abs2, d[:x] .- [1.0, 2.0]))
    addvar!(model, :x, [-Inf, -Inf], [Inf, Inf]; init=[0.0, 0.0])
    r = optimize(model, UnoAlg(); options=UnoOptions(logger="SILENT"))
    @test r.minimizer[:x] ≈ [1.0, 2.0] atol = 1e-5
end

@testset "SLP and funnel presets" begin
    for (preset, hessian_model) in (("filterslp", "zero"), ("funnelsqp", "exact"))
        model = Model(x -> -sum(x))
        addvar!(model, zeros(2), ones(2))
        r = optimize(
            model,
            UnoAlg(),
            [0.5, 0.5];
            options=UnoOptions(; preset, hessian_model, logger="SILENT"),
        )
        @test r.minimizer ≈ ones(2) atol = 1e-5
        @test r.minimum ≈ -2.0 atol = 1e-5
    end
end

@testset "JuMP objective-only maximization and scaling" begin
    jm = JuMP.Model()
    JuMP.@variable(jm, x, start=0.5)
    JuMP.@objective(jm, Max, -(x - 2)^2)
    dm = NonconvexCore.DictModel(jm)
    vm, x0, _ = NonconvexCore.tovecmodel(dm, NonconvexCore.getinit(dm))
    NonconvexCore.getobjective(vm).multiple[] = 3.0
    w = Workspace(vm, UnoAlg(), x0; options=UnoOptions(first_order=false, logger="SILENT"))
    p = w.problem
    g, h = zeros(1), zeros(1)
    p.eval_gradient(p.user_model, g, [0.0])
    p.eval_hessian(p.user_model, h, [0.0], Float64[], 2.0)
    @test g ≈ [-12.0]
    @test h ≈ [12.0]
    @test optimize!(w).minimizer ≈ [2.0] atol = 1e-5
    calls = Ref(0)
    hessian = (x, sigma, lambda) -> begin
        calls[] += 1
        reshape([6sigma], 1, 1)
    end
    w = Workspace(
        vm,
        UnoAlg(),
        x0;
        options=UnoOptions(first_order=false, logger="SILENT"),
        lagrangian_hessian=hessian,
    )
    @test optimize!(w).minimizer ≈ [2.0] atol = 1e-5
    @test calls[] > 0
    w = Workspace(
        vm,
        UnoAlg(),
        x0;
        options=UnoOptions(logger="SILENT"),
        lagrangian_hessian=(args...) -> error("First-order mode requested a Hessian"),
    )
    @test optimize!(w).minimizer ≈ [2.0] atol = 1e-5
end

@testset "Unsupported variable types are rejected" begin
    model = Model(x -> sum(x))
    addvar!(model, [0.0], [2.0]; integer=[true])
    @test_throws ArgumentError optimize(model, UnoAlg(), [1.0])
end

@testset "Modified JuMP models use all constraints" begin
    jm = JuMP.Model()
    JuMP.@variable(jm, x, start=0.5)
    JuMP.@objective(jm, Min, x^2)
    dm = NonconvexCore.DictModel(jm)
    add_eq_constraint!(dm, d -> d[:x], 1.0)
    vm, _, _ = NonconvexCore.tovecmodel(dm, NonconvexCore.getinit(dm))
    @test NonconvexUno.jump_evaluator(vm) === nothing
    r = optimize(dm, UnoAlg(); options=UnoOptions(logger="SILENT"))
    @test r.minimizer[:x] ≈ 1.0 atol = 1e-5
    @test r.minimum ≈ 1.0 atol = 1e-5
    # A nonlinear JuMP objective with no JuMP constraint rows, plus a
    # Nonconvex inequality appended after conversion: the evaluator's empty
    # constraint channel rejects the fast path instead of dropping the row.
    jm2 = JuMP.Model()
    JuMP.@variable(jm2, x, start=0.5)
    JuMP.@objective(jm2, Min, x^2 + 0.1 * sin(x))
    dm2 = NonconvexCore.DictModel(jm2)
    add_ineq_constraint!(dm2, d -> 2 - d[:x])
    vm2, _, _ = NonconvexCore.tovecmodel(dm2, NonconvexCore.getinit(dm2))
    @test NonconvexUno.jump_evaluator(vm2) === nothing
    r2 = optimize(dm2, UnoAlg(); options=UnoOptions(logger="SILENT"))
    @test r2.minimizer[:x] ≈ 2.0 atol = 1e-4
    @test r2.minimum ≈ 4.0 + 0.1 * sin(2.0) atol = 1e-4
end
