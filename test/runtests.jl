using Xerxes
using Test

@testset "Xerxes.jl" begin
    prog = Xerxes.Program(
        [
            1.0 -1 1 0 0
            2.0 -1 0 1 0
            0    1 0 0 1
        ],
        [1.0, 3, 5],
        [4.0, 3, 0, 0, 0],
        [1, 2],
        [3, 4, 5],
        [0.0, 0, 1, 3, 5],
        [-4.0, -3, 0, 0, 0]
    )

    status = Xerxes.primal_iterate!(prog)
    @test status == Xerxes.Nonoptimal
    @test issetequal(prog.basic_vars, [1, 4, 5])
    @test issetequal(prog.nonbasic_vars, [3, 2])
    @test prog.A[:, sort(prog.basic_vars)] ==
        [ 1 0 0
          2 1 0
          0 0 1 ]
    @test prog.A[:, sort(prog.nonbasic_vars)] ==
        [ -1 1
          -1 0
           1 0 ]
    @test prog.x[sort(prog.basic_vars)] == [1, 1, 5]
    @test prog.z[sort(prog.nonbasic_vars)] == [-7, 4]

    status = Xerxes.primal_iterate!(prog)
    @test status == Xerxes.Nonoptimal
    @test issetequal(prog.basic_vars, [1, 2, 5])
    @test issetequal(prog.nonbasic_vars, [3, 4])
    @test prog.A[:, sort(prog.basic_vars)] ==
        [ 1 -1 0
          2 -1 0
          0  1 1 ]
    @test prog.A[:, sort(prog.nonbasic_vars)] ==
        [ 1 0
          0 1
          0 0 ]
    @test prog.x[sort(prog.basic_vars)] == [2, 1, 4]
    @test prog.z[sort(prog.nonbasic_vars)] == [-10, 7]

    status = Xerxes.primal_iterate!(prog)
    @test status == Xerxes.Nonoptimal
    @test issetequal(prog.basic_vars, [1, 2, 3])
    @test issetequal(prog.nonbasic_vars, [5, 4])
    @test prog.A[:, sort(prog.basic_vars)] ==
        [ 1 -1 1
          2 -1 0
          0  1 0 ]
    @test prog.A[:, sort(prog.nonbasic_vars)] ==
        [ 0 0
          1 0
          0 1 ]
    @test prog.x[sort(prog.basic_vars)] == [4, 5, 2]
    @test prog.z[sort(prog.nonbasic_vars)] == [2, 5]

    status = Xerxes.primal_iterate!(prog)
    @test status == Xerxes.Optimal
end
