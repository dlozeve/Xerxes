module Xerxes

using Logging

mutable struct Program
    A::Matrix
    b::Vector
    c::Vector
    nonbasic_vars::Vector
    basic_vars::Vector
    x::Vector # current primal solution
    z::Vector # current dual solution
end

@enum Status Optimal Nonoptimal Unbounded Infeasible

function primal_iterate!(prog::Program)
    B = view(prog.A, :, prog.basic_vars)
    N = view(prog.A, :, prog.nonbasic_vars)
    @debug "Submatrices" B N

    # Check for optimality
    if all(prog.z[prog.nonbasic_vars] .>= 0)
        @info "Program is optimal"
        return Optimal
    end
    @debug "z has negative values, therefore the program is not optimal"

    # Select entering variable
    j = argmin(prog.z)
    if prog.z[j] >= 0
        @debug "No entering variable"
        return Infeasible
    end
    @debug "Entering variable" j

    # Compute primal step direction
    ej = zeros(length(prog.x))
    ej[j] += 1
    primal_step = B \ (N * ej[prog.nonbasic_vars])
    @debug "Primal step direction" primal_step

    # Compute primal step length
    inv_primal_step_length = -Inf
    i = 0
    for k in eachindex(prog.basic_vars)
        if primal_step[k] == 0 && prog.x[prog.basic_vars[k]] == 0
            t = 0
        else
            t = primal_step[k] / prog.x[prog.basic_vars[k]]
        end
        if inv_primal_step_length < t
            inv_primal_step_length = t
            i = prog.basic_vars[k]
        end
    end
    if inv_primal_step_length <= 0
        return Unbounded
    end
    primal_step_length = 1 / inv_primal_step_length
    @debug "Primal step length" primal_step_length
    primal_step *= primal_step_length

    # Select leaving variable (already done, i above)
    @debug "Leaving variable" i

    # Compute dual step direction
    ei = zeros(length(prog.z))
    ei[i] += 1
    dual_step = - N' * (B' \ ei[prog.basic_vars])
    @debug "Dual step direction" dual_step

    # Compute dual step length
    idx = indexin(j, prog.nonbasic_vars)[1]
    dual_step_length = prog.z[j] / dual_step[idx]
    @debug "Dual step length" dual_step_length
    dual_step *= dual_step_length

    # Update current dual and primal solutions
    prog.x[j] = primal_step_length
    prog.x[prog.basic_vars] = prog.x[prog.basic_vars] - primal_step
    prog.z[i] = dual_step_length
    prog.z[prog.nonbasic_vars] = prog.z[prog.nonbasic_vars] - dual_step

    # Update basis
    prog.basic_vars = union(setdiff(prog.basic_vars, [i]), [j])
    prog.nonbasic_vars = union(setdiff(prog.nonbasic_vars, [j]), [i])

    return Nonoptimal
end

end
