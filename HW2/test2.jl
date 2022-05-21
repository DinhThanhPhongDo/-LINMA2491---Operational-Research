using JuMP, GLPK

m = Model(GLPK.Optimizer)
@variable(m, x[1:2])
@constraint(m, c[i=1:2], x[i] >= 2)
# Uncommenting the next line leads to an optimal solution instead of an unbounded problem...
# @constraint(m, c2[i=1:2], x[i] <= 5)
@objective(m, Max, sum(x))
optimize!(m)

if termination_status(m) == DUAL_INFEASIBLE && primal_status(m) == INFEASIBILITY_CERTIFICATE
    # Retrieving an extreme ray if the primal problem is unbounded
    an_extreme_ray = value.(x)
    println("extreme ray = ", an_extreme_ray)
    # Check the values returned. Could Other extreme rays  have been retuned ?
elseif termination_status(m) == OPTIMAL
    # Retrieving an optimal solution
    an_optimal_point = value.(x)
end
