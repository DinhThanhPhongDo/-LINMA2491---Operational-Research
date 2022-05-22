# Implementation of the Single-Cut L-Shaped Algorithm for Finite Support R.V.
# LINMA2491 - Operations Research - Homework 2
# Group Members: Phong Do, Julien Donadello, Matthias Geurts, Eduardo Vannini
# Date: 20-04-2021


export L_Shape_Algorithm

########################### Parameters ##############################

# Optimizer name, choose among "GLPK", "Gurobi", "Ipopt" and "CPLEX"
global OPTIMIZE = "GLPK"

######################################################################

using JuMP, Printf

if OPTIMIZE == "Gurobi"
    using Gurobi
elseif OPTIMIZE == "Ipopt"
    using Ipopt
elseif OPTIMIZE == "CPLEX"
    using CPLEX
else
    using GLPK
end

function get_model()
    if OPTIMIZE == "Gurobi"
        return Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    elseif OPTIMIZE == "Ipopt"
        return Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    elseif OPTIMIZE == "CPLEX"
        return Model(optimizer_with_attributes(CPLEX.Optimizer, "print_level" => 0))
    else
        return Model(optimizer_with_attributes(GLPK.Optimizer))
    end
end

function step_1(c, A_ge, b_ge, A_le, b_le, A_eq, b_eq, D, d, E, e)
    # c is a 1D vector
    # A_* are 2D matrices
    # b_* are 1D vectors
    # D is a 2D matrix: element D[l, :] is a 1D vector for l = 1..r
    # d is a 1D vector: element d[l] is a scalar for l = 1..r
    # E is a 2D array: element E[l, :] is a 1D vector for l = 1..s
    # e is a 1D vector: element e[l] is a scalar for l = 1..s

    m = get_model()

    n = size(A_ge, 2)  # Number of first-stage variables in vector 'x'
    r = size(D, 1)  # Number of feasibility cuts
    s = size(E, 1)  # Number of optimality cuts

    if s == 0
        s_coeff = 0
    else
        s_coeff = 1
    end

    # Note: Is there a need to add "x >= 0" if it is not part of the constraints
    # of the main problem ? Let's say no: if this constraint must be added, then
    # it must be done explicitly through A_ge and b_ge
    #
    @variable(m, x[1:n] >= 0)
    #@variable(m, x[1:n])

    @variable(m, theta)

    @constraint(m, A_ge * x .>= b_ge)
    @constraint(m, A_le * x .<= b_le)
    @constraint(m, A_eq * x .== b_eq)

    @constraint(m, D * x .>= d)

    # Not using matrix formulation for this constraint due to a strange
    # StackOverflowError bug in JuMP
    # https://stackoverflow.com/questions/67027176/julia-jump-stackoverflow-error-when-defining-a-constraint-on-linear-matrix-inequ
    @constraint(m, [l = 1:s], E[l, :]' * x + theta * s_coeff >= e[l])

    @objective(m, Min, c' * x + theta * s_coeff)

    #println()
    #println("    [Step 1] Model:")
    #println(m)

    optimize!(m)

    #println()
    #@printf("    [Step 1] Optimal objective value: z = %f\n", objective_value(m))

    if s == 0
        return value.(x), -Inf, objective_value(m)
    else
        return value.(x), value.(theta), objective_value(m)
    end
end


function step_2(W_ge, T_ge, h_ge, W_le, T_le, h_le, W_eq, T_eq, h_eq, K, xv, D, d)
    # W_* are 2D matrices
    # T_* are 3D arrays such that T_*[k, :, :] are 2D matrices for k = 1..K
    # h_* are 2D matrices such that h_*[k, :] are 1D vectors for k = 1..K
    # K is a scalar
    # xv is a 1D vector
    # D is a 2D matrix: element D[l, :] is a 1D vector for l = 1..r
    # d is a 1D vector: element d[l] is a scalar for l = 1..r

    n_ge = size(W_ge, 1)  # Number of >= constraints
    n_le = size(W_le, 1)  # Number of <= constraints
    n_eq = size(W_eq, 1)  # Number of == constraints
    n_cons = n_ge + n_le + n_eq  # Total number of second-stage constraints
    n_yvars = size(W_ge, 2)  # Number of second-stage decision variables (y)
    r = size(D, 1)

    k = 1
    while k <= K
        m = get_model()

        @variable(m, v_plus[1:n_cons] >= 0)
        @variable(m, v_minus[1:n_cons] >= 0)

        # Note: Is there a need to add "y >= 0" if it is not part of the constraints
        # of the main problem ? Let's say no: if this constraint must be added, then
        # it must be done explicitly through W_ge and h_ge
        #
        @variable(m, y[1:n_yvars] >= 0)
        #@variable(m, y[1:n_yvars])

        #@printf("\nW_ge = \n")
        #show(stdout, "text/plain", W_ge)
        #@printf("\nT_ge[k, :, :] = \n")
        #show(stdout, "text/plain", T_ge[k, :, :])
        #@printf("\nxv = \n")
        #show(stdout, "text/plain", xv)
        #@printf("\nh_ge[k, :] = \n")
        #show(stdout, "text/plain", h_ge[k, :])
        #@printf("\n")
        #@printf("n_cons = %i\n",n_cons)
        #@printf("n_yvars = %i\n",n_yvars)
        #@printf("n_yvars = %i\n",n_ge)
        #@printf("\n")
        
        @constraint(m, cons_ge, W_ge * y + v_plus[1:n_ge] - v_minus[1:n_ge] .>= h_ge[k, :] - T_ge[k, :, :] * xv)
        @constraint(m, cons_le, W_le * y + v_plus[(n_ge+1):(n_ge+n_le)] - v_minus[(n_ge+1):(n_ge+n_le)] .<= h_le[k, :] - T_le[k, :, :] * xv)

        # println(h_le)
        # println(h_le[k, :])
        # println(T_le[k, :, :])
        # println(xv)

        @constraint(m, cons_eq, W_eq * y + v_plus[(n_ge+n_le+1):n_cons] - v_minus[(n_ge+n_le+1):n_cons] .== h_eq[k, :] - T_eq[k, :, :] * xv)

        @objective(m, Min, ones(n_cons)' * v_plus + ones(n_cons)' * v_minus)

        #println()
        #println("    [Step 2] Model:")
        #println(m)

        optimize!(m)

        #println()
        #@printf("    [Step 2] Objective value: %f\n", objective_value(m))
        #println()
        #@printf("    [Step 2] Variables values: \n    v_plus: %s\n    v_minus: %s\n    y: %s\n", value.(v_plus), value.(v_minus), value.(y))

        if objective_value(m) > 0
            # Create D_{r + 1} and d_{r + 1}
            D_next = dual.(cons_ge)' * T_ge[k, :, :] + dual.(cons_le)' * T_le[k, :, :] + dual.(cons_eq)' * T_eq[k, :, :]
            d_next = dual.(cons_ge)' * h_ge[k, :] + dual.(cons_le)' * h_le[k, :] + dual.(cons_eq)' * h_eq[k, :]

            # Rebuild D and d and add D_{r + 1} and d_{r + 1} at the end
            # new_D = zeros(r + 1, size(D, 2))
            # new_d = zeros(r + 1)
            # l = 1
            # while l <= r
            #     new_D[l, :] = D[l, :]
            #     new_d[l] = d[l]
            #     l += 1
            # end
            # new_D[r + 1, :] = D_next
            # new_d[r + 1] = d_next
            new_D = [D; D_next]
            new_d = [d; d_next]

            # Return value indicates if a new feasibility cut was added or not
            #println()
            #@printf("    [Step 2] Added feasibility cut %s * x >= %s\n", D_next, d_next)
            return true, new_D, new_d
        end
        k += 1
    end
    return false, D, d
end


function step_3(W_ge, T_ge, h_ge, W_le, T_le, h_le, W_eq, T_eq, h_eq, q, p, K, xv, theta_v, E, e)
    s = size(E, 1)
    n_yvars = size(W_ge, 2)

    # Allocating memory for E_{s + 1} and e_{s + 1}
    E_next = zeros(size(E, 2))
    e_next = 0

    # Using while loop because I don't get the scope issues with for loops in Julia..
    k = 1
    while k <= K
        m = get_model()

        # Note: Is there a need to add "y >= 0" if it is not part of the constraints
        # of the main problem ? Let's say no: if this constraint must be added, then
        # it must be done explicitly through W_ge and h_ge
        #
        @variable(m, y[1:n_yvars] >= 0)
        #@variable(m, y[1:n_yvars])

        @constraint(m, cons_ge, W_ge * y .>= h_ge[k, :] - T_ge[k, :, :] * xv)
        @constraint(m, cons_le, W_le * y .<= h_le[k, :] - T_le[k, :, :] * xv)
        @constraint(m, cons_eq, W_eq * y .== h_eq[k, :] - T_eq[k, :, :] * xv)

        @objective(m, Min, q[k, :]' * y)

        #println()
        #@printf("    [Step 3] Computing scenario k = %d\n", k)
        #@printf("    [Step 3, k = %d] Model:\n", k)
        # println(m)

        optimize!(m)

        @printf("    [Step 3, k = %d] Dual variables of 'ge' constraints: %s\n", k, dual.(cons_ge))
        @printf("    [Step 3, k = %d] Dual variables of 'le' constraints: %s\n", k, dual.(cons_le))
        #@printf("    [Step 3, k = %d] Dual variables of 'eq' constraints: %s\n", k, dual.(cons_eq))
        #
        #@printf("    [Step 3, k = %d] h_k vector for 'ge' constraints: %s\n", k, h_ge[k, :])
        #@printf("    [Step 3, k = %d] h_k vector for 'le' constraints: %s\n", k, h_le[k, :])
        #@printf("    [Step 3, k = %d] h_k vector for 'eq' constraints: %s\n", k, h_eq[k, :])
        #
        #@printf("    [Step 3, k = %d] T_k matrix for 'ge' constraints: %s\n", k, T_ge[k, :, :])
        #@printf("    [Step 3, k = %d] T_k matrix for 'le' constraints: %s\n", k, T_le[k, :, :])
        #@printf("    [Step 3, k = %d] T_k matrix for 'eq' constraints: %s\n", k, T_eq[k, :, :])
        #
        @printf("    [Step 3, k = %d] Optimal objective value for the sub-problem: %f\n", k, objective_value(m))
        @printf("    [Step 3, k = %d] x for this subproblem: %s\n", k, xv)

        #@printf("\nT_ge[%i, :, :] = %s\n",k,T_ge[k, :, :])
        #@printf("\nT_ge[%i, :, :] = %s\n",k,dual.(cons_ge)')
        #@printf("\nT_le[%i, :, :] = %s\n",k,T_le[k, :, :])
        #@printf("\nT_ge[%i, :, :] = %s\n",k,dual.(cons_le)')
        #@printf("\nT_eq[%i, :, :] = %s\n",k,T_eq[k, :, :])
        #@printf("\nT_ge[%i, :, :] = %s\n",k,dual.(cons_eq)')
        #@printf("\np[%i] = %s\n",k,p[k])

        E_next += p[k] * (dual.(cons_ge)' * T_ge[k, :, :] + dual.(cons_le)' * T_le[k, :, :] + dual.(cons_eq)' * T_eq[k, :, :])'
        e_next += p[k] * (dual.(cons_ge)' * h_ge[k, :] + dual.(cons_le)' * h_le[k, :] + dual.(cons_eq)' * h_eq[k, :])

        k += 1
    end

    # Check for optimality
    wv = e_next - E_next' * xv
    @printf("    [Step 3] theta_v = %f , wv = %f \n", theta_v, wv)
    if theta_v >= wv || abs(theta_v - wv) <= 1e-3
        # Return value says if an optimality cut was added or not
        # If not, it means xv is an optimal solution and L-shaped algo can stop
        return false, E, e
    end

    #println()
    @printf("    [Step 3] Added optimality cut %s * x + theta >= %s\n", E_next, e_next)

    # Rebuild E and e and add E_{s + 1} and e_{s + 1} at the end
    # new_E = zeros(s + 1, size(E, 2))
    # new_e = zeros(s + 1)
    # l = 1
    # while l <= s
    #     new_E[l, :] = E[l, :]
    #     new_e[l] = e[l]
    #     l += 1
    # end
    # new_E[s + 1, :] = E_next
    # new_e[s + 1] = e_next
    new_E = [E; E_next']
    new_e = [e; e_next]
    # println(new_E)
    # println(new_e)

    return true, new_E, new_e
end

# L-Shaped Algorithm
function L_Shape_Algorithm(c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le)

    # L-Shaped Algorithm parameters
    D = zeros(0, size(c, 1))
    d = zeros(0)

    E = zeros(0, size(c, 1))
    e = zeros(0)

    z_opt = 0
    x_opt = -1
    theta_opt = -1
    done = false

    # Max iter for debug purposes
    use_max_iter = false
    max_iter = 60

    iter = 0
    while !done && (iter < max_iter || !use_max_iter)

        #println("----------------------------------------------------")
        #@printf("L-Shaped Algorithm: Starting iteration %d\n", iter + 1)

        # Step 1
        @printf("[LSA iter %d] Running step 1\n", iter + 1)
        x_opt, theta_opt, z_opt = step_1(c, A_ge, b_ge, A_le, b_le, A_eq, b_eq, D, d, E, e)

        #print()
        #@printf("    x_opt: %s\n", x_opt)
        #@printf("    theta_opt: %s\n", theta_opt)

        # Step 2
        #println("\n\n")
        @printf("[LSA iter %d] Running step 2\n", iter + 1)
        flag_2, D, d = step_2(W_ge, T_ge, h_ge, W_le, T_le, h_le, W_eq, T_eq, h_eq, K, x_opt, D, d)

        if flag_2
            @printf("    Added feasibility cut: going back to step 1\n")
        else
            @printf("    No feasibility cut added: going to step 3\n")
        end
        #@printf("    D matrix: %s\n", D)
        #@printf("    d vector: %s\n", d)

        if !flag_2
            # Step 3
            #println("\n\n")
            @printf("[LSA iter %d] Running step 3\n", iter + 1)
            flag_3, E, e = step_3(W_ge, T_ge, h_ge, W_le, T_le, h_le, W_eq, T_eq, h_eq, q, p, K, x_opt, theta_opt, E, e)
            done = !flag_3

            if flag_3
                @printf("    Added optimality cut: going back to step 1\n")
            else
                @printf("    No optimality cut added: end of algorithm\n")
            end
            #@printf("    E matrix: %s\n", E)
            #@printf("    e vector: %s\n", e)
        end

        iter += 1
    end

    if use_max_iter && iter >= max_iter
        println()
        println("Warning: max_iter reached, result could be suboptimal")
    end

    println()
    @printf("Result obtained in %d iterations\n", iter)
    @printf("Optimal solution is x_opt = %s\n", x_opt)

    return x_opt, z_opt
end