using GLPK,JuMP, Dualization
using Random
using Distributions
using MathOptInterface
const MOI = MathOptInterface

function getVariables(n,m)
    # n = 2 #number of supliers
    # m = 3 #number of customers
    c = rand(Uniform(0,1), n*m)
    f = rand(Uniform(0,1), n*m)

    s = rand(Uniform(0,1),n)
    d = rand(Uniform(0,1),m)
    sum_s = sum(s)
    sum_d = sum(d)

    s = s* sum_d
    d = d* sum_s
    M = [min(s[div(k, m, RoundUp)],d[if mod(k,m)!= 0 mod(k,m) else m end]) for k in 1:n*m]
    # println(s)
    # println(d)
    # println(M)
    return n,m,c,f,s,d,M
end

function Q4(n,m,c,f,s,d,M,display)
    println("=============== GLPK Solver ========================")
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:(n*m)] >= 0)
    @variable(model, y[1:n*m], Bin)

    @objective(model, Min, sum(c .* x) + sum(f .*y))


    @constraint(model, c1, [sum([x[(i-1)*m+j] for j in 1:m ]) for i in 1:n] .== s)
    @constraint(model, c2, [sum([x[(i-1)*m+j] for i in 1:n ]) for j in 1:m] .== d)
    @constraint(model, c3, x .<= M.*y)

    optimize!(model)
    #println(Lb)
    if display==true
        #print(model)
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
        @show objective_value(model)
        @show value.(x) #  @show value.(x[1:n*m])
        @show value.(y) #  @show value.(y[1:n*m])
    end
    
    return 
end
function Q5(n,m,c,f,s,d,M,tol,display)
    #Benders Decomposition
    println("=============== Benders Decomposition ========================")
    #Master Problem Description
    master = Model(GLPK.Optimizer)

    @variable(master, y[1:n*m], Bin)
    @variable(master,x>=0)
    @objective(master, Min, x + sum(f .*y))

    k = 1 #iteration

    Ub = Inf
    Lb = - Inf
    Optimality_Cut_Counter = 0
    Feasability_Cut_Counter= 0

    while(true)
        # println("===============iter k = ",k,"========================")
        # println("Lb = ",Lb,"    Ub = ",Ub)

        optimize!(master)
        t_status = termination_status(master)# == MOI.Success (if the call was succesful or stopped during)
        p_status = primal_status(master) # == MOI.FEASIBLE_POINT (result if not interupted)

        # println("--Master--")
        # println(t_status)
        # println(p_status)
        #Case where primal is executed and not terminated (add rays and vertices (unbounded problem-> Benders Decomposition))
        if t_status == MOI.INFEASIBLE_OR_UNBOUNDED
            # println("--> unbounded")
            y0 = value.(y)
            x0 = value.(x)
            Lb = - Inf

        elseif p_status == MOI.FEASIBLE_POINT
            # println("--> feasible")
            Lb = objective_value(master)
            y0 = value.(y)
            x0 = value.(x)
        #Cases where the primal is executed and terminated
        else p_status == MOI.INFEASIBLE_POINT
            println("--> stop: problem is infeasible")
            break
        end

        #Benders subproblem (Dual, LP)
        subProblem = Model(GLPK.Optimizer)

        @variable(subProblem, v[1:n])
        @variable(subProblem, u[1:m])
        @variable(subProblem, w[1:(n*m)]>=0)

        @objective(subProblem, Max, sum(v .*s) + sum(u .*d) - sum(M .* w .* y0))
        @constraint(subProblem,c1, [ (v[div(k,m,RoundUp)] + u[if mod(k,m)!= 0 mod(k,m) else m end] - w[k]) for k in 1:n*m] .<= c)

        optimize!(subProblem)  

        t_status_sub = termination_status(subProblem)# == MOI.Success
        p_status_sub = primal_status(subProblem) # == MOI.FEASIBLE_POINT
        # println("--Subproblem--")
        # println(t_status_sub)
        # println(p_status_sub)

        # if unbounded, add the feasibility cut (There is an  extreme ray, adding the corresponding constraint)
        if t_status_sub == MOI.DUAL_INFEASIBLE && p_status_sub == MOI.INFEASIBILITY_CERTIFICATE
            # println("-->extreme rays")
            ve = value.(v)
            ue = value.(u)
            we = value.(w)
            @constraint(master, sum(ve .*s) + sum(ue .*d) - sum(M .* we .* y ) <= 0) #TODO
            Feasability_Cut_Counter += 1
        end
        #if bounded and add optimality cut (add a vertex)

        if p_status_sub == MOI.FEASIBLE_POINT && x0 < objective_value(subProblem) 
            #  println("-->vertice")
            #  println(x0," <" ,  objective_value(subProblem) )
            vv = value.(v)
            uv = value.(u)
            wv = value.(w)
            @constraint(master, sum(vv .*s ) + sum(uv .*d ) - sum(M .* wv .* y)<= x )#TODO
            Ub = sum(f.*y0)+ objective_value(subProblem)
            Optimality_Cut_Counter +=1
        end 
        # we are done
        if p_status_sub == MOI.FEASIBLE_POINT &&  x0 >= objective_value(subProblem) 
            println("finish in ",k," iterations")
            println("y0=",y0)
            println("x0=",x0)
            println("Lb=",Lb," f*= ",objective_value(master)," Ub= ",Ub)
            println("Optimality_Cut_Counter =",Optimality_Cut_Counter)
            println("Feasability_Cut_Counter=",Feasability_Cut_Counter)
            break
        end  

        if abs(Lb-Ub)<1e-2
            println("other breaking condition")
            println("f*=",objective_value(master))
            break
        end

        k += 1
        
    end

    return
end
n,m,c,f,s,d,M = getVariables(2,3)

Q5(n,m,c,f,s,d,M,1e-2,true)
Q4(n,m,c,f,s,d,M,true)