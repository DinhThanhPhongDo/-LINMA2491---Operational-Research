using GLPK,JuMP, Gurobi
using Random, Distributions
using MathOptInterface
using LinearAlgebra
using Plots
const MOI = MathOptInterface

function getVariables(n,m)
    """
    return an random instance with n suppliers and m customers.
        parameters
        ----------
        n   : Integer
              Number of suppliers
        m   : Integer
              Number of customers
        returns
        -------
        n   : Integer
              Number of suppliers
        m   : Integer
              Number of customers
        c   : Array of size n*m
              Cost of transporting a quantity x between suppliers and customers
        f   : Array of size n*M
              Cost for opening a routes between suppliers and customers
        s   : Array of size n
              Quantity asked by each customer
        d   : Array of size m
              Quantity produced by each supplier
        M   : Array of size n*m 
              Upper bound of the quantity x between suppliers and customers
    """
    c = rand(Uniform(1,5), n*m)
    f = rand(Uniform(10,50), n*m)

    s = rand(Uniform(50,1000),n)
    d = rand(Uniform(1,10),m)
    sum_s = sum(s)
    sum_d = sum(d)
    s = s* sum_d
    d = d* sum_s

    M = [min(s[div(k, m, RoundUp)],d[mod(k-1,m)+1]) for k in 1:n*m]

    return n,m,c,f,s,d,M
end

function Q4(n,m,c,f,s,d,M;display=false)
    """
        Direct Solver for an instance of the Fixed Charge Trasportation Problem
        
        parameters:
        -----------
        n   : Integer
              Number of suppliers
        m   : Integer
              Number of customers
        c   : Array of size n*m
              Cost of transporting a quantity x between suppliers and customers
        f   : Array of size n*M
              Cost for opening a routes between suppliers and customers
        s   : Array of size n
              Quantity asked by each customer
        d   : Array of size m
              Quantity produced by each supplier
        M   : Array of size n*m 
              Upper bound of the quantity x between suppliers and customers
        display : Boolean
                  display or not the results (termination status, primal status, dual status and objective value)
        returns:
        --------
        objective_value(model) : Return the objective value of this instance of the Fixed Charge Trasportation Problem
    """
    
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
    @variable(model, x[1:(n*m)] >= 0)
    @variable(model, y[1:n*m], Bin)

    @objective(model, Min, sum(c.*x)+ sum(f.*y))


    @constraint(model, c1, [sum([x[(i-1)*m+j] for j in 1:m ]) for i in 1:n] .== s)
    @constraint(model, c2, [sum([x[(i-1)*m+j] for i in 1:n ]) for j in 1:m] .== d)
    @constraint(model, c3, x .<= M.*y)

    optimize!(model)
    x_opt = value.(x)
    y_opt = value.(y)

    if display==true
        #print(model)
        println("======================== Direct Solver ========================")
        @show termination_status(model)
        @show primal_status(model)
        @show dual_status(model)
        @show objective_value(model)
    end
    
    return objective_value(model), x_opt, y_opt
end
function Q5(n,m,c,f,s,d,M;tol=1e-2,max_iter=1e3,display=false)
    """
        Benders Decomposition for an instance of the Fixed Charge Trasportation Problem
        
        parameters:
        -----------
        n   : Integer
              Number of suppliers
        m   : Integer
              Number of customers
        c   : Array of size n*m
              Cost of transporting a quantity x between suppliers and customers
        f   : Array of size n*M
              Cost for opening a routes between suppliers and customers
        s   : Array of size n
              Quantity asked by each customer
        d   : Array of size m
              Quantity produced by each supplier
        M   : Array of size n*m 
              Upper bound of the quantity x between suppliers and customers
        tol : Real number 
              stoping conditions if the distance between the upper and lower bound is smaller than the parameter tol
        max_iter: Integer
                  stopping condition if the number of iteration performed exceed the parameter max_iter
        display : Boolean
                  display iterations, bounds, and final result.
        returns:
        --------
        objective_value(model) : Return the objective value of this instance of the Fixed Charge Trasportation Problem
    """
    gurobi_env = Gurobi.Env()

    #Master Problem
    master = Model(optimizer_with_attributes( () -> Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
    set_optimizer_attribute(master,"presolve",2)
    set_optimizer_attribute(master,"MIPFocus",1)
    @variable(master, y[1:n*m], Bin)
    @variable(master,x>=0)
    @objective(master, Min, x + sum(f .*y))




    k = 1 #iterations

    Ub = Inf
    Lb = - Inf
    Optimality_Cut_Counter = 0
    Feasability_Cut_Counter= 0

    keepgoing = true
    x_opt = 0
    y_opt = 0
    timeMaster = []
    timeSubProblem = []
    while(keepgoing)

        push!(timeMaster, @elapsed optimize!(master))
        t_status = termination_status(master)
        p_status = primal_status(master) 

        #If master is executed and not terminated
        if t_status == MOI.INFEASIBLE_OR_UNBOUNDED

            y_opt = value.(y)
            theta = value.(x)
            Lb = - Inf
        #If master is executed and an optimal solution is found
        elseif p_status == MOI.FEASIBLE_POINT

            Lb = objective_value(master)
            y_opt = value.(y)
            theta = value.(x)
            
        ##If master is executed and is infeasible
        else p_status == MOI.INFEASIBLE_POINT

            println("=======================================")
            println("--> stop: problem is infeasible")
            break
        end

        

        #Subproblem
        subProblem = Model(optimizer_with_attributes( () -> Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
        set_optimizer_attribute(subProblem,"presolve",0);
        @variable(subProblem, v[1:n])
        @variable(subProblem, u[1:m])
        @variable(subProblem, w[1:(n*m)]>=0)
        @objective(subProblem, Max, sum(v .*s) + sum(u .*d) - sum(M .* w .* y_opt))
        @constraint(subProblem,c1, [ (v[div(k,m,RoundUp)] + u[mod(k-1,m)+1] - w[k]) for k in 1:n*m] .<= c)

        push!(timeSubProblem, @elapsed optimize!(subProblem))

        t_status_sub = termination_status(subProblem)
        p_status_sub = primal_status(subProblem) 



        # stopping conditions
        if p_status_sub == MOI.FEASIBLE_POINT &&  theta >= objective_value(subProblem) - 1e-10 

            if(has_duals(subProblem))

                x_opt =  -dual.(c1)
            end
            if display

                println("=======================================")
                println("Optimal solution found after ",k," iterations")
                println("Number of feasability cuts=",Feasability_Cut_Counter)
                println("Number of optimality cuts =",Optimality_Cut_Counter)
                println("Lb=",Lb," f*= ",objective_value(master)," Ub= ",Ub)
            end
            keepgoing = false 
 
        elseif (abs(Lb-Ub)<tol || k>max_iter)
            
            if(has_duals(subProblem))

                x_opt =  -dual.(c1)
            end

            if display

                println("=======================================")
                println("stopped after",k," iterations because tolerances or numbers of iteration exceed max iteration")
                println("Optimality_Cut_Counter =",Optimality_Cut_Counter)
                println("Feasability_Cut_Counter=",Feasability_Cut_Counter)
                println("Lb=",Lb," f*= ",objective_value(master)," Ub= ",Ub)
            end

            keepgoing = false
        end

        # if unbounded, add the feasibility cut (There is an  extreme ray, adding the corresponding constraint)
        if  keepgoing && (t_status_sub == MOI.DUAL_INFEASIBLE && p_status_sub == MOI.INFEASIBILITY_CERTIFICATE)

            ve = value.(v)
            ue = value.(u)
            we = value.(w)
            @constraint(master, sum(ve .*s) + sum(ue .*d) - sum(M .* we .* y ) <= 0) #TODO
            Feasability_Cut_Counter += 1
        
        #if bounded and add optimality cut (add a vertex)
        elseif keepgoing && (p_status_sub == MOI.FEASIBLE_POINT && theta < objective_value(subProblem))

            vv = value.(v)
            uv = value.(u)
            wv = value.(w)
            @constraint(master, sum(vv .*s ) + sum(uv .*d ) - sum(M .* wv .* y)<= x )
            Ub = sum(f.*y_opt)+ objective_value(subProblem)
            Optimality_Cut_Counter +=1 
        end

        if display
            if k == 1

                println("======================== Benders Decomposition ========================")
                println("===============iter k = ",k,"========================")
                println("Lb = ",Lb,"    Ub = ",Ub)
            end
            i=20
            if mod(k,i)==0

                println("===============iter k = ",k,"========================")
                println("Lb = ",Lb,"    Ub = ",Ub)

                println("time for master = ",sum(timeMaster[end-i+1:end])/i,"    time for subProblem =",sum(timeSubProblem[end-i+1:end])/i)
            end
        end

        k += 1
        
    end

    return objective_value(master), x_opt, y_opt, timeMaster,timeSubProblem
end


function compareTime(n,m)

    n,m,c,f,s,d,M = getVariables(n,m)


    time1 = @elapsed obj1, x_opt1, y_opt1 = Q4(n,m,c,f,s,d,M,display=false);

    n,m,c,f,s,d,M = getVariables(n,m)

    time3 = @elapsed obj1, x_opt1, y_opt1= Q4(n,m,c,f,s,d,M,display=true);

    time4 = @elapsed obj2, x_opt2, y_opt2,timeMaster,timeSubProblem = Q5(n,m,c,f,s,d,M,tol=1e-2,display=true);
    x = range(1,length(timeMaster))
    plot(x,timeMaster)
    plot!(x,timeSubProblem)
    display(plot)
    return time1, time3, time4, abs(obj1 -obj2), norm(x_opt2 .- x_opt1), norm(y_opt2 .- y_opt1)
end


time1, time3, time4, diff_f, diff_x,diff_y = compareTime(5,5);
println("time1 = ",time1)
println("time3 = ",time3)
println("time4 = ",time4)
println("diff_f = ",diff_f)
println("diff_x = ",diff_x)
println("diff_y = ",diff_y)






