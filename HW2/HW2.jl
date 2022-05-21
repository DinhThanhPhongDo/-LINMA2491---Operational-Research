using JuMP, Gurobi
using MathOptInterface
using LinearAlgebra
const MOI = MathOptInterface


function Lshaped(c, A, b, q, W, h, T, p;tol=1e-2,max_iter=1e3,display=false)
    """
        L-shaped Decomposition for an instance of the Farmer's Problem
        
        parameters:
        -----------

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

    dim_x = length(c)

    nb_scen = length(p)


    gurobi_env = Gurobi.Env()

    #Master Problem
    master = Model(optimizer_with_attributes( () -> Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
    set_optimizer_attribute(master,"presolve",2)
    @variable(master, x[1:dim_x]>= 0)
    @variable(master,theta>=0)
    @objective(master, Min, theta + sum(c .*x))




    k = 1 #iterations

    Ub = Inf
    Lb = - Inf
    Optimality_Cut_Counter = 0
    Feasability_Cut_Counter= 0

    keepgoing = true
    x_opt = 0
    theta_opt = 0
    timeMaster = []
    timeSubProblem = []

    u_opt = zeros(length(q), 1)

    while(keepgoing)

        push!(timeMaster, @elapsed optimize!(master))
        t_status = termination_status(master)
        p_status = primal_status(master) 

        #If master is executed and not terminated
        if t_status == MOI.INFEASIBLE_OR_UNBOUNDED

            println("=======================================")
            println("--> stop: problem is infeasible or unbounded")
            break
        #If master is executed and an optimal solution is found
        elseif p_status == MOI.FEASIBLE_POINT

            Lb = objective_value(master)
            x_opt = value.(x)
            theta_opt = value.(theta)
            println("============== MASTER ==============")
            println("x_opt = ", x_opt)
            println("theta_opt = ", theta_opt)


        ##If master is executed and is infeasible
        else p_status == MOI.INFEASIBLE_POINT

            println("=======================================")
            println("--> stop: problem is infeasible")
            break
        end

        



        #Subproblem
        #-------------------------------------------------------------------------------------------------------------------

        u_opt = zeros(length(q), 1)
        subfeasible = true
        sub_objval = 0
        scen = 1
        dim_pi = length(h)

        pi_vertex= zeros(nb_scen, dim_pi)

        while ((scen <= nb_scen) && (subfeasible))
            subProblem = Model(optimizer_with_attributes( () -> Gurobi.Optimizer(gurobi_env),"OutputFlag"=>0))
            set_optimizer_attribute(subProblem,"presolve",0);
            @variable(subProblem, pi[1:dim_pi]>= 0)


            prod = h - (T[scen,:,:]*x_opt)

            println("prod = ", prod)
            @objective(subProblem, Max,sum(pi.*(prod)))
            @constraint(subProblem,c1, transpose(pi)*W .<= transpose(q))

            push!(timeSubProblem, @elapsed optimize!(subProblem))

            pi_opt = value.(pi)

            println("pi_opt = ", pi_opt)

            t_status_sub = termination_status(subProblem)
            p_status_sub = primal_status(subProblem) 


            u_scen_opt = zeros(length(q))


            # if unbounded, add the feasibility cut (There is an  extreme ray, adding the corresponding constraint)
            if  keepgoing && (t_status_sub == MOI.DUAL_INFEASIBLE && p_status_sub == MOI.INFEASIBILITY_CERTIFICATE)

                sigma = value.(pi)

                @constraint(master, transpose(sigma)*(h - T[scen,:,:]*x) <= 0) 
                Feasability_Cut_Counter += 1
                subfeasible = false

                println("-------------- FEASIBILITY CUT-------------------")
                println("number = ", Feasability_Cut_Counter)

                println("scenario :", scen)
                println("sigma = ", sigma)
                
            
            #if feasible remember pi
            elseif keepgoing && (p_status_sub == MOI.FEASIBLE_POINT)

                pi_vertex[scen,:] = pi_opt
                if(has_duals(subProblem))
                    u_scen_opt = transpose(-dual.(c1))  
                    println("u_scen_opt = ", u_scen_opt)
                    println("size u_opt = ",size(u_opt))
                    u_opt = u_opt + p[scen].*u_scen_opt

                end


                sub_objval = sub_objval + p[scen]*objective_value(subProblem)
                println("sub_objval_temp = ", sub_objval)

            end
        scen = scen + 1

        end
        #----------------------------------------------------------------------------------


        # Add optimality cut if no unfeasible subprobs and theta <= obj val of subprob (S)
        println("subobj = ", sub_objval)
        println("theta_opt = ", theta_opt)

        if ((subfeasible == true) &&  (theta_opt <= sub_objval))

            sum1 = 0
            sum2 = 0
            for sc in 1: nb_scen


                sum1 = sum1 .+ p[sc].*transpose(pi_vertex[sc,:])*h
                sum2 = sum2 .+ p[sc].*transpose(pi_vertex[sc,:])*T[sc,:,:]
            end
            @constraint(master, sum1 - sum2*x .<= theta)

            println("sum1 = ", sum1, " ; sum2 = ", sum2)
            Ub = sum(c.*x_opt)+ sub_objval
            Optimality_Cut_Counter +=1 

            println("-------------- OPTIMALITY CUT-------------------")
            println("number = ", Optimality_Cut_Counter)

            for sc in 1:nb_scen
                println("pi_vertex [", sc, "] = ", pi_vertex[sc,:])
            end


        # stopping conditions : Optimality or too may iterations
        elseif subfeasible &&  theta_opt >= sub_objval - 1e-10 


            if display

                println("=======================================")
                println("Optimal solution found after ",k," iterations")
                println("Number of feasability cuts=",Feasability_Cut_Counter)
                println("Number of optimality cuts =",Optimality_Cut_Counter)
                println("Lb=",Lb," f*= ",objective_value(master)," Ub= ",Ub)
            end
            keepgoing = false 
 
        elseif (abs(Lb-Ub)<tol || k>max_iter)
            


            if display

                println("=======================================")
                println("stopped after",k," iterations because tolerances or numbers of iteration exceed max iteration")
                println("Optimality_Cut_Counter =",Optimality_Cut_Counter)
                println("Feasability_Cut_Counter=",Feasability_Cut_Counter)
                println("Lb=",Lb," f*= ",objective_value(master)," Ub= ",Ub)
            end

            keepgoing = false
        end

        if display
            if k == 1

                
                println("===============iter k = ",k,"========================")
                println("Lb = ",Lb,"    Ub = ",Ub)
                println("=====================================================")
            end
            i=20
            if mod(k,i)==0

                println("===============iter k = ",k,"========================")
                println("Lb = ",Lb,"    Ub = ",Ub)
                println("=====================================================")


                println("time for master = ",sum(timeMaster[end-i+1:end])/i,"    time for subProblem =",sum(timeSubProblem[end-i+1:end])/i)
            end
        end

        k += 1
        
    end

    return objective_value(master), x_opt, u_opt, timeMaster,timeSubProblem
end

function compareTime()


    c = [150; 230; 260]

    A = [1 1 1]

    b = 500

    q = [238; 210; -170; -150; -36; -10]

    W = [1 0 -1 0 0 0; 0 1 0 -1 0 0; 0 0 0 0 -1 -1; 0 0 0 0 -1 0]

    h = [200; 240; 0; -6000]

    T1 = [3 0 0; 0 3.6 0; 0 0 24; 0 0 0]
    T2 = [2.5 0 0; 0 3 0; 0 0 20; 0 0 0]
    T3 = [2 0 0; 0 2.4 0; 0 0 16; 0 0 0]

    T = zeros(3, 4, 3)
    T[1,:,:] = T1
    T[2,:,:] = T2
    T[3,:,:] = T3

    p = [1/3 1/3 1/3]

    println("BEGIN")

    timeLshaped = @elapsed obj2, x_opt2, u_opt2,timeMaster,timeSubProblem = Lshaped(c, A, b, q, W, h, T, p;tol=1e-2,max_iter=1e3,display=true);
    println("timeLshaped = ",timeLshaped)
    println("objective value optimal = ", obj2)
    println("x optimal = ", x_opt2)
    println("u optimal = ", u_opt2)



end




compareTime()


