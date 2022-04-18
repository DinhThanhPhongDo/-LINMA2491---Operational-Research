#ligne 20
import Pkg; Pkg.add("Plots")
using LinearAlgebra
function f(k,m)

	if mod(k,m)!= 0 y = mod(k,m) else y = m end
	y
end

function f2(k,m)
	y = mod(k-1,m)+1
	# y = ((k-1)%m) + 1
end


#Ligne 33
# dot(c,x)



a = [1 1 1]
b = [3 4 5]
println(sum(a.*b))
println(dot(a,b))

# in benders, optimal x?
# why k<3 in while loop ?
# tol instead of 10e-2

#=

For the larger case initially in scope, you can also  on the convergence status you obtain, to document how slow is the Benders decomposition for your problem, the final gap you have, or the quality of the solution you obtain when stopping the process before convergence (comparing objective values with the objective for the optimal one).

Last but not least, as you know, floating-point arithmetic issues require you to be cautious when you compare numbers to trigger actions such as adding Benders cuts, etc. Consider for example taking into account very small tolerances "epsilon" (e.g. ~1e-6 or 1e-7, etc, depending on the context) when checking conditions such as "a >= b".

=#