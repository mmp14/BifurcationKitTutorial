using Revise, Plots
using BifurcationKit

const BK = BifurcationKit

# Define the function σ
function σ(v, e0 = 2.5, v0 = 6, r=0.56)
    2 * e0 / (1 + exp(r * (v0 - v)))
end

# vector field
function jansenrit(z, param)
	(;A, a, B, b, C1, C2, C3, C4, p) = param #These are the parameters. p will be our control paramter
	y0, y3, y1, y4, y2, y5 = z #These are our functions
	[
        #We define the derivative functions in the same order as above
	    y3
        A * a * σ(y1 - y2) - 2 * a * y3 - a^2 * y0
        y4
        A * a * (p + C2 * σ(C1 * y0)) - 2 * a * y4 - a^2 * y1
        y5
        B * b * C4 * σ(C3 * y0) - 2 * b * y5 - b^2 * y2
	]
end

# parameter values
C= 135
par = (A = 3.25, a = 100, B = 22, b = 50, C1 = C, C2 = 0.8 * C, C3 = 0.25 * C, C4 = 0.25 * C, p = 50.0)

# initial condition
z0 = [0.0, 0.0, 15.0, 0.0, 10.0, 0.0]

prob = BifurcationProblem(jansenrit, z0, par, (@optic _.p);
	record_from_solution = (x, p; k...) -> (y0 = x[1], y3 = x[2], y1 = x[3], y4 = x[4], y2 = x[5], y5 = x[6], y = x[3]-x[5]),)



optnewton = NewtonPar(tol = 1e-11, verbose = true)

 # continuation options, we limit the parameter range for p
opts_br = ContinuationPar(p_min = -50.0, p_max = 400.0, max_steps = 8000, newton_options = optnewton)

# continuation of equilibria
br = continuation(prob, PALC(), opts_br;
	# we want to compute both sides of the branch of the initial
	# value of p = 50
	bothside = true, normC = norminf)


###############
# Plot wrt y0 #
###############
colour_stability = [stable ? :green : :red for stable in br.stable]

# We get the bifurcation points and Hopf points
bp =  filter(x -> x.type==:bp, br.specialpoint)
hopf = filter(x -> x.type==:hopf, br.specialpoint)

# We get the x and y coordinate for each bifurcation point
y_bp = map(p -> p.x, bp) # This returns a vector of vectors of dimensions number of points * number of values (1 for each function we have recorded)
params_bp =  map(p -> p.param, bp)

# We get the x and y coordinate for each Hopf point
y_hopf = map(p -> p.x, hopf)
params_hopf =  map(p -> p.param, hopf)

plot(br.param, br.y0, xlabel = "p", ylabel = "y0", color = colour_stability, primary = false)
scatter!(params_bp, getindex.(y_bp,1), c="blue", ms=4, label = "Bifurcation Point" ) #we use getindex.(y_bp,1) because we want the coordinate for y0
scatter!(params_hopf, getindex.(y_hopf,1), c="red", ms=4, label = "Hopf Point" )

##############
# Plot wrt y #
##############
# We get the bifurcation points and Hopf points
bp =  filter(x -> x.type==:bp, br.specialpoint)
hopf = filter(x -> x.type==:hopf, br.specialpoint)

# We get the x and y coordinate for each bifurcation point
y_bp = map(p -> p.x, bp) # This returns a vector of vectors of dimensions number of points * number of values (1 for each function we have recorded)
params_bp =  map(p -> p.param, bp)

# We get the x and y coordinate for each Hopf point
y_hopf = map(p -> p.x, hopf)
params_hopf =  map(p -> p.param, hopf)

plot(br.param, br.y, xlabel = "p", ylabel = "y = y1 - y2", color = colour_stability, primary = false)
scatter!(params_bp, getindex.(y_bp, 3) - getindex.(y_bp, 5), c="blue", ms=4, label = "Bifurcation Point" ) #we use getindex.(y_bp,1) because we want the coordinate for y0
scatter!(params_hopf, getindex.(y_hopf,3) - getindex.(y_hopf, 5), c="red", ms=4, label = "Hopf Point" )

###################
# Periodic orbits #
###################

opts_br_po = ContinuationPar(p_min =-50.0, p_max = 400.0, max_steps = 3000, dsmin = 1e-7, ds = 1e-3, newton_options = NewtonPar(tol = 1e-9, verbose = true), detect_bifurcation = 3, tol_stability = 1e-4)

br_po = @time continuation(br, 4, opts_br_po,
                    PeriodicOrbitOCollProblem(50,5; jacobian = BK.DenseAnalytical(), meshadapt = true);
                    verbosity = 0,
                    alg = PALC(tangent = Bordered()),
                    linear_algo = COPBLS(),
                    normC = norminf,
                    callback_newton = BK.cbMaxNorm(1e2), #limit residual to avoid Inf or NaN
)

opts_br_po_2 = ContinuationPar(opts_br, p_min =-100.0, p_max = 400.0, max_steps = 3000, dsmax = 0.3, newton_options = NewtonPar(tol = 1e-9, verbose = false), detect_bifurcation = 3, tol_stability = 1e-4)

br_po_2 = @time continuation(br, 5, opts_br_po_2,
                    PeriodicOrbitOCollProblem(50,5; jacobian = BK.DenseAnalytical(), meshadapt = true); #We change 4 to 5
                    verbosity = 0,
                    alg = PALC(tangent = Bordered()),
                    linear_algo = COPBLS(),
                    normC = norminf,
                    callback_newton = BK.cbMaxNorm(1e2), #limit residual to avoid Inf or NaN
)

##Plotting everything
stability_po = [stable ? :green : :red for stable in br_po.stable]
stability_po_2 = [stable ? :green : :red for stable in br_po_2.stable]

plot(br)
plot!(br_po, vars=(:param, :max), c = stability_po, primary = false)
plot!(br_po, vars=(:param, :min), c = stability_po, primary = false)
plot!(br_po_2, vars=(:param, :max), c = stability_po_2, primary = false)
plot!(br_po_2, vars=(:param, :min), c = stability_po_2, primary = false)
