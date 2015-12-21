fit.variogram.gls <-
function(formula, data, model, maxiter = 30, 
		eps = .01, trace = TRUE, ignoreInitial = TRUE, cutoff = Inf,
		plot = FALSE) {
	v = as.data.frame(variogram(formula, data, cloud = TRUE, cutoff = cutoff))
	i = v$left
	j = v$right
	y = v$gamma
	h0 = v$dist
	dists = spDists(data)
	n = length(i)
	iter = 0
	converged = FALSE
	if (ignoreInitial)
		init = c(mean(y)/2, mean(y)/2, median(h0)/4)
	else
		init = c(model$psill, model$range[2])
	gamfn = function(h, th, m = as.character(model$model[2]))
	    variogramLine(vgm(th[2], m, th[3], th[1]), dist_vector=h)$gamma
	minfuncols = function(theta=rbind(1,1,1)) {
		res = y-gamfn(h0,theta)
		sum(res^2)
	}
	th = th0 = th.ols = optim(init, minfuncols,gr=NULL, method="L-BFGS-B",
		lower=c(0,1e-9,1e-9), upper=c(max(y),max(y),max(h0)))$par
	if (trace)
		print(th)
	while (!converged && iter < maxiter) {
		comb = function(i,j) cbind(rep(i,length(j)), rep(j,each=length(i)))
		cov = matrix(
			variogramLine(model, dist_vector = dists[comb(i,i)])$gamma
			+ variogramLine(model, dist_vector = dists[comb(j,j)])$gamma
			- variogramLine(model, dist_vector = dists[comb(i,j)])$gamma
			- variogramLine(model, dist_vector = dists[comb(j,i)])$gamma,
			length(j), length(j))
		cov = 0.5 * cov ^ 2
		#cov = solve(cov)
		cov = qr(cov)
		minfuncrange = function(range) {
			res = y-gamfn(h0,c(th[1],th[2],range))
			t(res) %*% solve(cov, res)
		}
		minfuncsill = function(sill) {
			res = y-gamfn(h0, c(th[1],sill,th[3]))
			t(res) %*% solve(cov, res)
		}
		minfuncnugget = function(nugget) {
			res = y - gamfn(h0, c(nugget, th[2], th[3]))
			t(res) %*% solve(cov, res)
		}
		th0 = th
		th[1] = optimize(minfuncnugget,lower=0,upper=max(y))$minimum
		th[2] = optimize(minfuncsill,lower=1e-9,upper=max(y))$minimum
		th[3] = optimize(minfuncrange,lower=1e-9,upper=max(h0))$minimum
		converged = sum(abs((th - th0)/th0)) < eps
		iter = iter + 1
		if (trace)
			print(th)
		model$psill = c(th[1], th[2])
		model$range[2] = th[3]
	}
	if (th[3] / max(h0) > .99)
		warning("range parameter at search space boundary")
	if (!converged) {
		warning("no convergence, returning OLS solution")
		th = th.ols
		model$psill = c(th[1], th[2])
		model$range[2] = th[3]
	}
	if (plot)
		plot(variogram(formula, data, cloud = TRUE, cutoff = cutoff), 
			model = model)
	else
		model
}

