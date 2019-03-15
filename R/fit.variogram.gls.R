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
	if (model$model[1] == 'Nug') {
		if (ignoreInitial)
			init = c(mean(y)/2, mean(y)/2, median(h0)/4)
		else
			init = c(model$psill, model$range[2])
		gamfn0 = function(h, th, m = as.character(model$model[2]))
	    	variogramLine(vgm(th[2], m, th[3], th[1]), dist_vector=h)$gamma
		minfuncols0 = function(theta=rbind(1,1,1)) {
			res = y-gamfn0(h0,theta)
			sum(res^2)
		}
	
		# If variogram is Power model, maximum range should be 2. 
		# Also if ignoreInitial change initial value of range to 1
		if (any(model$model == 'Pow')) { 
	  		upperOptim <- c(max(y), max(y), 2)
	  		if (ignoreInitial)
				init[3] <- 1
		} else 
			upperOptim = c(max(y), max(y), max(h0)) 
	
		th = th0 = th.ols = optim(init, minfuncols0, gr=NULL, method="L-BFGS-B",
			lower=c(0,1e-9,1e-9), upper=upperOptim)$par
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
				res = y-gamfn0(h0,c(th[1],th[2],range))
				t(res) %*% solve(cov, res)
			}
			minfuncsill = function(sill) {
				res = y-gamfn0(h0, c(th[1],sill,th[3]))
				t(res) %*% solve(cov, res)
			}
			minfuncnugget = function(nugget) {
				res = y - gamfn0(h0, c(nugget, th[2], th[3]))
				t(res) %*% solve(cov, res)
			}
			th0 = th
			# Avoid calculating max every time
			th[1] = optimize(minfuncnugget,lower=0,upper=upperOptim[1])$minimum
			th[2] = optimize(minfuncsill,lower=1e-9,upper=upperOptim[2])$minimum
			th[3] = optimize(minfuncrange,lower=1e-9,upper=upperOptim[3])$minimum
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
	} else {
	  # No Nugget component. Use only psill = th[1] and range = th[2] parameters
	  if (ignoreInitial)
	    init = c(mean(y)/2, median(h0)/4)
	  else
	    init = c(model$psill, model$range)
	  gamfn = function(h, th, m = as.character(model$model))
	    variogramLine(vgm(psill = th[1], model = m, range = th[2]), dist_vector=h)$gamma
	  minfuncols = function(theta=rbind(1,1)) {
	    res = y-gamfn(h0,theta)
	    sum(res^2)
	  }
	  
	  # If variogram is Power model, maximum range should be 2. 
	  # Also if ignoreInitial change initial value of range to 1
	  if (any(model$model == 'Pow')) { 
	    upperOptim <- c(max(y), 2)
	    if (ignoreInitial) init[2] <- 1
	  } else { upperOptim = c(max(y), max(h0)) }
	  
	  th = th0 = th.ols = optim(init, minfuncols,gr=NULL, method="L-BFGS-B",
	                            lower=c(1e-9,1e-9), upper=upperOptim)$par
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
	      res = y-gamfn(h0,c(th[1],range))
	      t(res) %*% solve(cov, res)
	    }
	    minfuncsill = function(sill) {
	      res = y-gamfn(h0, c(sill,th[2]))
	      t(res) %*% solve(cov, res)
	    }
	    th0 = th
	    th[1] = optimize(minfuncsill,lower=1e-9,upper=upperOptim[1])$minimum
	    th[2] = optimize(minfuncrange,lower=1e-9,upper=upperOptim[2])$minimum
	    converged = sum(abs((th - th0)/th0)) < eps
	    iter = iter + 1
	    if (trace)
	      print(th)
	    model$psill = th[1]
	    model$range = th[2]
	  }
	  if (th[2] / max(h0) > .99)
	    warning("range parameter at search space boundary")
	  if (!converged) {
	    warning("no convergence, returning OLS solution")
	    th = th.ols
	    model$psill = th[1]
	    model$range = th[2]
	  }
	}
	
	if (plot)
		plot(variogram(formula, data, cloud = TRUE, cutoff = cutoff), 
			model = model)
	else
		model
}

