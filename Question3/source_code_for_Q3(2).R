## All of the procedure code necessary for Chapter 6 - Regression Testing

##############################################################################

## Ullah (1985) Test

##############################################################################

	## For OLS estimation
	## Do not include a column of ones (intercept included already)
	## If you want higher-order terms, you need to specify the full matrix in x.p
	## If you do not specify estimation procedure it will be LCLS
	## If you do not specify bandwidths, it will be LSCV
	ullah.ls <- function(y,x,x.p=NULL,np.est=NULL,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)

		if(is.null(x.p)==TRUE){
			x.p<-x
		}

		if(is.null(np.est)==TRUE){
			np.est = "lcls"
		}

		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)

  			if(np.est=="lcls") {
				bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="llls") {
				bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="lqls") {
				bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}

			bw <- bw.optim$par

		}

		## Parametric Estimation
		p.m <- lm(y~x.p)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m) 
 
		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Nonparametric estimation
  		if(np.est=="lcls") {
			np.model <- lcls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="llls") {
			np.model <- llls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="lqls") {
			np.model <- lqls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}

		## Calculate test statistic
		p.ssr <- t(p.resid)%*%p.resid
		np.ssr <- t(np.resid)%*%np.resid
		
		tstat <- (p.ssr-np.ssr)/np.ssr

		## Bootstrap test statistic
		boot.tstat <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			## Parametric Estimation
			p.m <- lm(y.star~x.p)
			p.resid.star <- residuals(p.m) 

			## Nonparametric estimation
  			if(np.est=="lcls") {
				np.model <- lcls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="llls") {
				np.model <- llls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="lqls") {
				np.model <- lqls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}

			## Calculate test statistic
			p.ssr <- t(p.resid.star)%*%p.resid.star
			np.ssr <- t(np.resid.star)%*%np.resid.star
		
			boot.tstat[jj] <- (p.ssr-np.ssr)/np.ssr

		}

		## Calculate p-value
		rank.tstat <- rank(c(tstat,boot.tstat))

		## p-values
		p.value <- 1-(rank.tstat[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat=tstat,p.value=p.value,boot.tstat=boot.tstat,bw=bw,nb=nb)

		return(return.list)

	}

	## For NLS estimates of a CES production function
	## Do not include a column of ones (intercept included already)
	## If you do not specify estimation procedure it will be LCLS
	## If you do not specify bandwidths, it will be LSCV
	ullah.ces <- function(y,x,np.est=NULL,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)

		x.name <- dimnames(x)[[2]]

		if(is.null(np.est)==TRUE){
			np.est = "lcls"
		}

		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)

  			if(np.est=="lcls") {
				bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="llls") {
				bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="lqls") {
				bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}

			bw <- bw.optim$par

		}

		## Parametric Estimation
		x.ces <- data.frame(y,x)
		p.m <- cesEst("y",x.name,x.ces,control=nls.lm.control(maxiter = 1000),vrs=TRUE)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m)  
 
		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Nonparametric estimation
  		if(np.est=="lcls") {
			np.model <- lcls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="llls") {
			np.model <- llls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="lqls") {
			np.model <- lqls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}

		## Calculate test statistic
		p.ssr <- t(p.resid)%*%p.resid
		np.ssr <- t(np.resid)%*%np.resid
		
		tstat <- (p.ssr-np.ssr)/np.ssr

		## Bootstrap test statistic
		boot.tstat <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			## Parametric Estimation
			x.ces <- data.frame(y.star,x)
			p.m <- cesEst("y.star",x.name,x.ces,control=nls.lm.control(maxiter = 1000),vrs=TRUE)
			p.resid.star <- residuals(p.m)  

			## Nonparametric estimation
  			if(np.est=="lcls") {
				np.model <- lcls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="llls") {
				np.model <- llls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="lqls") {
				np.model <- lqls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}

			## Calculate test statistic
			p.ssr <- t(p.resid.star)%*%p.resid.star
			np.ssr <- t(np.resid.star)%*%np.resid.star
		
			boot.tstat[jj] <- (p.ssr-np.ssr)/np.ssr

		}

		## Calculate p-value
		rank.tstat <- rank(c(tstat,boot.tstat))

		## p-values
		p.value <- 1-(rank.tstat[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat=tstat,p.value=p.value,boot.tstat=boot.tstat,bw=bw,nb=nb)

		return(return.list)

	}

	## NLS estimation
	## You must provide the formula (parametric function - pf)
	## You must provide the starting values (sv)
	## If you do not specify estimation procedure it will be LCLS
	## If you do not specify bandwidths, it will be LSCV
	ullah.nls <- function(y,x,pf,sv,np.est=NULL,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)

		if(is.null(np.est)==TRUE){
			np.est = "lcls"
		}
 
		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)

  			if(np.est=="lcls") {
				bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="llls") {
				bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}
  			if(np.est=="lqls") {
				bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
					control = list(npt = 4*n), y=y, x=x)
			}

			bw <- bw.optim$par

		}

		## Parametric Estimation
		p.m <- nls(pf,start=sv)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m) 
 
		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Nonparametric estimation
  		if(np.est=="lcls") {
			np.model <- lcls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="llls") {
			np.model <- llls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="lqls") {
			np.model <- lqls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}

		## Calculate test statistic
		p.ssr <- t(p.resid)%*%p.resid
		np.ssr <- t(np.resid)%*%np.resid
		
		tstat <- (p.ssr-np.ssr)/np.ssr

		## Bootstrap test statistic
		boot.tstat <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			pf[[2]] <- y.star

			## Parametric Estimation
			p.m <- nls(pf,start=sv)
			p.resid.star <- residuals(p.m) 

			## Nonparametric estimation
  			if(np.est=="lcls") {
				np.model <- lcls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="llls") {
				np.model <- llls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="lqls") {
				np.model <- lqls(y=y.star,x=x,bw)
				np.resid.star <- y.star - np.model[,1]
  			}

			## Calculate test statistic
			p.ssr <- t(p.resid.star)%*%p.resid.star
			np.ssr <- t(np.resid.star)%*%np.resid.star
		
			boot.tstat[jj] <- (p.ssr-np.ssr)/np.ssr

		}

		## Calculate p-value
		rank.tstat <- rank(c(tstat,boot.tstat))

		## p-values
		p.value <- 1-(rank.tstat[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat=tstat,p.value=p.value,boot.tstat=boot.tstat,bw=bw,nb=nb)

		return(return.list)

	}

##############################################################################

## Li and Wang (1998) Test

##############################################################################

	## For OLS estimation
	## Do not include a column of ones (intercept included already)
	## If you want higher-order terms, you need to specify the full matrix in x.p
	lw.ls <- function(y,x,x.p=NULL,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)

		if(is.null(x.p)==TRUE){
			x.p<-x
		}
 
		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)
			bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
				control = list(npt = 4*n), y=y, x=x)
			bw <- bw.optim$par

		}

		## Parametric Estimation
		p.m <- lm(y~x.p)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m) 
 
		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Calculate test statistic
		## Creating and storing the kernel matricies
		KK.store <- matrix(0,nrow=n,ncol=n)

		## starting the loop for the storage
		for (j in 1:n){

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x[j,jj])/bw[jj]

				KK <- KK*kk(dx)

			}

			KK.store[,j] <- KK
		
		}

		## creating product bandwidths
		bw.prod <- 1

		for (j in 1:q){

			bw.prod <- bw.prod*bw[j]
		
		}

		## Calculating the test statistics
		e.store <- numeric(n)
		e2.store <- numeric(n)
	
		for (j in 1:n){

			temp.store <- KK.store[j,]
			temp.store[j] <- 0

			e.store[j] <- sum(p.resid*temp.store)
			e2.store[j] <- sum(p.resid^2*temp.store^2)

		}

		## scaled and non-scaled test statitics
		tstat.un <- solve(n*(n-1)*bw.prod)*(sum(p.resid*e.store))
		sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid^2*e2.store))
		tstat.norm <- n*sqrt(bw.prod)*tstat.un/sqrt(sn)

		## Bootstrap test statistic
		boot.tstat.un <- numeric(nb)
		boot.tstat.norm <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			## Parametric Estimation
			p.m <- lm(y.star~x.p)
			p.resid.star <- residuals(p.m) 

			## Calculate the test statistic
			e.store <- numeric(n)
			e2.store <- numeric(n)
	
			for (j in 1:n){

				temp.store <- KK.store[j,]
				temp.store[j] <- 0

				e.store[j] <- sum(p.resid.star*temp.store)
				e2.store[j] <- sum(p.resid.star^2*temp.store^2)

			}
		
			boot.tstat.un[jj] <- solve(n*(n-1)*bw.prod)*(sum(p.resid.star*e.store))
			sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid.star^2*e2.store))
			boot.tstat.norm[jj] <- n*sqrt(bw.prod)*boot.tstat.un[jj]/sqrt(sn)

		}

		## Calculate p-value
		rank.tstat.un <- rank(c(tstat.un,boot.tstat.un))
		rank.tstat.norm <- rank(c(tstat.norm,boot.tstat.norm))

		## p-values
		p.value.un <- 1-(rank.tstat.un[1]/(nb+1))
		p.value.norm <- 1-(rank.tstat.norm[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat.un=tstat.un,tstat.norm=tstat.norm,p.value.un=p.value.un,
						p.value.norm=p.value.norm,boot.tstat.un=boot.tstat.un,
						boot.tstat.norm=boot.tstat.norm,bw=bw,nb=nb)

		return(return.list)

	}

	## For NLS estimates of a CES production function
	## Do not include a column of ones (intercept included already)
	lw.ces <- function(y,x,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)

		x.name <- dimnames(x)[[2]]
 
		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)
			bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
				control = list(npt = 4*n), y=y, x=x)
			bw <- bw.optim$par

		}

		## Parametric Estimation
		x.ces <- data.frame(y,x)
		p.m <- cesEst("y",x.name,x.ces,control=nls.lm.control(maxiter = 1000),vrs=TRUE)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m)  			

		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Calculate test statistic
		## Creating and storing the kernel matricies
		KK.store <- matrix(0,nrow=n,ncol=n)

		## starting the loop for the storage
		for (j in 1:n){

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x[j,jj])/bw[jj]

				KK <- KK*kk(dx)

			}

			KK.store[,j] <- KK
		
		}

		## creating product bandwidths
		bw.prod <- 1

		for (j in 1:q){

			bw.prod <- bw.prod*bw[j]
		
		}

		## Calculating the test statistics
		e.store <- numeric(n)
		e2.store <- numeric(n)
	
		for (j in 1:n){

			temp.store <- KK.store[j,]
			temp.store[j] <- 0

			e.store[j] <- sum(p.resid*temp.store)
			e2.store[j] <- sum(p.resid^2*temp.store^2)

		}

		## scaled and non-scaled test statitics
		tstat.un <- solve(n*(n-1)*bw.prod)*(sum(p.resid*e.store))
		sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid^2*e2.store))
		tstat.norm <- n*sqrt(bw.prod)*tstat.un/sqrt(sn)

		## Bootstrap test statistic
		boot.tstat.un <- numeric(nb)
		boot.tstat.norm <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			## Parametric Estimation
			x.ces <- data.frame(y.star,x)
			p.m <- cesEst("y.star",x.name,x.ces,control=nls.lm.control(maxiter = 1000),vrs=TRUE)
			p.resid.star <- residuals(p.m)  

			## Calculate the test statistic
			e.store <- numeric(n)
			e2.store <- numeric(n)
	
			for (j in 1:n){

				temp.store <- KK.store[j,]
				temp.store[j] <- 0

				e.store[j] <- sum(p.resid.star*temp.store)
				e2.store[j] <- sum(p.resid.star^2*temp.store^2)

			}
		
			boot.tstat.un[jj] <- solve(n*(n-1)*bw.prod)*(sum(p.resid.star*e.store))
			sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid.star^2*e2.store))
			boot.tstat.norm[jj] <- n*sqrt(bw.prod)*boot.tstat.un[jj]/sqrt(sn)

		}

		## Calculate p-value
		rank.tstat.un <- rank(c(tstat.un,boot.tstat.un))
		rank.tstat.norm <- rank(c(tstat.norm,boot.tstat.norm))

		## p-values
		p.value.un <- 1-(rank.tstat.un[1]/(nb+1))
		p.value.norm <- 1-(rank.tstat.norm[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat.un=tstat.un,tstat.norm=tstat.norm,p.value.un=p.value.un,
						p.value.norm=p.value.norm,boot.tstat.un=boot.tstat.un,
						boot.tstat.norm=boot.tstat.norm,bw=bw,nb=nb)

		return(return.list)

	}

	## NLS estimation
	## You must provide the formula (parametric function - pf)
	## You must provide the starting values (sv)
	lw.nls <- function(y,x,pf,sv,bw=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)
 
		## Estimate bandwidths if not given
		if(is.null(bw)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)
			bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
				control = list(npt = 4*n), y=y, x=x)
			bw <- bw.optim$par

		}

		## Parametric Estimation
		p.m <- nls(pf,start=sv)
		p.fitted <- fitted(p.m)
		p.resid <- residuals(p.m) 

		## Centered residuals
		p.resid.c <- p.resid - mean(p.resid)

		## Calculate test statistic
		## Creating and storing the kernel matricies
		KK.store <- matrix(0,nrow=n,ncol=n)

		## starting the loop for the storage
		for (j in 1:n){

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x[j,jj])/bw[jj]

				KK <- KK*kk(dx)

			}

			KK.store[,j] <- KK
		
		}

		## creating product bandwidths
		bw.prod <- 1

		for (j in 1:q){

			bw.prod <- bw.prod*bw[j]
		
		}

		## Calculating the test statistics
		e.store <- numeric(n)
		e2.store <- numeric(n)
	
		for (j in 1:n){

			temp.store <- KK.store[j,]
			temp.store[j] <- 0

			e.store[j] <- sum(p.resid*temp.store)
			e2.store[j] <- sum(p.resid^2*temp.store^2)

		}

		## scaled and non-scaled test statitics
		tstat.un <- solve(n*(n-1)*bw.prod)*(sum(p.resid*e.store))
		sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid^2*e2.store))
		tstat.norm <- n*sqrt(bw.prod)*tstat.un/sqrt(sn)

		## Bootstrap test statistic
		boot.tstat.un <- numeric(nb)
		boot.tstat.norm <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(p.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			p.resid.boot <- ifelse(unif <= prob.a, a*p.resid.c, b*p.resid.c)

    			y.star <- p.fitted + p.resid.boot

			pf[[2]] <- y.star

			## Parametric Estimation
			p.m <- nls(pf,start=sv)
			p.resid.star <- residuals(p.m) 
  
			## Calculate the test statistic
			e.store <- numeric(n)
			e2.store <- numeric(n)
	
			for (j in 1:n){

				temp.store <- KK.store[j,]
				temp.store[j] <- 0

				e.store[j] <- sum(p.resid.star*temp.store)
				e2.store[j] <- sum(p.resid.star^2*temp.store^2)

			}
		
			boot.tstat.un[jj] <- solve(n*(n-1)*bw.prod)*(sum(p.resid.star*e.store))
			sn <- 2*solve(n*(n-1)*bw.prod)*(sum(p.resid.star^2*e2.store))
			boot.tstat.norm[jj] <- n*sqrt(bw.prod)*boot.tstat.un[jj]/sqrt(sn)

		}

		## Calculate p-value
		rank.tstat.un <- rank(c(tstat.un,boot.tstat.un))
		rank.tstat.norm <- rank(c(tstat.norm,boot.tstat.norm))

		## p-values
		p.value.un <- 1-(rank.tstat.un[1]/(nb+1))
		p.value.norm <- 1-(rank.tstat.norm[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat.un=tstat.un,tstat.norm=tstat.norm,p.value.un=p.value.un,
						p.value.norm=p.value.norm,boot.tstat.un=boot.tstat.un,
						boot.tstat.norm=boot.tstat.norm,bw=bw,nb=nb)

		return(return.list)
 
 	}

##############################################################################

## Lavergne and Vuong (2000) Test

##############################################################################

	## Things people need to know
	lv.test <- function(y,x,w,np.est=NULL,np.bw=NULL,bw.x=NULL,bw.w=NULL,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		w <- as.matrix(w)	
		qw <- ncol(w)
		ones <- rep(1,n)
 
		if(is.null(np.est)==TRUE){
			np.est = "lcls"
		}

		if(is.null(np.bw)==TRUE){
			np.bw = "lscv"
		}

		## Estimate bandwidths for x if not given
		if(is.null(bw.x)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)

			if(np.bw=="lscv"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
			}

			if(np.bw=="aicc"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}

			}

			bw.x <- bw.optim$par

		}

		## Estimate bandwidths for w if not given
		if(is.null(bw.w)==TRUE){

			## Rule of Thumb
			bw.start <- numeric(qw)

			for (j in 1:qw){

				bw.start[j] <- 1.06*sd(w[,j])*n^(-1/(4+qw))

			}

			lower <- rep(0,qw)
			upper <- 5*sd(w)

			if(np.bw=="lscv"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}

			}

			if(np.bw=="aicc"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=w)
				}

			}

			bw.w <- bw.optim$par

		}

		## Restricted model estimation
  		if(np.est=="lcls") {
			np.model <- lcls(y=y,x=w,h=bw.w,loo=1)
			np.fitted <- np.model[,1]
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="llls") {
			np.model <- llls(y=y,x=w,h=bw.w,loo=1)
			np.fitted <- np.model[,1]
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="lqls") {
			np.model <- lqls(y=y,x=w,h=bw.w,loo=1)
			np.fitted <- np.model[,1]
			np.resid <- y - np.model[,1]
  		}

		## Centered residuals
		np.resid.c <- np.resid - mean(np.resid)

		## Calculate test statistic
		## Creating and storing the kernel matricies
		KK.store <- matrix(0,nrow=n,ncol=n)
		KK.store.w <- matrix(0,nrow=n,ncol=n)

		## starting the loop for the storage (x and w)
		for (j in 1:n){

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x[j,jj])/bw.x[jj]
				KK <- KK*kk(dx)

			}

			KK.store[,j] <- KK

			KK <- ones

			for (jj in 1:qw){

				dx <- (w[,jj] - w[j,jj])/bw.w[jj]
				KK <- KK*kk(dx)

			}

			KK.store.w[,j] <- KK
			KK.store.w[j,j] <- 0	# leave-one-out
		
		}

		## creating product bandwidths
		bw.prod <- 1
		bw.prod.w <- 1

		for (j in 1:q){
			bw.prod <- bw.prod*bw.x[j]
		}
		for (j in 1:qw){
			bw.prod.w <- bw.prod.w*bw.w[j]
		}

		## Estimating the density f_w (leave-one-out)
		fw <- solve((n-1)*bw.prod.w)*colSums(KK.store.w)

		## Calculating the test statistics
		e.store <- numeric(n)
		e2.store <- numeric(n)
	
		for (j in 1:n){

			temp.store <- KK.store[j,]
			temp.store[j] <- 0

			e.store[j] <- sum(np.resid*fw*temp.store)
			e2.store[j] <- sum(np.resid^2*fw^2*temp.store^2)

		}

		## scaled and non-scaled test statitics
		tstat.un <- solve(n*(n-1)*bw.prod)*(sum(np.resid*fw*e.store))
		sn <- 2*solve(n*(n-1)*bw.prod)*(sum(np.resid^2*fw^2*e2.store))
		tstat.norm <- n*sqrt(bw.prod)*tstat.un/sqrt(sn)

		## Bootstrap test statistic
		boot.tstat.un <- numeric(nb)
		boot.tstat.norm <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(np.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			np.resid.boot <- ifelse(unif <= prob.a, a*np.resid.c, b*np.resid.c)

    			y.star <- np.fitted + np.resid.boot

			## Restricted model estimation
  			if(np.est=="lcls") {
				np.model <- lcls(y=y.star,x=w,h=bw.w,loo=1)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="llls") {
				np.model <- llls(y=y.star,x=w,h=bw.w,loo=1)
				np.resid.star <- y.star - np.model[,1]
  			}
  			if(np.est=="lqls") {
				np.model <- lqls(y=y.star,x=w,h=bw.w,loo=1)
				np.resid.star <- y.star - np.model[,1]
  			}
  
			## Calculate the test statistic
			e.store <- numeric(n)
			e2.store <- numeric(n)
	
			for (j in 1:n){

				temp.store <- KK.store[j,]
				temp.store[j] <- 0

				e.store[j] <- sum(np.resid.star*fw*temp.store)
				e2.store[j] <- sum(np.resid.star^2*fw^2*temp.store^2)

			}
		
			boot.tstat.un[jj] <- solve(n*(n-1)*bw.prod)*(sum(np.resid.star*fw*e.store))
			sn <- 2*solve(n*(n-1)*bw.prod)*(sum(np.resid.star^2*fw^2*e2.store))
			boot.tstat.norm[jj] <- n*sqrt(bw.prod)*boot.tstat.un[jj]/sqrt(sn)

		}

		## Calculate p-value
		rank.tstat.un <- rank(c(tstat.un,boot.tstat.un))
		rank.tstat.norm <- rank(c(tstat.norm,boot.tstat.norm))

		## p-values
		p.value.un <- 1-(rank.tstat.un[1]/(nb+1))
		p.value.norm <- 1-(rank.tstat.norm[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat.un=tstat.un,tstat.norm=tstat.norm,p.value.un=p.value.un,
						p.value.norm=p.value.norm,boot.tstat.un=boot.tstat.un,
						boot.tstat.norm=boot.tstat.norm,bw.x=bw.x,bw.w=bw.w,nb=nb)

		return(return.list)
 
 	}

##############################################################################

## Zheng (2009) Test

##############################################################################

	## If you do not specify np.est and np.bw, they will be lcls and lscv, respectively
	## If you do not specify bw, it will calculate bandwidths via lcls and lscv
	zheng.test <- function(y,x,np.est=NULL,np.bw=NULL,bw,nb){

		n <- length(y)
		x <- as.matrix(x)
		q <- ncol(x)
		ones <- rep(1,n)
 
		if(is.null(np.est)==TRUE){
			np.est = "lcls"
		}

		if(is.null(np.bw)==TRUE){
			np.est = "lscv"
		}
 
		## Estimate bandwidths for x if not given
		if(is.null(bw)==TRUE){ 

			## Rule of Thumb
			bw.start <- numeric(q)

			for (j in 1:q){

				bw.start[j] <- 1.06*sd(x[,j])*n^(-1/(4+q))

			}

			lower <- rep(0,q)
			upper <- 5*sd(x)

			if(np.bw=="lscv"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.lscv, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
			}

			if(np.bw=="aicc"){

  				if(np.est=="lcls") {
					bw.optim <- bobyqa(bw.start,lcls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="llls") {
					bw.optim <- bobyqa(bw.start,llls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}
  				if(np.est=="lqls") {
					bw.optim <- bobyqa(bw.start,lqls.aicc, lower, upper, 
						control = list(npt = 4*n), y=y, x=x)
				}

			}

			bw <- bw.optim$par

		}

		## Nonparametric estimation
  		if(np.est=="lcls") {
			np.model <- lcls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="llls") {
			np.model <- llls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}
  		if(np.est=="lqls") {
			np.model <- lqls(y=y,x=x,bw)
			np.resid <- y - np.model[,1]
  		}

		## Homoskedastic variance
		sn.homo <- solve(n)*t(np.resid)%*%np.resid

		## difference between squared residual and homoskedastic variance
		resid2.homo <- np.resid^2 - sn.homo

		## Creating and storing the kernel matricies
		KK.store <- matrix(0,nrow=n,ncol=n)

		## starting the loop for the storage
		for (j in 1:n){

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x[j,jj])/bw[jj]
				KK <- KK*kk(dx)

			}

			KK.store[,j] <- KK
		
		}

		## creating product bandwidths
		bw.prod <- 1

		for (j in 1:q){
			bw.prod <- bw.prod*bw[j]
		}

		## Calculating the test statistics
		e.store <- numeric(n)
		e2.store <- numeric(n)
	
		for (j in 1:n){

			temp.store <- KK.store[j,]
			temp.store[j] <- 0

			e.store[j] <- sum(resid2.homo*temp.store)
			e2.store[j] <- sum(resid2.homo^2*temp.store^2)

		}

		## scaled and non-scaled test statitics
		tstat.un <- solve(n*(n-1)*bw.prod)*(sum(resid2.homo*e.store))
		sn <- 2*solve(n*(n-1)*bw.prod)*(sum(resid2.homo^2*e2.store))
		tstat.norm <- n*sqrt(bw.prod)*tstat.un/sqrt(sn)

		## Bootstrap test statistic
		boot.tstat.un <- numeric(nb)
		boot.tstat.norm <- numeric(nb)

		for (jj in 1:nb){

    			## Wild bootstrap
			unif <- runif(length(np.resid),0,1)

			a <- (1-sqrt(5))/2
			prob.a <- (sqrt(5)+1)/(2*sqrt(5))
			b <- (1+sqrt(5))/2

			resid2.homo.star <- ifelse(unif <= prob.a, a*resid2.homo, b*resid2.homo)
  
			## Calculate the test statistic
			e.store <- numeric(n)
			e2.store <- numeric(n)
	
			for (j in 1:n){

				temp.store <- KK.store[j,]
				temp.store[j] <- 0

				e.store[j] <- sum(resid2.homo.star*temp.store)
				e2.store[j] <- sum(resid2.homo.star^2*temp.store^2)

			}
		
			boot.tstat.un[jj] <- solve(n*(n-1)*bw.prod)*(sum(resid2.homo.star*e.store))
			sn <- 2*solve(n*(n-1)*bw.prod)*(sum(resid2.homo.star^2*e2.store))
			boot.tstat.norm[jj] <- n*sqrt(bw.prod)*boot.tstat.un[jj]/sqrt(sn)

		}

		## Calculate p-value
		rank.tstat.un <- rank(c(tstat.un,boot.tstat.un))
		rank.tstat.norm <- rank(c(tstat.norm,boot.tstat.norm))

		## p-values
		p.value.un <- 1-(rank.tstat.un[1]/(nb+1))
		p.value.norm <- 1-(rank.tstat.norm[1]/(nb+1))

		## Listing items to send back
		return.list <- list(tstat.un=tstat.un,tstat.norm=tstat.norm,p.value.un=p.value.un,
						p.value.norm=p.value.norm,boot.tstat.un=boot.tstat.un,
						boot.tstat.norm=boot.tstat.norm,bw=bw,nb=nb)

		return(return.list)
 
 	}





