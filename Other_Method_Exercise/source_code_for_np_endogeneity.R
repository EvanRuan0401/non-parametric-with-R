## Source code for Chapter 10 of Henderson and Parmeter
## Procedure code necessary for local polynomial least squares estimation
## and bandwidth selection

require(minqa)
require(np)

## Function to avoid division by zero
NZD <- function(a){
  sapply(1:NROW(a),function(i) {if(a[i]<0) min(-.Machine$double.xmin,a[i]) else max(.Machine$double.xmin,a[i])})
}

## Function to calculate bounds on first stage bandwidths for 
## Su and Ullah Estimator.
SUh <- function(p1,p2,d1,d2=1,dx=1){

	optimal <- 2*(p1+1)+d1+d2
	gamma   <- 2*(p2+1)+dx+d1
	lower   <- max(c((p2+1)/(p1+1),(p2+3)/(2*(p1+1))))
	upper   <- (p2+dx+d1-1)/(d1+d2)

	return(list(opt=1/optimal,upper=upper,lower=lower,gamma=gamma))

}

###############################################################################
## Local Polynomial Least Squares 
###############################################################################

lpls <- function(y,x.eval=NULL,x.train,h,p,loo=FALSE,kr="gaussian"){

	## Call lcls if p=0
	if(p==0){return(lcls1(y=y,x.eval=x.eval,x.train=x.train,h=h,loo=loo,kr=kr))}

	## Set x.eval equal to x.train if x.eval is NULL
	if(is.null(x.eval)) x.eval <- x.train

	X.train <- as.data.frame(x.train)
	X.eval  <- as.data.frame(x.eval)

	n.train <- nrow(X.train)
	n.eval  <- nrow(X.eval)

	X.col.numeric <- sapply(1:ncol(X.train),function(i){is.numeric(X.train[,i])})

  	## k represents the number of numeric regressors, this will return
  	## zero if there are none

  	q <- ncol(as.data.frame(X.train[,X.col.numeric]))
	q1 <- q+1
	
	## Number of terms in the local polynomial estimate
	## Determine the number of cross-product terms. It depends on the
  	## number of numeric regressors and on the order of the
  	## polynomial. The cross products can max out when p > q, p must be
  	## greater than one.
	pq1 <- p*q+1
	num.cp <- 0
	
	if(q > 0){
		X.train.numeric <- as.data.frame(X.train[,X.col.numeric])
		X.eval.numeric  <- as.data.frame(X.eval[,X.col.numeric])
   		if(p>=2 && q>=2){for(i in 2:min(q,p)) num.cp <- num.cp + ncol(combn(1:q,i))}
   	}
	
	ones <- rep(1,n.train)
	epsilon <- 0.00001

	XX <- matrix(1,n.train,pq1+num.cp)
	thetahat <- matrix(0,ncol=ncol(XX),nrow=0)

	for(j in 1:n.eval){ 

		KK <- ones

		## First construct W matrix for lpls estimation
		## differences, then squares, then cubic, then 
		## cross products

		if(q > 0 && p > 0) {

			l.beg <- 1

			## Fill in polynomial terms
			for (l in 1:p){

				XX[,((l-1)*q+2):(l*q+1)] <- (1/factorial(l))*(as.matrix(X.train.numeric)-as.matrix(X.eval.numeric[rep(j,n.train),]))^l

				## Fill in cross-product terms -- this requires q>1

				if(q >= 2 && l >= 2 && l <= q){

          			l.end <- l.beg + ncol(as.matrix(combn(1:q,l)))-1
      
          			## Initialize the cross-product terms if l <= k
          			## ncol(combn(1:k,l))) columns taken which shrinks as l
          			## increases.

					col.use  <- combn(1:q,l)[1,]+1
					col.stor <- (pq1+l.beg):(pq1+l.end)

          			XX[,col.stor] <- as.matrix(XX[,col.use])

          			## Generate cross-products
          			for(jj in 2:l){

						col.jj <- combn(1:q,l)[jj,]+1

          				XX[,col.stor] <- XX[,col.stor]*as.matrix(XX[,col.jj])
          			}

          			l.beg <- l.end + 1

				}
			}

		}
		
		KK <- npksum(exdat=X.train,
					 txdat=X.eval[j,],
					 bws=h,
					 ckertype=kr,
					 ukertype="liracine",
					 okertype="liracine")$ksum

		if(loo){KK[j]<- 0}

		XK <- XX
		
		for (jj in 1:ncol(XK)){XK[,jj] <- XX[,jj]*KK}

		XK <- t(XK)

		## Ridging (if necessary)
		ridge <- 0

		qx1 <- nrow(XK)

		while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,qx1)))),
				error = function(e){
				return(matrix(FALSE,qx1,qx1))
				})[1,1]==FALSE) {
			ridge <- ridge + epsilon
		}

		# warning.immediate <- FALSE
		# warning(paste("Ridging obs. ", j, ", ridge = ", 
						# signif(ridge,6),sep=""),
				# immediate.=warning.immediate,
				# call.=!warning.immediate)

		XKX <- XK%*%XX + diag(rep(ridge,qx1))
		XKXI <- solve(XKX)

		XKY <- XK%*%y

      	thetahat.estimate <- XKXI%*%XKY	
		thetahat.estimate <- t(thetahat.estimate)
		thetahat <- rbind(thetahat,thetahat.estimate)

	}

	return(thetahat)

}

## LCLS (loo = TRUE -> leave-one-out estimator)
lcls1 <- function(y,x.train,x.eval,h,loo=FALSE,kr="gaussian"){

	## Set x.eval equal to x.train if x.eval is NULL
	if(is.null(x.eval)) x.eval <- x.train

	X.train <- as.data.frame(x.train)
	X.eval  <- as.data.frame(x.eval)

	n.train <- nrow(X.train)
	n.eval  <- nrow(X.eval)
	
	X.col.numeric <- sapply(1:ncol(X.train),function(i){is.numeric(X.train[,i])})
	id <- which(X.col.numeric)
	
  	## q represents the number of numeric regressors, this will return
  	## zero if there are none

  	q <- ncol(as.data.frame(X.train[,X.col.numeric]))
	q1 <- q+1

	ones <- rep(1,n.train)

	thetahat <- matrix(0,ncol=q1,nrow=n.eval)

	for(j in 1:n.eval){ 

		KK <- ones
    	DKK <- ones

		KK <- npksum(exdat=X.train,
					 txdat=X.eval[j,],
					 bws=h,
					 ckertype=kr,
					 ukertype="liracine",
					 okertype="liracine")$ksum
			
		if(loo){KK[j]<- 0}

		thetahat[j,1] <- sum(y*KK)/NZD(sum(KK))

		operator <- rep("normal",ncol(X.train))

		## Now calculate derivatives for continuous covariates
		if(q>0){

			for(jj in 1:q){
    			
    			operator <- rep("normal",ncol(X.train))
    			operator[id[jj]] <- "derivative"
    			dkdx1 <- -npksum(exdat=X.train,
					 			  txdat=X.eval[j,],
					 			  bws=h,
					 			  ckertype=kr,
					 			  ukertype="liracine",
					 			  okertype="liracine",
					 			  operator=operator)$ksum/h[jj]
				
				if(loo){dkdx1[j]<- 0}

				thetahat[j,jj+1] <- ((sum(y*dkdx1)*sum(KK))-(sum(y*KK)*sum(dkdx1)))/(NZD((sum(KK))^2))

			}
		}
	}

	return(thetahat)

}

############################################################################
## Least Squares Cross Validation Bandwidth Selection
############################################################################
lscv.fun <- function(h,y,x,poly=0,kr="gaussian"){
 
	if(min(h) <= 0)return(.Machine$double.xmax)

	mhat <- lpls(y=y,x.train=x,h=h,p=poly,loo=TRUE,kr=kr)[,1]
 	
	return(mean((y - mhat)^2))
	
}

bw.lscv <- function(y,x,poly=0,kr="gaussian"){

	x <- as.data.frame(x)
	x.numeric <- sapply(1:ncol(x),function(i){is.numeric(x[,i])})

  	## First initialize initial search values of the vector of
  	## bandwidths to lie in [0,1]
  
  	bw.start <- runif(ncol(x))

  	## Next, figure out which are of type numeric and use initial values
  	## that lie in a draw from [.5, 1.5] times \sigma n^{-1/(4+k)}

  	k <- length(which(x.numeric))
  	lower <- rep(0,ncol(x))
	upper <- rep(1,ncol(x))

  	for(i in 1:ncol(x)) {
    	if(x.numeric[i]==TRUE) {
      		bw.start[i] <- runif(1,.5,1.5)*(IQR(x[,i])/1.349)*nrow(x)^{-1/(4+k)}
      		upper[i]    <- 3*IQR(x[,i])   ## Can switch to this once sims finish to .Machine$double.xmax
    	}
  	}

	bw.optim <- bobyqa(bw.start,lscv.fun, lower, upper, 
						control=list(npt=4*n,maxfun=100), 
						y=y,x=x,poly=poly,kr=kr)

	return(list(bw=bw.optim$par,fun=bw.optim$fval))

}

############################################################################
## Gradient Based Cross Validation, 
############################################################################
##CURRENTLY SETUP FOR VECTOR X NEED TO MAKE THIS FOR MATRIX X

gbcv.fun <- function(h,x,y,poly=0,loo.ll=TRUE,loo.lp=FALSE,
						kr="gaussian",beta.tr=NULL,id=NULL,
						trim=FALSE,alpha=0.05,freq=FALSE){

	x <- as.data.frame(x)

	x.numeric <- sapply(1:ncol(x),function(i){is.numeric(x[,i])})
	qid <- which(x.numeric)
	qq  <- length(qid)

	if(freq){
		
		if(any((h[which(!x.numeric)])>0)){return(.Machine$double.xmax)}
		
	}

	if(is.null(id)){id <- 2:(qq+1)}

	if(min(h) <= .Machine$double.xmin)	return(.Machine$double.xmax)

	beta.ll <- lpls(y=y,x=x,p=1,h=h,loo=loo.ll,kr=kr)[,id]

	if(is.null(beta.tr)){
			beta.tr <- lpls(y=y,x=x,h=h,p=poly,loo=loo.lp,kr=kr)[,id]
	}

	M <- 1

	if(trim && qq>0){

		M <- matrix(1,nrow(x),qq)

		for(j in 1:qq){
				M[,j] <- ifelse(abs(x[,qid[j]])<=quantile(abs(x[,qid[j]]),1-alpha),1,0)

		}

		M <- apply(M,1,FUN=prod)

	}
	
	A <- M*(beta.tr-beta.ll)^2
	
	if(ncol(as.data.frame(A))>1){A<-rowSums(A)}

	CV <- mean(A)
	return(CV)

}

## This function allows for a range of gradient-based cross
## validation procedures. One can pass in gradients ahead of
## time, via beta.tr, or compare to gradients from a local
## estimator via poly.

bw.gbcv <- function(y,x,poly=NULL,kr="gaussian",
					loo.ll=TRUE,loo.lp=FALSE,
					beta.tr=NULL,id=NULL,
					trim=FALSE,alpha=0.05,freq=FALSE){

	x <- as.data.frame(x)
	x.numeric <- sapply(1:ncol(x),function(i){is.numeric(x[,i])})

  	## First initialize initial search values of the vector of
  	## bandwidths to lie in [0,1]
  
  	bw.start <- runif(ncol(x))

  	## Next, figure out which are of type numeric and use initial values
  	## that lie in a draw from [.5, 1.5] times \sigma n^{-1/(4+k)}

  	k <- length(which(x.numeric))
  	lower <- rep(0,ncol(x))
	upper <- rep(1,ncol(x))

  	for(i in 1:ncol(x)) {
    	if(x.numeric[i]==TRUE) {
      		bw.start[i] <- runif(1,.5,1.5)*(IQR(x[,i])/1.349)*nrow(x)^{-1/(4+k)}
      		upper[i]    <- 3*IQR(x[,i])   ## Can switch to this once sims finish to .Machine$double.xmax
    	}
  	}

	## Finding the minimum for the lscv function

	bw.optim <- bobyqa(bw.start,gbcv.fun,lower,upper, 
						control=list(maxfun=250), 
						y=y,x=x,poly=poly,kr=kr,
						beta.tr=beta.tr,id=id,
						loo.ll=loo.ll,loo.lp=loo.lp,
						trim=trim,alpha=alpha,freq=freq)

	return(list(bw=bw.optim$par,fun=bw.optim$fval,num.eval=bw.optim$feval))

}


lscv.fun <- function(h,y,x,poly=0,kr="gaussian"){
 
	if(min(h) <= 0)return(.Machine$double.xmax)

	mhat <- lpls(y=y,x.train=x,h=h,p=poly,loo=TRUE,kr=kr)[,1]
 	
	return(mean((y - mhat)^2))
	
}


############################################################################
## Gradient Based Cross Validation, 
############################################################################
## Currently coded for a single endogenous variable
## Also assumes that we have exact identification, i.e. d_x=d_2

bw.lscv.su <- function(y,x,z1,z2,p1=3,p2=1,c=1,kr="gaussian"){

	xx <- data.frame(x,z1)
	xx.numeric <- sapply(1:ncol(xx),function(i){is.numeric(xx[,i])})

  	## First initialize initial search values of the vector of
  	## bandwidths to lie in [0,1]
  
  	bw.start <- runif(ncol(xx))

  	## Next, figure out which are of type numeric and use initial values
  	## that lie in a draw from [.5, 1.5] times \sigma n^{-1/(4+k)}

  	k <- length(which(xx.numeric))
  	lower <- rep(0,ncol(xx))
	upper <- rep(1,ncol(xx))

  	for(i in 1:ncol(xx)) {
    	if(xx.numeric[i]==TRUE) {
      		bw.start[i] <- runif(1,.5,1.5)*(IQR(xx[,i])/1.349)*nrow(xx)^{-1/(4+k)}
      		upper[i]    <- 3*IQR(xx[,i])
    	}
  	}
  	
  	n <- length(y)
  	id <- sample(1:n,ceiling(c*n))

	bw.optim <- bobyqa(bw.start,lscv.fun.su, lower, upper, 
						control=list(npt=4*n,maxfun=100), 
						y=y,x=x,z1=z1,z2=z2,p1=p1,p2=p2,id=id,kr=kr)

	return(list(bw=bw.optim$par,fun=bw.optim$fval))

}

lscv.fun.su <- function(h2,y,x,z1,z2,p1=3,p2=1,id=NULL,kr="gaussian"){

	d1 <- dim(as.data.frame(z1))[2]
	d2 <- dim(as.data.frame(z2))[2]

	n <- length(y)

	su2 <- SUh(p1=p1,p2=p2,d1=1,d2=1,dx=1)

	h1 <- h2*n^(-(su2$upper-su2$lower)/(2*su2$gamma))

	if(min(h1) <= 0)return(.Machine$double.xmax)
	if(min(h2) <= 0)return(.Machine$double.xmax)

	z <- data.frame(z1,z2)

	## Step 1
	ghat <- lpls(y=x,x.train=z,h=h1,p=p1,loo=TRUE,kr=kr)[,1]

	## Residuals for control function
	uhat <- x-ghat
	
	xx  <- data.frame(x)
	zz1 <- data.frame(z1)

	xx.train <- data.frame(xx,zz1,uhat)

	mhat <- as.numeric()

	for (j in 1:dim(xx.train)[1]){

		xx.eval <- expand.grid(xx=xx[j,],zz1=zz1[j,],sort(uhat)[id])
		mhat.cf <- lpls(y=y,x.train=xx.train,x.eval=xx.eval,
						h=c(h2,h2[1]),p=p2,kr=kr,loo=TRUE)[,1]

		mhat[j] <- mean(mhat.cf)

	}

	return(mean((y - mhat)^2))
	
}

