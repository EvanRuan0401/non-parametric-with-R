## All of the procedure code necessary for Chapter 5 - Regression 

##############################################################################

## Kernel functions

##############################################################################

  	#if(kernel=="epanechnikov") {#
    		#kk <- function(g) ifelse(abs(g)<=1, 0.75*(1-g^2), 0)#
  	#}#
  	#if(kernel=="biweight") {#
    		#kk <- function(g) ifelse(abs(g)<=1, (15/16)*(1-g^2)^2, 0)#
  	#}#
	#if(kernel=="triweight") {#
    		#kk <- function(g) ifelse(abs(g)<=1, (35/32)*(1-g^2)^3, 0)#
  	#}#
  	#if(kernel=="gaussian") {#
    		kk <- function(g) (1/sqrt(2*pi))*exp(-0.5*g^2)
  	#}#

###############################################################################

## Regression types (LCLS, LLLS, LQLS)

###############################################################################

	## LCLS (loo = 0 -> leave-one-out estimator)
	## loo = 0 by default and you only need to specify if you want a loo
	lcls <- function(y,x,h,loo=0){

		n <- length(y)
		q <- ncol(x)
		q1 <- q+1
		ones <- rep(1,n)
		ones1 <- rep(1,n-1)
		epsilon <- 0.00001
    
    		h.prod <- 1
    
    		for(jj in 1:q){

      		h.prod <- h.prod*h[jj]

    		}
      
		thetahat <- matrix(0,ncol=q1,nrow=n)

		for(j in 1:n){ 

 			if (loo==1){

				KK <- ones1
 
				for (jj in 1:q){
       
					dx <- (x[-j,jj] - x[j,jj])/h[jj]
					id <- which(abs(dx)<=1)
					KK <- KK*kk(dx)

				}

				thetahat[j,1] <- sum(y[-j]*KK)/sum(KK)

 			} else {

				KK <- ones

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					KK <- KK*kk(dx)
				}

				thetahat[j,1] <- sum(y*KK)/sum(KK)

				dkdx <- matrix(0,ncol=q,nrow=n)

				for (jj in 1:q){

					dx2 <- (x[,jj] - x[j,jj])/(h[jj]^2)
					dkdx[,jj] <- KK*dx2
				}

			}


			## Do not do gradients for loo because we only use it for bandwidth selection
			## Gradients are only correctly coded for Gaussian kernel
			if (loo==0){

				for (jj in 1:q){

					thetahat[j,jj+1] <- ((sum(y*dkdx[,jj])*sum(KK)) - (sum(y*KK)*sum(dkdx[,jj])))/((sum(KK))^2)

				}

			}


		}

		return(thetahat)

	}

	## LLLS (loo = 0 -> leave-one-out estimator)
	## loo = 0 by default and you only need to specify if you want a loo
	llls <- function(y,x,h,loo=0){

		n <- length(y)
		q <- ncol(x)
		q1 <- q+1
		ones <- rep(1,n)
		ones1 <- rep(1,n-1)
		epsilon <- 0.00001

		thetahat <- matrix(0,ncol=q1,nrow=0)

		for(j in 1:n){ 

 			if (loo==1){

				XX <- ones1
				KK <- ones1
 
				for (jj in 1:q){
       
					dx <- (x[-j,jj] - x[j,jj])/h[jj]
					KK <- KK*kk(dx)

					xt <- x[-j,jj] - x[j,jj]
					XX <- cbind(XX,xt)
				}

				XK <- ones1*KK

				for (jj in 1:q){

					XK <-cbind(XK,XX[,jj+1]*KK)

				}

 
 			} else {

				XX <- ones
				KK <- ones

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					xt <- x[,jj] - x[j,jj]
					XX <- cbind(XX,xt)
	
					KK <- KK*kk(dx)
				}

				XK <- ones*KK

				for (jj in 1:q){

					XK <-cbind(XK,XX[,jj+1]*KK)

				}

			}

			XK <- t(XK)

			## Ridging (if necessary)
			ridge <- 0

			while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,q1)))),
					error = function(e){
					return(matrix(FALSE,q1,q1))
					})[1,1]==FALSE) {
				ridge <- ridge + epsilon
			}

			XKX <- XK%*%XX + diag(rep(ridge,q1))
			XKXI <- solve(XKX)

			if (loo==1){

				XKY <- XK%*%y[-j]

			} else {

				XKY <- XK%*%y

			}

 	      	thetahat.estimate <- XKXI%*%XKY	
			thetahat.estimate <- t(thetahat.estimate)
			thetahat <- rbind(thetahat,thetahat.estimate)
		}

		return(thetahat)

	}

	## LQLS (loo = 0 -> leave-one-out estimator)
	## loo = 0 by default and you only need to specify if you want a loo
	lqls <- function(y,x,h,loo=0){

		n <- length(y)
		q <- ncol(x)
		q1 <- q+1
		qx <- 2*q + (factorial(q)/(factorial(2)*factorial(q-2))) 
		qx1 <- qx+1
		ones <- rep(1,n)
		ones1 <- rep(1,n-1)
		epsilon <- 0.00001

		thetahat <- matrix(0,ncol=qx1,nrow=0)

		for(j in 1:n){ 

 			if (loo==1){

				XX <- ones1
				KK <- ones1
 
				for (jj in 1:q){
       
					dx <- (x[-j,jj] - x[j,jj])/h[jj]
					KK <- KK*kk(dx)

					xt <- x[-j,jj] - x[j,jj]
					XX <- cbind(XX,xt)
				}

				## constructing additional components of XX (LQLS)
				## squares first, then interactions
				XX <- cbind(XX,XX[,2:q1]^2)

				for (jj in 2:q){

					jj1 <- jj+1
					XX <- cbind(XX,XX[,jj]*XX[,jj1:q1])

				}

				XK <- ones1*KK

				for (jj in 2:qx1){

					XK <-cbind(XK,XX[,jj]*KK)

				}

 
 			} else {

				XX <- ones
				KK <- ones

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					xt <- x[,jj] - x[j,jj]
					XX <- cbind(XX,xt)
	
					KK <- KK*kk(dx)
				}

				## constructing additional components of XX (LQLS)
				## squares first, then interactions
				XX <- cbind(XX,XX[,2:q1]^2)

				for (jj in 2:q){

					jj1 <- jj+1
					XX <- cbind(XX,XX[,jj]*XX[,jj1:q1])

				}

				XK <- ones*KK

				for (jj in 2:qx1){

					XK <-cbind(XK,XX[,jj]*KK)

				}

			}

			XK <- t(XK)

			## Ridging (if necessary)
			ridge <- 0

			while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,qx1)))),
					error = function(e){
					return(matrix(FALSE,qx1,qx1))
					})[1,1]==FALSE) {
				ridge <- ridge + epsilon
			}

			XKX <- XK%*%XX + diag(rep(ridge,qx1))
			XKXI <- solve(XKX)

			if (loo==1){

				XKY <- XK%*%y[-j]

			} else {

				XKY <- XK%*%y

			}

 	      	thetahat.estimate <- XKXI%*%XKY	
			thetahat.estimate <- t(thetahat.estimate)
			thetahat <- rbind(thetahat,thetahat.estimate)
		}

		return(thetahat)

	}

############################################################################

## Bandwidth Selectors (Each estimator for each type)

############################################################################

	## LSCV - LCLS
	lcls.lscv <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			thetahat.loo<-lcls(y,x,h,loo=1)
     			mhat<-thetahat.loo[,1]
 			
 			return(sum((y - mhat)^2)/n)
 
 		}
 
 	}

	## LSCV - LLLS
	llls.lscv <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			thetahat.loo<-llls(y,x,h,loo=1)
     			mhat<-thetahat.loo[,1]
 			
 			return(sum((y - mhat)^2)/n)
 
 		}
 
 	}

	## LSCV - LQLS
	lqls.lscv <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			thetahat.loo<-lqls(y,x,h,loo=1)
     			mhat<-thetahat.loo[,1]
 			
 			return(sum((y - mhat)^2)/n)
 
 		}
 
 	}

############################################################################

	## AICC for each

############################################################################

	## AICc - LCLS
	lcls.aicc <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			n <- length(y)
			q <- ncol(x)

			thetahat <- matrix(0,ncol=1,nrow=n)

			thetahat<-lcls(y,x,h)
     			mhat<-thetahat[,1]

			## creating the H matrix, faster if in lcls function
			HH <- matrix(0,ncol=n,nrow=n)

			for (j in 1:n){

				KK <- rep(1,n)

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					KK <- KK*kk(dx)

				}

				HH[,j] <- KK/sum(KK)

			}

			zih <- (1+(sum(diag(HH))/n))/(1-((sum(diag(HH))+2)/n))
 			
 			return((1/n)*sum(((y - mhat)^2)*exp(zih)))
 
 		}
 
 	}

	## AICc - LLLS
	llls.aicc <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			n <- length(y)
			q <- ncol(x)
			q1 <- q+1
			zqmat <- rep(0,q)
			ed <- cbind(1,t(zqmat))
			ones <- rep(1,n)
			epsilon <- 0.00001

			thetahat <- matrix(0,ncol=q1,nrow=n)

			thetahat<-llls(y,x,h)
     			mhat<-thetahat[,1]

			## creating the H matrix, faster if in llls function
			HH <- matrix(0,ncol=n,nrow=n)

			for (j in 1:n){

				XX <- ones 
				KK <- ones 

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					xt <- x[,jj] - x[j,jj]
					XX <- cbind(XX,xt)
	
					KK <- KK*kk(dx)
				}

				XK <- ones*KK

				for (jj in 1:q){

					XK <-cbind(XK,XX[,jj+1]*KK)

				}

				XK <- t(XK)

				## Ridging (if necessary)
				ridge <- 0

				while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,q1)))),
						error = function(e){
						return(matrix(FALSE,q1,q1))
						})[1,1]==FALSE) {
					ridge <- ridge + epsilon
				}

				XKX <- XK%*%XX + diag(rep(ridge,q1))
				XKXI <- solve(XKX)

				HH[,j] <- t(ed%*%(XKXI%*%XK))

			}

			zih <- (1+(sum(diag(HH))/n))/(1-((sum(diag(HH))+2)/n))
 			
 			return((1/n)*sum(((y - mhat)^2)*exp(zih)))
 
 		}
 
 	}

	## AICc - LQLS
	lqls.aicc <- function(h,y,x){
 
 		if(min(h) <= 0){
 
 			return(.Machine$double.xmax)
 
 		} else {

			n <- length(y)
			q <- ncol(x)
			q1 <- q+1
			qx <- 2*q + (factorial(q)/(factorial(2)*factorial(q-2))) 
			qx1 <- qx+1
			zqmat <- rep(0,qx)
			ed <- cbind(1,t(zqmat))
			ones <- rep(1,n)
			epsilon <- 0.00001

			thetahat <- matrix(0,ncol=qx1,nrow=0)

			thetahat<-lqls(y,x,h)
     			mhat<-thetahat[,1]

			## creating the H matrix, faster if in llls function
			HH <- matrix(0,ncol=n,nrow=n)

			for (j in 1:n){

				XX <- ones
				KK <- ones

				for (jj in 1:q){
       
					dx <- (x[,jj] - x[j,jj])/h[jj]
					xt <- x[,jj] - x[j,jj]
					XX <- cbind(XX,xt)
	
					KK <- KK*kk(dx)
				}

				## constructing additional components of XX (LQLS)
				## squares first, then interactions
				XX <- cbind(XX,XX[,2:q1]^2)

				for (jj in 2:q){

					jj1 <- jj+1
					XX <- cbind(XX,XX[,jj]*XX[,jj1:q1])

				}

				XK <- ones*KK

				for (jj in 2:qx1){

					XK <-cbind(XK,XX[,jj]*KK)

				}

				XK <- t(XK)

				## Ridging (if necessary)
				ridge <- 0

				while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,qx1)))),
						error = function(e){
						return(matrix(FALSE,qx1,qx1))
						})[1,1]==FALSE) {
					ridge <- ridge + epsilon
				}

				XKX <- XK%*%XX + diag(rep(ridge,qx1))
				XKXI <- solve(XKX)

				HH[,j] <- t(ed%*%(XKXI%*%XK))

			}

			zih <- (1+(sum(diag(HH))/n))/(1-((sum(diag(HH))+2)/n))
 			
 			return((1/n)*sum(((y - mhat)^2)*exp(zih)))
 
 		}
 
 	}

###############################################################################

## Counterfactual regressions (LCLS, LLLS, LQLS)

###############################################################################

	## LCLS 
	lcls.counter <- function(y,x,x.eval,h){

		n <- length(y)
		n.eval <- nrow(x.eval)
		q <- ncol(x)
		q1 <- q+1
		ones <- rep(1,n)
		ones.eval <- rep(1,n.eval)
		epsilon <- 0.00001
    
    h.prod <- 1
    
    for(jj in 1:q){
      h.prod <- h.prod*h[jj]
    }

		thetahat <- matrix(0,ncol=q1,nrow=n.eval)

		for(j in 1:n.eval){ 

			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x.eval[j,jj])/h[jj]
				KK <- KK*kk(dx)
			}

			thetahat[j,1] <- sum(y*KK)/sum(KK)

			dkdx <- matrix(0,ncol=q,nrow=n)

			for (jj in 1:q){

				dx2 <- (x[,jj] - x.eval[j,jj])/(h[jj]^2)
				dkdx[,jj] <- KK*dx2
			}

			for (jj in 1:q){

				thetahat[j,jj+1] <- ((sum(y*dkdx[,jj])*sum(KK)) - (sum(y*KK)*sum(dkdx[,jj])))/((sum(KK))^2)

			}

		}

		return(thetahat)

	}

	## LLLS 
	llls.counter <- function(y,x,x.eval,h){

		n <- length(y)
		n.eval <- nrow(x.eval)
		q <- ncol(x)
		q1 <- q+1
		ones <- rep(1,n)
		ones.eval <- rep(1,n.eval)
		epsilon <- 0.00001

		thetahat <- matrix(0,ncol=q1,nrow=0)

		for(j in 1:n.eval){ 

			XX <- ones
			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x.eval[j,jj])/h[jj]
				xt <- x[,jj] - x.eval[j,jj]
				XX <- cbind(XX,xt)
	
				KK <- KK*kk(dx)
			}

			XK <- ones*KK

			for (jj in 1:q){

				XK <-cbind(XK,XX[,jj+1]*KK)

			}

			XK <- t(XK)

			## Ridging (if necessary)
			ridge <- 0

			while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,q1)))),
					error = function(e){
					return(matrix(FALSE,q1,q1))
					})[1,1]==FALSE) {
				ridge <- ridge + epsilon
			}

			XKX <- XK%*%XX + diag(rep(ridge,q1))
			XKXI <- solve(XKX)

			XKY <- XK%*%y

 	      	thetahat.estimate <- XKXI%*%XKY	
			thetahat.estimate <- t(thetahat.estimate)
			thetahat <- rbind(thetahat,thetahat.estimate)
		}

		return(thetahat)

	}

	## LQLS 
	lqls.counter <- function(y,x,x.eval,h){

		n <- length(y)
		n.eval <- nrow(x.eval)
		q <- ncol(x)
		q1 <- q+1
		qx <- 2*q + (factorial(q)/(factorial(2)*factorial(q-2))) 
		qx1 <- qx+1
		ones <- rep(1,n)
		ones.eval <- rep(1,n.eval)
		epsilon <- 0.00001

		thetahat <- matrix(0,ncol=qx1,nrow=0)

		for(j in 1:n.eval){ 

			XX <- ones
			KK <- ones

			for (jj in 1:q){
       
				dx <- (x[,jj] - x.eval[j,jj])/h[jj]
				xt <- x[,jj] - x.eval[j,jj]
				XX <- cbind(XX,xt)
	
				KK <- KK*kk(dx)
			}

			## constructing additional components of XX (LQLS)
			## squares first, then interactions
			XX <- cbind(XX,XX[,2:q1]^2)

			for (jj in 2:q){

				jj1 <- jj+1
				XX <- cbind(XX,XX[,jj]*XX[,jj1:q1])

			}

			XK <- ones*KK

			for (jj in 2:qx1){

				XK <-cbind(XK,XX[,jj]*KK)

			}

			XK <- t(XK)

			## Ridging (if necessary)
			ridge <- 0

			while(tryCatch(as.matrix(solve(XK%*%XX+diag(rep(ridge,qx1)))),
					error = function(e){
					return(matrix(FALSE,qx1,qx1))
					})[1,1]==FALSE) {
				ridge <- ridge + epsilon
			}

			XKX <- XK%*%XX + diag(rep(ridge,qx1))
			XKXI <- solve(XKX)

			XKY <- XK%*%y

 	      	thetahat.estimate <- XKXI%*%XKY	
			thetahat.estimate <- t(thetahat.estimate)
			thetahat <- rbind(thetahat,thetahat.estimate)
		}

		return(thetahat)

	}




