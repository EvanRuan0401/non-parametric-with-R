
dat70 <- read.csv("/Users/ellisruan/Desktop/hw/non_parametric/Homework/homework_1/data_for_homework_1(1).csv",header=TRUE)
attach(dat70)

rgdpch70 <- as.matrix(rgdpch70/mean(rgdpch70))
rgdpch70 <- rgdpch70[order(rgdpch70)]

#################### codes for questions 1 and 2 ####################

n70 <- length(rgdpch70)
bw70 <- 1.059*sd(rgdpch70)*n70^(-1/5)
#bw70 <- 2.345*sd(rgdpch70)*n70^(-1/5)#

KK.store.70 <- matrix(0,nrow=n70,ncol=n70)

for(j in 1:n70){

dx <- (rgdpch70-rgdpch70[j])/bw70
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)
#KK <- ifelse(abs(dx)>1,0,0.75*(1-dx^2))#
KK.store.70[,j] <- KK

}

est.den.70 <- solve(n70*bw70)*colSums(KK.store.70)

##########

dat00 <- read.csv("/Users/ellisruan/Desktop/hw/non_parametric/Homework/homework_1/data_for_homework_1(2).csv",header=TRUE)
attach(dat00)

rgdpch00 <- as.matrix(rgdpch00/mean(rgdpch00))
rgdpch00 <- rgdpch00[order(rgdpch00)]

n00 <- length(rgdpch00)
bw00 <- 1.059*sd(rgdpch00)*n00^(-1/5)
#bw00 <- 2.345*sd(rgdpch00)*n00^(-1/5)#

KK.store.00 <- matrix(0,nrow=n00,ncol=n00)

for(j in 1:n00){

dx <- (rgdpch00-rgdpch00[j])/bw00
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)
#KK <- ifelse(abs(dx)>1,0,0.75*(1-dx^2))#
KK.store.00[,j] <- KK

}

est.den.00 <- solve(n00*bw00)*colSums(KK.store.00)

##########

plot(rgdpch70,est.den.70,type="l",col="black",xlab="RGDP/POP",ylab="Density",main=" Kernel (Rule-of-Thumb Bandwidth)",xlim=c(0,5),ylim=c(0,0.6),lwd=2)
lines(rgdpch00,est.den.00,lty=5,col="black",lwd=2)
legend("topright",c("1970","2000"),col=c("black","black"),lty=c(1,5),bty="n")

#################### codes for questions 3 and 4 ####################

n <- n00

KK.store.loo <- matrix(0,nrow=n,ncol=n)
KK.store.convolution  <- matrix(0,nrow=n,ncol=n)

h <- seq(0.001,0.5,0.001)

lscv.loo <- matrix(0,nrow=length(h),ncol=1)

for(jj in 1:length(h)){

for(j in 1:n){

dx <- (rgdpch00-rgdpch00[j])/h[jj]
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)
KK.convolution <- solve(sqrt(4*pi))*exp(-0.25*dx^2)

KK.store.loo[,j] <- KK
KK.store.loo[j,j] <- 0
KK.store.convolution[,j] <- KK.convolution

}

lscv.loo[jj] <- solve(n^2*h[jj])*sum(KK.store.convolution)-2*solve(n*(n-1)*h[jj])*sum(KK.store.loo)

}

min.lscv.loo <- which(lscv.loo==min(lscv.loo))

##########

plot(h,lscv.loo,type="l",col="black",xlab="Bandwidth",ylab="LSCV",main="LSCV Function",xlim=c(0,0.5),ylim=c(-1,1),lwd=2)
plot(h,lscv.loo,type="l",col="black",xlab="Bandwidth",ylab="LSCV",main="LSCV Function",xaxt = "n", yaxt ="n",lwd=2) #沒苛度
#plot(h,lscv.loo,type="l",col="black",xlab="Bandwidth",ylab="LSCV",main="LSCV Function",xlim=c(0,0.5),ylim=c(-5,0),lwd=2)#

##########

KK.store <- matrix(0,nrow=n,ncol=n)

for(j in 1:n){

dx <- (rgdpch00-rgdpch00[j])/h[min.lscv.loo]
#dx <- (rgdpch00-rgdpch00[j])/#
KK <- solve(sqrt(2*pi))*exp(-0.5*dx^2)
KK.store[,j] <- KK

}

est.den <- solve(n*h[min.lscv.loo])*colSums(KK.store)
#est.den <- solve(n*)*colSums(KK.store)#

##########

plot(rgdpch00,est.den,type="l",lty=5,col="black",xlab="RGDP/POP",ylab="Density",main="Gaussian Kernel (Cross-Validated Bandwidth)",xlim=c(0,5),ylim=c(0,2.5),lwd=2)
legend("topright",c("2000"),col=c("black"),lty=c(5),bty="n")

#################### answers to questions 1, 2, and 3 ####################

# 1 the bimodal feature is apparent in either year. #

# 2 the bimodal feature is apparent in either year. #

# 3 the bimodal feature is not apparent. #

