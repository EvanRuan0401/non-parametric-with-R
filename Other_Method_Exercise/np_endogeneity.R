library(AER)
library(minqa)
library(np)

data(CigarettesSW)

CigarettesSW$rprice <- with(CigarettesSW,price/cpi)
CigarettesSW$rincome <- with(CigarettesSW,income/population/cpi)
CigarettesSW$tdiff <- with(CigarettesSW,(taxs-tax)/cpi)
c1995 <- subset(CigarettesSW,year=="1995")

attach(c1995)

model.iv <- ivreg(formula=log(packs)~log(rprice)+log(rincome)|log(rincome)+tdiff)

summary(model.iv)

###

model.lm <- lm(formula=log(packs)~log(rprice)+log(rincome))

summary(model.lm)

#####

Y <- log(packs)
X  <- log(rprice)
Z1 <- log(rincome)
Z2 <- tdiff
Z  <- data.frame(Z1,Z2)
W <- data.frame(X,Z1)
n <- length(Y)

p1 <- 3
p2 <- 1
d1 <- 1
d2 <- 1
dx <- 1
su2 <- SUh(p1=p1,p2=p2,d1=d1,d2=d2,dx=dx)

bs <- npregbw(formula=Y~X+Z1,regtype="ll")
bw2 <- bs$bw
bw1 <- bw2*n^(-(su2$upper-su2$lower)/(2*su2$gamma))

lpls.X <- lpls(y=X,x.train=Z,h=bw1,p=p1,kr="gaussian")
XX <- as.data.frame(X)
ZZ1 <- as.data.frame(Z1)
Uhat <- X-lpls.X[,1]
X.train <- data.frame(X=XX,Z1=ZZ1,U=Uhat)

ghat <- matrix(0,n,1)
ghat.prime <- matrix(0,n,2)

for (j in 1:n){

X.eval <- expand.grid(X=XX[j,],Z1=ZZ1[j,],U=sort(Uhat))
X.eval <- as.matrix(X.eval)
lpls.Y <- lpls(y=Y,x.train=X.train,x.eval=X.eval,h=c(bw2,1.059*sd(Uhat)*n^(-1/7)),p=p2,kr="gaussian")
ghat[j] <- mean(lpls.Y[,1])
ghat.prime[j,] <- colMeans(lpls.Y[,2:3])

}

beta.dist.IV <- ghat.prime[,1]
mean(beta.dist.IV)
median(beta.dist.IV)

###

lpls.Y <- lpls(y=Y,x.train=W,h=bw2,p=p2,kr="gaussian")

beta.dist.noIV <- lpls.Y[,2]
mean(beta.dist.noIV)
median(beta.dist.noIV)

