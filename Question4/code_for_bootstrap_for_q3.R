
library("np")

dat <- read.csv("/Users/ellisruan/Desktop/hw/non_parametric/Homework/homework_4/data_for_homework_4.csv",header=TRUE)
attach(dat)

dat <- dat[order(dat[,2]),]
sell <- log(dat[,1])
lot <- log(dat[,2])
bdms <- dat[,3]
reg <- dat[,12]

n <- length(sell)

#################### code for question 1 ####################

par.model <- lm(sell~lot+bdms+reg,x=TRUE,y=TRUE)

hlr.test <- npcmstest(formula=sell~lot+ordered(bdms)+factor(reg),okertype="wangvanryzin",model=par.model,boot.method="wild",boot.num=199)
hlr.test

#################### code for question 2 ####################

llls.bw <- npregbw(formula=sell~lot+ordered(bdms)+factor(reg),okertype="wangvanryzin",regtype="ll",bwmethod="cv.ls")
llls.bw

#llls.bw$bw#
#c(2*sd(lot),1,0.5)#

rhl.test <- npsigtest(bws=llls.bw,index=c(1,2,3),boot.method="wild",boot.num=199)
#rhl.test <- npsigtest(bws=llls.bw,index=c(1,2,3),joint=TRUE,boot.method="wild",boot.num=199)#
rhl.test

#################### code for question 3 ####################
# plot(bws=llls.bw,plot.errors.method="bootstrap",plot.errors.boot.num=199,neval=100)
# lot 的時候 房間數大概 3 or 4 好區域 0

plot.out.r <- plot(bws=llls.bw,plot.errors.method="bootstrap",plot.errors.boot.num=199,neval=100,plot.behavior=c("data"))
# get regrssion data
plot.out.rg <- plot(bws=llls.bw,gradients=TRUE,plot.errors.method="bootstrap",plot.errors.boot.num=199,neval=100,plot.behavior=c("data"))
# get gradient data
# plot behavior 因為我們只要data而已 不畫出圖

########## a ##########

fit <- fitted(plot.out.r$r1)
fit.se <- se(plot.out.r$r1)
fit.lower <- fit+1.96*fit.se[,1]
fit.upper <- fit+1.96*fit.se[,2]
# se有兩個一正一負
lot.star <- seq(from=min(lot),to=max(lot),length=100)

plot(lot.star,fit,type="l",col="black",xlab="Log(Lot Size)",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(10,12),lwd=2)
lines(lot.star,fit.lower,lty=5,col="black",lwd=2)
lines(lot.star,fit.upper,lty=5,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

##########

grad <- gradients(plot.out.rg$rg1)
grad.se <- gradients(plot.out.rg$rg1,errors=TRUE)
grad.lower <- grad+1.96*grad.se[,1]
grad.upper <- grad+1.96*grad.se[,2]

zero <- array(0,dim=c(100,1))

plot(lot.star,grad,type="l",col="black",xlab="Log(Lot Size)",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(-1,2),lwd=2)
lines(lot.star,grad.lower,lty=5,col="black",lwd=2)
lines(lot.star,grad.upper,lty=5,col="black",lwd=2)
lines(lot.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

########## b ##########

fit <- fitted(plot.out.r$r2)
fit.se <- se(plot.out.r$r2)
fit.lower <- fit+1.96*fit.se[,1]
fit.upper <- fit+1.96*fit.se[,2]

bdms.star <- sort(unique(bdms))

plot(bdms.star,fit,type="p",pch=16,col="black",xlab="Number of Bedrooms",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(10,12),lwd=2)
points(bdms.star,fit.lower,col="black",lwd=2)
points(bdms.star,fit.upper,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")

##########

grad <- gradients(plot.out.rg$rg2)
grad.se <- gradients(plot.out.rg$rg2,errors=TRUE)
grad.lower <- grad+1.96*grad.se[,1]
grad.upper <- grad+1.96*grad.se[,2]

#grad <- c(0,grad[2:length(grad)])#
#grad.lower <- c(0,grad.lower[2:length(grad.lower)])#
#grad.upper <- c(0,grad.upper[2:length(grad.upper)])#

zero <- array(0,dim=c(6,1))

plot(bdms.star,grad,type="p",pch=16,col="black",xlab="Number of Bedrooms",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(-0.5,0.5),lwd=2)
points(bdms.star,grad.lower,col="black",lwd=2)
points(bdms.star,grad.upper,col="black",lwd=2)
lines(bdms.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")

########## c ##########

fit <- fitted(plot.out.r$r3)
fit.se <- se(plot.out.r$r3)
fit.lower <- fit+1.96*fit.se[,1]
fit.upper <- fit+1.96*fit.se[,2]

reg.star <- sort(unique(reg))

plot(reg.star,fit,type="p",pch=16,col="black",xlab="Housing Location",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(10,12),lwd=2)
points(reg.star,fit.lower,col="black",lwd=2)
points(reg.star,fit.upper,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")

##########

grad <- gradients(plot.out.rg$rg3)
grad.se <- gradients(plot.out.rg$rg3,errors=TRUE)
grad.lower <- grad+1.96*grad.se[,1]
grad.upper <- grad+1.96*grad.se[,2]

zero <- array(0,dim=c(2,1))

plot(reg.star,grad,type="p",pch=16,col="black",xlab="Housing Location",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Bootstrap Variance)",ylim=c(-0.5,0.5),lwd=2)
points(reg.star,grad.lower,col="black",lwd=2)
points(reg.star,grad.upper,col="black",lwd=2)
lines(reg.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")

#################### answers to questions 1, 2, and 3 ####################

# 1 We reject the null hypothesis that the linear parametric specification is appropriate. #

# 2 We reject the null hypothesis that the regressors are irrelevant, regardless of using an individual test or a joint test. #

# 3a The fitted couve shows that the housing price increases as the lot size increases. However, the housing price decreases as the lot size increases after the lot size is around 9.2. #
#    The gradients shows the marginal effects are positive and stable but become insignificant after the lot size is around 9.2. #

# 3b The fitted points show that essentially the housing price increases as the number of bedrooms increases. The gradients shows that the marginal effects are unstable. #

# 3c The fitted points show that a better housing location has a slightly higher housing price and a positive marginal effect. #

