library("np")
load(file = "/Users/ellisruan/Desktop/hw/non_parametric/npdata/mroz.rda")

mroz <- subset(mroz, lfp==1)
attach(mroz)

lnwage <- log(wage)

lnhours <- log(hours)
educ <- educ
largecity <- largecity
n <- length(lnwage)

# 
# cor(lnhours, lnwage)
# cor.test(lnhours, lnwage)
# plot(x = lnhours, y = lnwage)
# abline(lm(lnwage~lnhours),col="red")
# legend('topleft', legend = c('r = 0.00146889', 'p — value = 0.9758'))

# parametric model

par.model <- lm(lnwage ~ lnhours + educ + largecity, x=TRUE, y=TRUE)
summary(par.model)

# Correct Specification Test (only for LC)

hlr.test <- npcmstest(formula = lnwage ~ lnhours + ordered(educ) + factor(largecity),
                      okertype="wangvanryzin",
                      model=par.model,
                      boot.method="wild",
                      boot.num=199)
hlr.test


# Irrevelant Regressor Test

llls.bw <- npregbw(formula = lnwage ~ lnhours + ordered(educ) + factor(largecity),
                   okertype="wangvanryzin",
                   regtype="ll",
                   bwmethod="cv.ls")
llls.bw

rhl.test <- npsigtest(bws=llls.bw,
                      index=c(1,2,3),
                      boot.method="wild",
                      boot.num=199) 

rhl.test_joint <- npsigtest(bws=llls.bw,
                            index=c(1,2,3),
                            joint=TRUE,
                            boot.method="wild",
                            boot.num=199)
rhl.test
rhl.test_joint

# Interval of fitted curve and gradient (Asymptotic)

 ## how do lnhours effect lnwage (educ level at 12 years and locate in large city)
  ### fitted curve

data.train <- data.frame(lnhours,
                         ordered(educ),
                         factor(largecity))

data.eval <- data.frame(lnhours=seq(from=min(lnhours),to=max(lnhours),length=200),
                        educ=ordered(rep(x=12,times=200)),
                        largecity=factor(rep(x=1,times=200)))

llls <- npreg(bws=llls.bw,
              residuals=TRUE,
              gradients=TRUE,
              data=data.train,
              newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se
lnhours.star <- seq(from=min(lnhours),to=max(lnhours),length=200)

plot(lnhours.star,
     fit,type="l",
     col="black",
     xlab="Log(Hours)",
     ylab="Log(Wage)",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-5,5),
     lwd=2)
lines(lnhours.star,fit.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,fit.upper,lty=5,col="black",lwd=2)
legend("topleft",
       c("Local Linear Estimates","95% Confidence Intervals"),
       col=c("black","black"),
       lty=c(1,5),
       bty="n")

   ### gradient

grad <- llls$grad[,1]
grad.se <- llls$gerr[,1]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(200,1))

plot(lnhours.star,
     grad,type="l",
     col="black",
     xlab="Log(Hours)",
     ylab="Marginal Effects",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-5,5),
     lwd=2)
lines(lnhours.star,grad.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,grad.upper,lty=5,col="black",lwd=2)
lines(lnhours.star,zero,lty=1,col="black",lwd=2)
legend("topleft",
       c("Gradients","95% Confidence Intervals"),
       col=c("black","black"),
       lty=c(1,5),
       bty="n")


  ## How do lnhours effect lnwage (educ level at 16 years and locate in large city)
   ### fitted curve

data.train <- data.frame(lnhours,
                         ordered(educ),
                         factor(largecity))

data.eval <- data.frame(lnhours=seq(from=min(lnhours),to=max(lnhours),length=200),
                        educ=ordered(rep(x=16,times=200)),
                        largecity=factor(rep(x=1,times=200)))

llls <- npreg(bws=llls.bw,
              residuals=TRUE,
              gradients=TRUE,
              data=data.train,
              newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se
lnhours.star <- seq(from=min(lnhours),to=max(lnhours),length=200)

plot(lnhours.star,
     fit,type="l",
     col="black",
     xlab="Log(Hours)",
     ylab="Log(Wage)",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-5,5),
     lwd=2)
lines(lnhours.star,fit.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,fit.upper,lty=5,col="black",lwd=2)
legend("topleft",
       c("Local Linear Estimates","95% Confidence Intervals"),
       col=c("black","black"),
       lty=c(1,5),
       bty="n")

### gradient

grad <- llls$grad[,1]
grad.se <- llls$gerr[,1]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(200,1))

plot(lnhours.star,
     grad,type="l",
     col="black",
     xlab="Log(Hours)",
     ylab="Marginal Effects",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-5,5),
     lwd=2)
lines(lnhours.star,grad.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,grad.upper,lty=5,col="black",lwd=2)
lines(lnhours.star,zero,lty=1,col="black",lwd=2)
legend("topleft",
       c("Gradients","95% Confidence Intervals"),
       col=c("black","black"),
       lty=c(1,5),
       bty="n")


  ## how do educ effect lnwage (lnhours level at 0.75 and locate in large city)
   ### fitted curve

data.eval <- data.frame(lnhours=rep(x=uocquantile(x=lnhours,prob=0.75),times=13),
                        educ=ordered(sort(unique(educ))),
                        largecity=factor(rep(x=1,times=13)))

llls <- npreg(bws=llls.bw,
              residuals=TRUE,
              gradients=TRUE,
              data=data.train,
              newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se

educ.star <- sort(unique(educ))

plot(educ.star,
     fit,type="p",
     pch=16,col="black",
     xlab="Years of education",
     ylab="Log(Wage)",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-1,5),
     lwd=2)
points(educ.star,fit.lower,col="black",lwd=2)
points(educ.star,fit.upper,col="black",lwd=2)
legend("topleft",
       c("Local Linear Estimates","95% Confidence Intervals"),
       col=c("black","black"),
       pch=c(16,1),
       bty="n")

   ### gradient

grad <- llls$grad[,2]
grad.se <- llls$gerr[,2]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

grad <- c(0,grad[2:length(grad)])
grad.lower <- c(0,grad.lower[2:length(grad.lower)])
grad.upper <- c(0,grad.upper[2:length(grad.upper)])

zero <- array(0,dim=c(13,1))

plot(educ.star,
     grad,
     type="p",
     pch=16,
     col="black",
     xlab="Years of Education",
     ylab="Marginal Effects",
     main="LLLS (Cross-Validated Bandwidth Selector)",
     ylim=c(-0.5,0.5),
     lwd=2)
points(educ.star,grad.lower,col="black",lwd=2)
points(educ.star,grad.upper,col="black",lwd=2)
lines(educ.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients"),col=c("black"),pch=c(16),bty="n")




  ## how do city effect lnwage (lnhours level at 0.75 and educ level at 12)
   ### fitted curve
data.eval <- data.frame(lnhours=rep(x=uocquantile(x=lnhours,prob=0.75),times=2),
                        educ=ordered(rep(x=12,times=2)),
                        largecity=factor(sort(unique(largecity))))

llls <- npreg(bws=llls.bw,
              residuals=TRUE,
              gradients=TRUE,
              data=data.train,
              newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se

largecity.star <- sort(unique(largecity))

plot(largecity.star,
     fit,
     type="p",
     pch=16,
     col="black",
     xlab="CityLocation",
     ylab="Log(Wage)",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(-1,2),
     lwd=2)
points(largecity.star,fit.lower,col="black",lwd=2)
points(largecity.star,fit.upper,col="black",lwd=2)
legend("topleft",
       c("Local Linear Estimates","95% Confidence Intervals"),
       col=c("black","black"),
       pch=c(16,1),
       bty="n")

   ### gradient

grad <- llls$grad[,3]
grad.se <- llls$gerr[,3]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(2,1))

plot(largecity.star,
     grad,
     type="p",
     pch=16,
     col="black",
     xlab="Housing Location",
     ylab="Marginal Effects",
     main="LLLS (Cross-Validated Bandwidth Selector)",
     ylim=c(-0.5,0.5),lwd=2)
points(largecity.star,grad.lower,col="black",lwd=2)
points(largecity.star,grad.upper,col="black",lwd=2)
lines(largecity.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients"),col=c("black"),pch=c(16),bty="n")

  ## how do city effect lnwage (lnhours level at mean and educ level at 16)
   ### fitted curve

data.eval <- data.frame(lnhours=rep(x=uocquantile(x=lnhours,prob=0.75),times=2),
                        educ=ordered(rep(x=16,times=2)),
                        largecity=factor(sort(unique(largecity))))

llls <- npreg(bws=llls.bw,
              residuals=TRUE,
              gradients=TRUE,
              data=data.train,
              newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se

largecity.star <- sort(unique(largecity))

plot(largecity.star,
     fit,
     type="p",
     pch=16,
     col="black",
     xlab="CityLocation",
     ylab="Log(Wage)",
     main="LLLS (Cross-Validated Bandwidth Selector) \n CI (Asymptotic Variance)",
     ylim=c(0,3),
     lwd=2)
points(largecity.star,fit.lower,col="black",lwd=2)
points(largecity.star,fit.upper,col="black",lwd=2)
legend("topleft",
       c("Local Linear Estimates","95% Confidence Intervals"),
       col=c("black","black"),
       pch=c(16,1),
       bty="n")

   ### gradient
grad <- llls$grad[,3]
grad.se <- llls$gerr[,3]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(2,1))

plot(largecity.star,
     grad,
     type="p",
     pch=16,
     col="black",
     xlab="Housing Location",
     ylab="Marginal Effects",
     main="LLLS (Cross-Validated Bandwidth Selector)",
     ylim=c(-0.5,0.5),lwd=2)
points(largecity.star,grad.lower,col="black",lwd=2)
points(largecity.star,grad.upper,col="black",lwd=2)
lines(largecity.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients"),col=c("black"),pch=c(16),bty="n")

