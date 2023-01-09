
library("np")

dat <- read.csv("/Users/ellisruan/Desktop/hw/non_parametric/Homework/homework_4/data_for_homework_4.csv",header=TRUE)
attach(dat)

dat <- dat[order(dat[,2]),]
sell <- log(dat[,1])
lot <- log(dat[,2])
bdms <- dat[,3]
reg <- dat[,12]
# 變數已經變成排序好的
n <- length(sell)

#################### code for question 1 ####################

par.model <- lm(sell~lot+bdms+reg,x=TRUE,y=TRUE)

# 檢定參數與非參數設定 
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

data.train <- data.frame(lot,ordered(bdms),factor(reg))

########## a ##########

data.eval <- data.frame(lot=seq(from=min(lot),to=max(lot),length=100),bdms=ordered(rep(x=3,times=100)),reg=factor(rep(x=0,times=100)))
# 變數命名要與training set 一樣
# seq創立一個等距的數列 rep重複  
llls <- npreg(bws=llls.bw,residuals=TRUE,gradients=TRUE,data=data.train,newdata=data.eval)
# 要計算信賴區間要把residuals=TRUE marginal effects要gradient=TRUE
fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se
# 95%信心下之上下界 以上四個皆有100個結果與樣本數不符合 因此需要創造lot.star
lot.star <- seq(from=min(lot),to=max(lot),length=100)

plot(lot.star,fit,type="l",col="black",xlab="Log(Lot Size)",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(10,12),lwd=2)
lines(lot.star,fit.lower,lty=5,col="black",lwd=2)
lines(lot.star,fit.upper,lty=5,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

##########

grad <- llls$grad[,1]
grad.se <- llls$gerr[,1]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(100,1))
# 創一條0的水平線
plot(lot.star,grad,type="l",col="black",xlab="Log(Lot Size)",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(-1,2),lwd=2)
lines(lot.star,grad.lower,lty=5,col="black",lwd=2)
lines(lot.star,grad.upper,lty=5,col="black",lwd=2)
lines(lot.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

########## b ##########

data.eval <- data.frame(lot=rep(x=uocquantile(x=lot,prob=0.5),times=6),bdms=ordered(sort(unique(bdms))),reg=factor(rep(x=0,times=6)))
# uocquantile中的prob=0.5 可以得到中位數
llls <- npreg(bws=llls.bw,residuals=TRUE,gradients=TRUE,data=data.train,newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se

bdms.star <- sort(unique(bdms))

plot(bdms.star,fit,type="p",pch=16,col="black",xlab="Number of Bedrooms",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(10,12),lwd=2)
points(bdms.star,fit.lower,col="black",lwd=2)
points(bdms.star,fit.upper,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")

##########

grad <- llls$grad[,2]
grad.se <- llls$gerr[,2]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

grad <- c(0,grad[2:length(grad)])#
grad.lower <- c(0,grad.lower[2:length(grad.lower)])#
grad.upper <- c(0,grad.upper[2:length(grad.upper)])#

zero <- array(0,dim=c(6,1))

plot(bdms.star,grad,type="p",pch=16,col="black",xlab="Number of Bedrooms",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector)",ylim=c(-0.5,0.5),lwd=2)
points(bdms.star,grad.lower,col="black",lwd=2)
points(bdms.star,grad.upper,col="black",lwd=2)
lines(bdms.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients"),col=c("black"),pch=c(16),bty="n")
# discrete variable 的第二個marginal effect 會等於第一個
# 第一個marginal effect=0應比較合理 若要畫圖請參考上面#三行
########## c ##########

data.eval <- data.frame(lot=rep(x=uocquantile(x=lot,prob=0.5),times=2),bdms=ordered(rep(x=3,times=2)),reg=factor(sort(unique(reg))))

llls <- npreg(bws=llls.bw,residuals=TRUE,gradients=TRUE,data=data.train,newdata=data.eval)

fit <- llls$mean
fit.se <- llls$merr
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se

reg.star <- sort(unique(reg))

plot(reg.star,fit,type="p",pch=16,col="black",xlab="Housing Location",ylab="Log(Price)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(10,12),lwd=2)
points(reg.star,fit.lower,col="black",lwd=2)
points(reg.star,fit.upper,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),pch=c(16,1),bty="n")
# sort()函數是對向量進行從小到大的排序
# order()函數返回的值表示位置，依次對應的是向量的最小值、次小值、第三小值……最大值等（位置索引）
##########

grad <- llls$grad[,3]
grad.se <- llls$gerr[,3]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(2,1))

plot(reg.star,grad,type="p",pch=16,col="black",xlab="Housing Location",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector)",ylim=c(-0.5,0.5),lwd=2)
points(reg.star,grad.lower,col="black",lwd=2)
points(reg.star,grad.upper,col="black",lwd=2)
lines(reg.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients"),col=c("black"),pch=c(16),bty="n")

#################### answers to questions 1, 2, and 3 ####################

# 1 We reject the null hypothesis that the linear parametric specification is appropriate. #

# 2 We reject the null hypothesis that the regressors are irrelevant, regardless of using an individual test or a joint test. #

# 3a The fitted couve shows that the housing price increases as the lot size increases. However, the housing price decreases as the lot size increases after the lot size is around 9.2. #
#    The gradients shows the marginal effects are positive and stable but become insignificant after the lot size is around 9.2. #

# 3b The fitted points show that essentially the housing price increases as the number of bedrooms increases. The gradients shows that the marginal effects are unstable. #

# 3c The fitted points show that a better housing location has a slightly higher housing price and a positive marginal effect. #

