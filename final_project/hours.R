
library("np")
load(file = "/Users/ellisruan/Desktop/npdata/mroz.rda")

mroz <- subset(mroz, lfp==1)
attach(mroz)

lnwage <- log(wage)

hours <- hours
educ <- educ
largecity <- largecity
n <- length(lnwage)


#################### code for OLS ####################


par.model <- lm(lnwage ~ hours + educ + largecity, x=TRUE, y=TRUE)
# x y =TRUE可以使變數相關資料進入par.model 再進入npcmstest
summary(par.model)

hlr.test <- npcmstest(formula = lnwage ~ hours + ordered(educ) + factor(largecity) ,
                      okertype="wangvanryzin",model=par.model,boot.method="wild",boot.num=199)
hlr.test
# 連續變數默認gaussian unorder默認AA order默認LR

#################### code for question 2 ####################

llls.bw <- npregbw(formula = lnwage ~ hours + ordered(educ) + factor(largecity)
                   ,okertype="wangvanryzin",regtype="ll",bwmethod="cv.ls")
llls.bw

#llls.bw$bw#show出解釋變數的bw
#c(2*sd(lot),1,0.5)#bw理論上限 連續 WR AA(根據c數目) 
# 都在上限內代表都以非線性方式影響y 

rhl.test <- npsigtest(bws=lcls.bw,index=c(1,2,3),boot.method="wild",boot.num=199) 
# index代表哪些要來做檢定
rhl.test_joint <- npsigtest(bws=lcls.bw,index=c(1,2,3),joint=TRUE,boot.method="wild",boot.num=199)#
rhl.test
rhl.test_joint

#################### code for question 3 ####################

data.train <- data.frame(lnhours,ordered(educ),factor(largecity))
# 僅表示educ是有序的資料

########## a ##########

data.eval <- data.frame(lnhours=seq(from=min(lnhours),to=max(lnhours),length=100),
                        educ=ordered(rep(x=12,times=100)),
                        largecity=factor(rep(x=1,times=100)))
# 變數命名要與training set 一樣
# seq創立一個等距的數列 rep重複  
llls <- npreg(bws=llls.bw,residuals=TRUE,gradients=TRUE,data=data.train,newdata=data.eval)
# 要計算信賴區間要把residuals=TRUE marginal effects要gradient=TRUE
fit <- llls$mean
fit.se <- llls$merr # 理論的se
fit.lower <- fit-1.96*fit.se
fit.upper <- fit+1.96*fit.se
# 95%信心下之上下界 以上四個皆有100個結果與樣本數不符合 因此需要創造lot.star
lnhours.star <- seq(from=min(lnhours),to=max(lnhours),length=100)

plot(lnhours.star,fit,type="l",col="black",xlab="Log(Hours)",ylab="Log(Wage)",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(-5,10),lwd=2)
lines(lnhours.star,fit.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,fit.upper,lty=5,col="black",lwd=2)
legend("topleft",c("Local Linear Estimates","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")

##########

grad <- llls$grad[,1]
grad.se <- llls$gerr[,1]
grad.lower <- grad-1.96*grad.se
grad.upper <- grad+1.96*grad.se

zero <- array(0,dim=c(100,1))
# 創一條0的水平線
plot(lnhours.star,grad,type="l",col="black",xlab="Log(Hours)",ylab="Marginal Effects",main="LLLS (Cross-Validated Bandwidth Selector) and CI (Asymptotic Variance)",ylim=c(-1,2),lwd=2)
lines(lnhours.star,grad.lower,lty=5,col="black",lwd=2)
lines(lnhours.star,grad.upper,lty=5,col="black",lwd=2)
lines(lnhours.star,zero,lty=1,col="black",lwd=2)
legend("topleft",c("Gradients","95% Confidence Intervals"),col=c("black","black"),lty=c(1,5),bty="n")
