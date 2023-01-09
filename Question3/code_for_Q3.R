
library("np")
library("minqa")
data("cps71")
attach(cps71)

y <- logwage
x <- age
x2 <- age^2

X <- cbind(x,x2)

n <- length(y)

#################### code for question 1 ####################

set.seed(1)

bw.x <- 1.059*sd(x)*n^(-1/5)
bw.x

lw <- lw.ls(y=y,x=x,x.p=X,bw=bw.x,nb=199)

#lw <- lw.ls(y=y,x=x,x.p=X,nb=199)# # LSCV
#bw.x <- lw$bw#
#bw.x#

boot.tstat <- lw$boot.tstat.norm

tstat <- lw$tstat.norm # stat from real data 
tstat

critical.value <- quantile(c(tstat,boot.tstat),0.95,type=4)
critical.value # reject h0

##rank.dist <- rank(c(tstat,boot.tstat))## 
# rank() return ordered index
##p.value <- 1-(rank.dist[1]/(199+1))##
p.value <- lw$p.value.norm
p.value

#################### code for question 3 ####################

set.seed(1)

v <- rnorm(n=n,mean=0,sd=1) # creat irrelevant

X <- cbind(x,v)

bw.X <- c(1.059*sd(x)*n^(-1/6),1.059*sd(v)*n^(-1/6))

lv <- lv.test(y=y,x=X,w=X[,1],bw.x=bw.X,bw.w=bw.X[1],nb=199) 
# x=全部regressor w=relevant regressor

#lv <- lv.test(y=y,x=X,w=X[,1],nb=199)# # LSCV

boot.tstat <- lv$boot.tstat.norm

tstat <- lv$tstat.norm
tstat

critical.value <- quantile(c(tstat,boot.tstat),0.95,type=4)
critical.value

##rank.dist <- rank(c(tstat,boot.tstat))##
##p.value <- 1-(rank.dist[1]/(199+1))##
p.value <- lv$p.value.norm
p.value
# fail to reject h0 v is irrelevant


###Q4比較h與兩倍標準差###
#lv$bw.x#
#c(2*sd(x),2*sd(v))#

#################### answers to questions 2 and 4 ####################

# 2 We reject the null hypothesis that the quadratic parametric specification is appropriate, and this implies the significance of the dip #
#   and is consistent with the changes in the confidence intervals. #

# 4 The cross-validated bandwidth of v is greater than two times the standard deviation of v, and this implies the irrevlance of v #
#   and is consistent with the fact that we fail to reject the null hypothesis that v is an irrelevant regressor. #


