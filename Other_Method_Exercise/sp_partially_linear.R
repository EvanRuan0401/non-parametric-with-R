library(np)

data(wage1)
attach(wage1)

model.p <- lm(formula=lwage~
              female+
              married+
              educ+
              tenure+
              exper+
              expersq)

summary(model.p)

bw <- npplregbw(formula=lwage~
                factor(female)+
                factor(married)+
                educ+
                tenure|exper)

model.sp <- npplreg(bws=bw)

summary(model.sp)

fit <- model.sp$mean
coef <- as.vector(coef(model.sp))
x <- cbind(female,married,educ,tenure)
fit.p <- x%*%coef
fit.np <- fit-fit.p

w <- cbind(fit.np,exper)
w <- w[order(w[,2]),]
fit.np <- w[,1]
exper <- w[,2]
plot(exper,fit.np,type="l",ylab="Marginal Effects",xlab="Experience")

