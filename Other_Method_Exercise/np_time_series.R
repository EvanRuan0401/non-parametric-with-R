library(urca)
library(np)

data(Raotbl1)
attach(Raotbl1)

y <- r10y[3:144]
y1 <- r10y[2:143]
x <- r3m[3:144]
x1 <- r3m[2:143]
t <- seq(1,142,1)
t <- seq(1953.5,1988.75,0.25)

model.p <- lm(formula=y~y1+x+x1)

summary(model.p)

fit.p <- fitted(model.p)

bw <- npregbw(formula=y~y1+x+x1)

model.np <- npreg(bws=bw)

summary(model.np)

fit.np <- model.np$mean

plot(t,fit.p,type="l",lty=5,ylab="Estimated Long-Term Interest Rate",xlab="1953Q3-1988Q4",xlim=c(1953.5,1988.75),lwd=2)
lines(t,fit.np,col="red",lwd=2)
legend("topleft",c("Parametric Estimates","Nonparametric Estimates"),col=c("black","red"),lty=c(5,1),bty="n")

