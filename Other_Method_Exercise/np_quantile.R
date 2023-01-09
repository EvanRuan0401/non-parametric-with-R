library(np)
library(quantreg)

data(Engel95)
attach(Engel95)

x <- cbind(leisure,logwages)
x <- x[order(x[,2]),]
Leisure <- x[57:1654,1]
Logwage <- x[57:1654,2]

plot(Logwage,Leisure,ylab="Leisure",xlab="Log(Wage)",col="blue")
taus <- c(0.25,0.50,0.75)
for (i in 1:length(taus)){
abline(rq(Leisure~Logwage,tau=taus[i]),col="red",lwd=2)
}
legend("topright",c("tau=0.25, 0.50, 0.75 (bottom to top)"),bty="n")

rot <- c(1.059*sd(Leisure)*length(Leisure)^(-1/6),1.059*sd(Logwage)*length(Logwage)^(-1/6))
bw <- npcdistbw(formula=Leisure~Logwage,bandwidth.compute=FALSE,bws=rot)

model.q0.25 <- npqreg(bws=bw,tau=0.25,itmax=5)
model.q0.50 <- npqreg(bws=bw,tau=0.50,itmax=5)
model.q0.75 <- npqreg(bws=bw,tau=0.75,itmax=5)

plot(Logwage,Leisure,col="blue")
lines(Logwage,model.q0.25$quantile,col="red",lwd=2)
lines(Logwage,model.q0.50$quantile,col="red",lwd=2)
lines(Logwage,model.q0.75$quantile,col="red",lwd=2)
legend("topright",c("tau=0.25, 0.50, 0.75 (bottom to top)"),bty="n")

