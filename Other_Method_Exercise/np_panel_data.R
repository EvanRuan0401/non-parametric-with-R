library(plm)
library(Ecdat)
library(np)

data(Airline)
attach(Airline)

model.p <- plm(formula=log(cost)~
           log(output)+
           log(pf)+
           lf,
           model="within",
           index=c("airline","year"),
           data=Airline)

summary(model.p)

lcost <- as.numeric(log(cost))
loutput <- as.numeric(log(output))
lpf <- as.numeric(log(pf))
lf <- as.numeric(lf)

bw <- npregbw(formula=lcost~
              loutput+
              lpf+
              lf+
              factor(airline)+
              ordered(year),
              okertype="wangvanryzin",
              bwmethod="cv.ls",
              regtype="ll")

summary(bw)

model.np <- npreg(bws=bw,gradients=TRUE)

summary(model.np)

plot(bw,gradients=TRUE,plot.errors.method="bootstrap")


