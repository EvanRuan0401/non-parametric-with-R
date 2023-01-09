library(np)
library(gam)

data(wage1)
attach(wage1)

model.p <- lm(formula=lwage~
              educ+
              tenure+
              exper+
              expersq)

summary(model.p)

model.sp <- gam(formula=lwage~
                s(educ)+
                s(tenure)+
                s(exper))

par(mfrow=c(2,2))
plot(model.sp,se=TRUE)

# s(.) is estimated using spline smoothing #

