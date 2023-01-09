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

bw <- npscoefbw(formula=lwage~
                educ+
                tenure+
                exper+
                expersq|factor(female)+factor(married))

model.sp <- npscoef(bws=bw,betas=TRUE)

summary(model.sp)

colMeans(coef(model.sp))

