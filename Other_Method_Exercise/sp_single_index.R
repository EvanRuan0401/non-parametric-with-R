library(np)
library("MASS")

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

bw <- npindexbw(formula=lwage~
                female+
                married+
                educ+
                exper+
                expersq+
                tenure)

model.sp <- npindex(bws=bw)

summary(model.sp)

#####

data("birthwt")
attach(birthwt)

model.p <- glm(formula=low~
               smoke+
               race+
               ht+
               ui+
               ftv+
               age+
               lwt,
               family=binomial(link=probit))

summary(model.p)

cc <- table(low,ifelse(fitted(model.p)>0.5,1,0))
ccr <- sum(diag(cc))/sum(cc)
ccr

bw <- npindexbw(formula=low~
                factor(smoke)+
                factor(race)+
                factor(ht)+
                factor(ui)+
                ordered(ftv)+
                age+
                lwt,
                method="kleinspady")

model.sp <- npindex(bws=bw)

summary(model.sp)

