
# Assigment 1 -------------------------------------------------------------
library(tidyverse)
library(mvtnorm)
library(mfx)


rm(list = ls())
## 1st Problem
## Simulate Probit Model with known Latent Variable

# set seed for comparibility
set.seed(101)
#*****************************************************************************
## Define Parameters
beta0 <- -30
beta1 <- 4
# Number of Observations 30.000 (use 400 for fast results)
obs = 30000
# Use package mvnorm for simulating multivariate Data with given vetcor of means
# and given Covariance
# method used described by Ripley (1987, p.98)
## Define Covariance Matrix V
V1 <- matrix(c(4,0,0,1), nrow=2, ncol=2)
V2 <- matrix(c(4,0,0,4), nrow=2, ncol=2)
V3 <- matrix(c(4,3,3,4), nrow=2, ncol=2)
distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
names(distr1) = c('xNull','error1')
distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame()
names(distr2) = c('xNull','error2')
distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame()
names(distr3) = c('xNull','error3')

# Latent Model for the joint multivariate Distribution
yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
## Expectation y = -30 + 4 * 10 = 10
yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
## transform Latent Variable into Bernoulli Variable
## y <- yLatent %>% 
y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
# Mean 0.92
y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
distr1 <- distr1 %>% mutate(Y = y1)
distr2 <- distr2 %>% mutate(Y = y2)
distr3 <- distr3 %>% mutate(Y = y3)


#*******************************************************************************
# Probit Model
probit1<-glm(distr1[,"Y"] ~ distr1[,"xNull"], family=binomial(link="probit")) 
probit1$coefficients[2]
summary(probit1)

probit2<-glm(distr2[,"Y"] ~ distr2[,"xNull"], family=binomial(link="probit")) 
probit2$coefficients[2]
summary(probit2)

probit3<-glm(distr3[,"Y"] ~ distr3[,"xNull"], family=binomial(link="probit")) 
probit3$coefficients[2]
summary(probit2)


# b.) Repeat estimation 400 times
n <- 400
# empty vector estimated ßs
betaEst1 <- NULL
betaEst2 <- NULL
betaEst3 <- NULL

for (i in 1:n){
  ## Distribution 1
  distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
  names(distr1) = c('xNull','error1')
  yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  probit1<-glm(distr1[,"Y"] ~ distr1[,"xNull"], family=binomial(link="probit")) 
  betaEst1[i] <- probit1$coefficients[2]
  ## Distribution 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame() 
  names(distr2) = c('xNull','error2')
  yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  probit2<-glm(distr2[,"Y"] ~ distr2[,"xNull"], family=binomial(link="probit")) 
  betaEst2[i] <- probit2$coefficients[2]
  ## Distribution 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame() 
  names(distr3) = c('xNull','error3')
  yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  probit3<-glm(distr3[,"Y"] ~ distr3[,"xNull"], family=binomial(link="probit")) 
  betaEst3[i] <- probit3$coefficients[2]
}

# Plot the Kernel Density
par(mfrow = c(1,3))
k1 <- density(betaEst1)
k2 <- density(betaEst2)
k3 <- density(betaEst3)

plot(k1$x,k1$y,type="l",xlim = c(0,12), xlab = "Estimated Beta j = 1", ylab = "Density", main ="Distribution 1")
plot(k2$x,k2$y,type="l",xlim = c(0,12), xlab = "Estimated Beta j = 2", ylab = "Density", main ="Distribution 2")
plot(k3$x,k3$y,type="l",xlim = c(0,12), xlab = "Estimated Beta j = 3", ylab = "Density", main ="Distribution 3")

# mean of betaEst2
mean(betaEst1)
mean(betaEst2)
mean(betaEst3)
summary(betaEst1)
summary(betaEst2)
summary(betaEst3)

# c.) Average marginal probability effect --------------------------------------
AMPEEst1 <- NULL
AMPEEst2 <- NULL
AMPEEst3 <- NULL


for (i in 1:n){
  ## Distribution 1
  distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
  names(distr1) = c('xNull','error1')
  yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  logit1 <- logitmfx(Y~xNull,data=distr1,atmean=FALSE)
  AMPEEst1[i] <- logit1$mfxest[1]
  ## Distribution 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame() 
  names(distr2) = c('xNull','error2')
  yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  logit2 <- logitmfx(Y~xNull,data=distr2,atmean=FALSE)
  AMPEEst2[i] <- logit2$mfxest[1]
  ## Distribution 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame() 
  names(distr3) = c('xNull','error3')
  yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  logit3 <- logitmfx(Y~xNull,data=distr3,atmean=FALSE)
  AMPEEst3[i] <- logit3$mfxest[1]
}

# Plot Densities for the vectors of AMPE
k1 <- density(AMPEEst1)
k2 <- density(AMPEEst2)
k3 <- density(AMPEEst3)


plot(k1$x,k1$y,type="l",xlim = c(0,0.2), xlab = "Estimated AMPE j = 1", ylab = "Density", main ="Distribution 1")
plot(k2$x,k2$y,type="l",xlim = c(0,0.2), xlab = "Estimated AMPE j = 2", ylab = "Density", main ="Distribution 2")
plot(k3$x,k3$y,type="l",xlim = c(0,0.2), xlab = "Estimated AMPE j = 3", ylab = "Density", main ="Distribution 3")

