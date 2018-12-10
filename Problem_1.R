
# Assigment 1 ------------------------------------------------------------------
library(tidyverse) # for pipe operator %>%
library(mvtnorm) # for joint multivariate Distribution
library(mfx) # for estimating AMPE

# clear environment
rm(list = ls())

# set seed for comparibility
set.seed(101)
# ------------------------------------------------------------------------------
## Define Parameters for the latent model
beta0 <- -30
beta1 <- 4
# Number of Observations 30.000 (use 400 for fast results)
obs = 30000

# Simulate given Distributions -------------------------------------------------

# Use package mvnorm for simulating multivariate Data with given vetcor of means
# and given Covariance Matrix V
# method used described by Ripley (1987, p.98)

## Define Covariance Matrix V for Population j = 1,2,3
V1 <- matrix(c(4,0,0,1), nrow=2, ncol=2)
V2 <- matrix(c(4,0,0,4), nrow=2, ncol=2)
V3 <- matrix(c(4,3,3,4), nrow=2, ncol=2)

# 1.a) Simulate three Datasets with Y, xNull and error -------------------

## Simulate distributions with given means and Covariance; round xNull 
distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
names(distr1) = c('xNull','error1')
distr1$xNull <- round(distr1$xNull,2)
distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame()
names(distr2) = c('xNull','error2')
distr2$xNull <- round(distr2$xNull,2)
distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame()
names(distr3) = c('xNull','error3')
distr3$xNull <- round(distr3$xNull,2)

# Latent Model for the joint multivariate Distributions
yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
## transform Latent Variable into Bernoulli Variable
y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
distr1 <- distr1 %>% mutate(Y = y1)
distr2 <- distr2 %>% mutate(Y = y2)
distr3 <- distr3 %>% mutate(Y = y3)



# 1.a) Save the estimate for Betaj, j = 1,2,3 -----------------------------

# Probit Model
probit1<-glm(distr1[,"Y"] ~ distr1[,"xNull"], family=binomial(link="probit")) 
probit1$coefficients[2]
## beta1 = 3.93
summary(probit1)
confint(probit1)
## CI for beta: 3.74 lower bound and 4.14 upper bound for 95% confidence

probit2<-glm(distr2[,"Y"] ~ distr2[,"xNull"], family=binomial(link="probit")) 
probit2$coefficients[2]
## beta2 = 1.98
summary(probit2)
confint(probit2)
## CI for beta: 1.91 lower bound and 2.05 upper bound for 95% confidence

probit3<-glm(distr3[,"Y"] ~ distr3[,"xNull"], family=binomial(link="probit")) 
probit3$coefficients[2]
## beta3 = 3.61
summary(probit3)
confint(probit3)
## CI for beta: 3.46 lower bound and 3.77 upper bound for 95% confidence

# b.) Repeat estimation 400 times -----------------------------------------
## Number of repeated iterations n
n <- 400
# empty vectors for estimated betas for each population
betaEst1 <- NULL
betaEst2 <- NULL
betaEst3 <- NULL

for (i in 1:n){
  ## Population 1
  distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
  names(distr1) = c('xNull','error1')
  distr1$xNull <- round(distr1$xNull,2)
  yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  probit1<-glm(distr1[,"Y"] ~ distr1[,"xNull"], family=binomial(link="probit")) 
  betaEst1[i] <- probit1$coefficients[2]
  ## Population 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame() 
  names(distr2) = c('xNull','error2')
  distr2$xNull <- round(distr2$xNull,2)
  yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  probit2<-glm(distr2[,"Y"] ~ distr2[,"xNull"], family=binomial(link="probit")) 
  betaEst2[i] <- probit2$coefficients[2]
  ## Population 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame() 
  names(distr3) = c('xNull','error3')
  distr3$xNull <- round(distr3$xNull,2)
  yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  probit3<-glm(distr3[,"Y"] ~ distr3[,"xNull"], family=binomial(link="probit")) 
  betaEst3[i] <- probit3$coefficients[2]
}


# Plot the kernel density estimates for beta based on the three populations ----

par(mfrow = c(1,3))


# Explain the choices estimating the density function ---------------------
## Using a standard Gaussian Kernel and the standard deviation of the kernel 
## for bandwidth (I am not sure if we have to elaborate more; you can play around 
## with the Kernel and the bandwidth but I didn't get much better results than 
## with the default settings)
denBeta1 <- density(betaEst1, bw = "nrd0", adjust = 1, kernel = "gaussian")
denBeta2 <- density(betaEst2, bw = "nrd0", adjust = 1, kernel = "gaussian")
denBeta3 <- density(betaEst3, bw = "nrd0", adjust = 1, kernel = "gaussian")

plot(denBeta1$x,denBeta1$y,type="l",xlim = c(min(betaEst1) ,max(betaEst1)), 
     xlab = "Estimated Beta j = 1", ylab = "Density", main ="Population 1")
abline(v=mean(betaEst1))
plot(denBeta2$x,denBeta2$y,type="l",xlim = c(min(betaEst2) ,max(betaEst2)), 
     xlab = "Estimated Beta j = 2", ylab = "Density", main ="Population 2")
abline(v=mean(betaEst2))
plot(denBeta3$x,denBeta3$y,type="l",xlim = c(min(betaEst3) ,max(betaEst3)), 
     xlab = "Estimated Beta j = 3", ylab = "Density", main ="Population 3")
abline(v=mean(betaEst3))
## Question: Should we plot them next to each other or within one Plot?
## In the final Plots ranges should be calculated manually

# i) Does distribution of beta1 confirm to your expectations? ------------------
# The expectation is that estimates of beta are asymptotically normal 
# (Woolbridge 2010, p.43) and/or given the Central Limit Theorem that the estimates
# converge approximately to the true value of Beta (which is defined in the latent
# linear model: Beta = 4)
# Estimates are expected to be consistent because of zero mean condition of errors
# and Indepence of errors (Sources ?) 
# And the upper and lower limits of the distribution correspond to the CI of beta
# which was estimated in the first probit model
# all assumptions check out more or less! (actually given that 400 samples are not 
# that many - Do we want to try with 4000? (I assume that would need more then two 
# hours to compute))


# ii) Compare. Why do the means of the distributions differ? -------------------

mean(betaEst1)
## mean of 4.0
# this mean is explained trough the true value of the latent model
mean(betaEst2)
## mean of 2.0
# thgis mean is explained trough the fact of neglected heterogeneity with 
# betaEst = beta / variance (Woolbridge 2010, p. 614) 
# Also see slide 12 of the first exercise sessions
# [I didn't completly understand the formular yet, but it seems to be fitting
# our sample] [this belongs into iii]
mean(betaEst3)
## mean of 3.6
# We have Heteroskedasrticity, what causes estimates to be smaller then the true
# value. This needs more elaboration

# iii.) Explain the value of the mean of beta2 ----------------------------
mean(betaEst2)
## mean of 2.0
# this mean is explained trough the fact of neglected heterogeneity with 
# betaEst = beta / variance (Woolbridge 2010, p. 614) 
# Also see slide 12 of the first exercise sessions
# [I didn't completly understand the formular yet, but it seems to be fitting
# our sample] 


# c.) Average marginal probability effect --------------------------------------

## create empty vectors for our estimates
AMPEEst1 <- NULL
AMPEEst2 <- NULL
AMPEEst3 <- NULL


for (i in 1:n){
  ## Population 1
  distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame() 
  names(distr1) = c('xNull','error1')
  distr1$xNull <- round(distr1$xNull,2)
  yLatent1 <- beta0 + beta1 * distr1[,"xNull"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  logit1 <- logitmfx(Y~xNull,data=distr1,atmean=FALSE)
  AMPEEst1[i] <- logit1$mfxest[1]
  ## Population 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame() 
  names(distr2) = c('xNull','error2')
  distr2$xNull <- round(distr2$xNull,2)
  yLatent2 <- beta0 + beta1 * distr2[,"xNull"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  logit2 <- logitmfx(Y~xNull,data=distr2,atmean=FALSE)
  AMPEEst2[i] <- logit2$mfxest[1]
  ## Population 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame() 
  names(distr3) = c('xNull','error3')
  distr3$xNull <- round(distr3$xNull,2)
  yLatent3 <- beta0 + beta1 * distr3[,"xNull"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  logit3 <- logitmfx(Y~xNull,data=distr3,atmean=FALSE)
  AMPEEst3[i] <- logit3$mfxest[1]
}


# c) Plot the kernel density estimates ------------------------------------

# adjust broder bandwidth to get a smooth function
denAMPE1 <- density(AMPEEst1, bw = "nrd0", adjust = 2, kernel = "gaussian")
denAMPE2 <- density(AMPEEst2, bw = "nrd0", adjust = 2, kernel = "gaussian")
denAMPE3 <- density(AMPEEst3, bw = "nrd0", adjust = 2, kernel = "gaussian")


plot(denAMPE1$x,denAMPE1$y,type="l",xlim = c(min(AMPEEst1),max(AMPEEst1)), 
     xlab = "Estimated AMPE j = 1", ylab = "Density", main ="Population 1")
abline(v=mean(AMPEEst1))
plot(denAMPE2$x,denAMPE2$y,type="l",xlim = c(min(AMPEEst2),max(AMPEEst2)), 
     xlab = "Estimated AMPE j = 2", ylab = "Density", main ="Population 2")
abline(v=mean(AMPEEst2))
plot(denAMPE3$x,denAMPE3$y,type="l",xlim = c(min(AMPEEst3),max(AMPEEst3)), 
     xlab = "Estimated AMPE j = 3", ylab = "Density", main ="Population 3")
abline(v=mean(AMPEEst3))
## when asked for specific values  is the mean optimal?


# i) calculate the relative difference between AMPE1 and AMPE2 ------------

(relDiff1 <- abs((mean(AMPEEst1) - mean(AMPEEst3)) / mean(AMPEEst1)))
## relative difference of 25.11 percent
(mean(AMPEEst1))
# average effect of X1 on Y is 0.09 (from population 1), which we would use
# because only here we find no violations of the standard assumptions of the
# of the probit model


# ii) Calculate the relative difference between AMPE1 and AMPE2 -----------

(relDiff2 <- abs((mean(AMPEEst1) - mean(AMPEEst2)) / mean(AMPEEst1)))
## relative difference of 1.2 percent


# iii.) Please provide a detailed explanation for the results in i --------

## Puh, that gives three points, so I guess we need some details here; I have 
## some ideas but would like to discuss them with you tomorrow
## for ii.):: estimates are more or less consistent because the marginal 
## effect is a ration and the effect for betas (in a.) cancels out for 
## beta / variance (woolbridge 2010, p. 583)
## for i.) it's probably the result of heterogenity which causes an overestimation
## of marginal effect size -> why is that in particular?


##------------------------------------------------------------------------------
# Additional question: Can we explain why the density of the estimates 
# (beta and AMPE) for population 2 is higher than for 1 and 3 at their means? 
