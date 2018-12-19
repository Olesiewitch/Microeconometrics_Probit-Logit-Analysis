# title: "First Assignment"
# author: "Malgorzata, Kevin, Johannes"
# date: "20th of December, 2018"

#-------------------------------------------------------------------------------
# platform       x86_64-pc-linux-gnu         
# arch           x86_64                      
# os             linux-gnu                   
# system         x86_64, linux-gnu           
# language       R                           
# version.string R version 3.4.4 (2018-03-15)
#------------------------------------------------------------------------------


## Microeconometrics 2018/2019 
# Assesment 1

### Group Members: 
# Johannes Wagner, ID: 598797, Msc Statistics, <wagnejoh@hu-berlin.de>
#   Malgorzata Paulina Olesiewicz, ID:598939, Msc Statistics, <malgorzata.paulina.olesiewicz@student.hu-berlin.de>
#   Kevin Hope, ID: 598247, Msc Statistics, <Kevin.Hoppe1@web.de>                                       
#   ----------------------------------------------------------------------------

# Prepare the environment
library(tidyverse, quietly = TRUE) # for pipe operator %>%
library(mvtnorm, quietly = TRUE) # for joint multivariate Distribution
library(mfx, quietly = TRUE) # for estimating AMPE
library(foreign, quietly = TRUE)
library(margins, quietly = TRUE)
library(texreg, quietly = TRUE) # for printing the probit model

# Assigment 1 ------------------------------------------------------------------
# Define specific Directory
dir <- "~/R/git/Micro_Ass1/Markdown"

setwd(dir)
# clear environment
rm(list = ls())

# set seed for comparibility
set.seed(101)
# ------------------------------------------------------------------------------
## Define Parameters for the latent model
beta0 <- -30
beta1 <- 4
# Number of Observations 30.000
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
probit2<-glm(distr2[,"Y"] ~ distr2[,"xNull"], family=binomial(link="probit")) 
probit3<-glm(distr3[,"Y"] ~ distr3[,"xNull"], family=binomial(link="probit")) 

estBeta1 <- round(probit1$coefficients[2], 2)
estBeta2 <- round(probit2$coefficients[2], 2)
estBeta3 <- round(probit3$coefficients[2], 2)

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
  names(distr1) = c('xOne','error1')
  distr1$xOne <- round(distr1$xOne,2)
  yLatent1 <- beta0 + beta1 * distr1[,"xOne"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  probit1<-glm(distr1[,"Y"] ~ distr1[,"xOne"], family=binomial(link="probit"))
  betaEst1[i] <- probit1$coefficients[2]
  ## Population 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame()
  names(distr2) = c('xOne','error2')
  distr2$xOne <- round(distr2$xOne,2)
  yLatent2 <- beta0 + beta1 * distr2[,"xOne"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  probit2<-glm(distr2[,"Y"] ~ distr2[,"xOne"], family=binomial(link="probit"))
  betaEst2[i] <- probit2$coefficients[2]
  ## Population 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame()
  names(distr3) = c('xOne','error3')
  distr3$xOne <- round(distr3$xOne,2)
  yLatent3 <- beta0 + beta1 * distr3[,"xOne"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  probit3<-glm(distr3[,"Y"] ~ distr3[,"xOne"], family=binomial(link="probit"))
  betaEst3[i] <- probit3$coefficients[2]
}


## Make Plots for the answers
# Plot the kernel density estimates for beta based on the three populations ----
denBeta1 <- density(betaEst1, bw = "nrd0", adjust = 1.5, kernel = "gaussian")
denBeta2 <- density(betaEst2, bw = "nrd0", adjust = 1.5, kernel = "gaussian")
denBeta3 <- density(betaEst3, bw = "nrd0", adjust = 1.5, kernel = "gaussian")

png(filename="Plot2.png", height = 370)

plot(denBeta1$x,denBeta1$y,type="l",xlim = c(min(betaEst1) ,max(betaEst1)),
     xlab = "Estimated Coefficients of Beta 1", ylab = "Density", main ="Plot 2: Density Plot for Population 1")
## Generate a theoretical normal distribution with the parameters of Beta hat
rv <- rnorm(n, mean=mean(betaEst1),sd=sd(betaEst1)) %>% sort()
y <- dnorm(rv, mean=mean(betaEst1), sd=sd(betaEst1))
lines(rv, y, col = "indianred1")
abline(v=mean(betaEst1))
dev.off()

png(filename="Plot1.png")

plot(denBeta1$x,denBeta1$y,type="l",xlim = c(1.85 ,4.4), ylim = c(0, 13),
     xlab = "Estimated Coefficients", ylab = "Density", main ="Plot 1: Density Plot of Estimates of Beta 1")
lines(denBeta2$x,denBeta2$y, lty = 2)
lines(denBeta3$x,denBeta3$y, lty = 3)
legend("topright", legend=c("Population 1", "Population 2", "Population 3"), lty=1:3)
dev.off()

## Repeat for AMPES
## create empty vectors for our estimates
AMPEEst1 <- NULL
AMPEEst2 <- NULL
AMPEEst3 <- NULL


for (i in 1:n){
  ## Population 1
  distr1 <- rmvnorm(obs, mean = c(10,0) ,sigma = V1 , method = 'eigen') %>% as.data.frame()
  names(distr1) = c('xOne','error1')
  distr1$xOne <- round(distr1$xOne,2)
  yLatent1 <- beta0 + beta1 * distr1[,"xOne"] + distr1[,"error1"]
  y1 <- yLatent1 %>% replace(yLatent1<=0,0) %>% replace(yLatent1>0,1)
  distr1 <- distr1 %>% mutate(Y = y1)
  probitAMPE1 <- probitmfx(Y~xOne,data=distr1,atmean=FALSE)
  AMPEEst1[i] <- probitAMPE1$mfxest[1]
  ## Population 2
  distr2 <- rmvnorm(obs, mean = c(10,0) ,sigma = V2 , method = 'eigen') %>% as.data.frame()
  names(distr2) = c('xOne','error2')
  distr2$xOne <- round(distr2$xOne,2)
  yLatent2 <- beta0 + beta1 * distr2[,"xOne"] + distr2[,"error2"]
  y2 <- yLatent2 %>% replace(yLatent2<=0,0) %>% replace(yLatent2>0,1)
  distr2 <- distr2 %>% mutate(Y = y2)
  probitAMPE2 <- probitmfx(Y~xOne,data=distr2,atmean=FALSE)
  AMPEEst2[i] <- probitAMPE2$mfxest[1]
  ## Population 3
  distr3 <- rmvnorm(obs, mean = c(10,0) ,sigma = V3 , method = 'eigen') %>% as.data.frame()
  names(distr3) = c('xOne','error3')
  distr3$xOne <- round(distr3$xOne,2)
  yLatent3 <- beta0 + beta1 * distr3[,"xOne"] + distr3[,"error3"]
  y3 <- yLatent3 %>% replace(yLatent3<=0,0) %>% replace(yLatent3>0,1)
  distr3 <- distr3 %>% mutate(Y = y3)
  probitAMPE3 <- probitmfx(Y~xOne,data=distr3,atmean=FALSE)
  AMPEEst3[i] <- probitAMPE3$mfxest[1]
}

# c) Plot the kernel density estimates ------------------------------------

# adjust broder bandwidth to get a smooth function
denAMPE1 <- density(AMPEEst1, bw = "nrd0", adjust = 2, kernel = "gaussian")
denAMPE2 <- density(AMPEEst2, bw = "nrd0", adjust = 2, kernel = "gaussian")
denAMPE3 <- density(AMPEEst3, bw = "nrd0", adjust = 2, kernel = "gaussian")


## Make Plot
png(filename="Plot3.png")

plot(denAMPE1$x,denAMPE1$y,type="l",xlim = c(min(AMPEEst1),max(AMPEEst3)), ylim = c(0, 300),
     xlab = "Estimated Coefficients", ylab = "Density", main ="Plot 3: Density Plot of AMPES")
lines(denAMPE2$x,denAMPE2$y, lty = 2)
lines(denAMPE3$x,denAMPE3$y, lty = 3)
legend("topright", legend=c("Population 1", "Population 2", "Population 3"), lty=1:3)

clip(0,0.120,0,max(denAMPE1$y))
abline(v=mean(AMPEEst1))
clip(0,0.120,0,max(denAMPE2$y))
abline(v=mean(AMPEEst2))
clip(0,0.120,0,max(denAMPE3$y))
abline(v=mean(AMPEEst3))
dev.off()


## Calculate the relative differences
relDiff1 <- round(abs((mean(AMPEEst1) - mean(AMPEEst3)) / mean(AMPEEst3)*100),2)
relDiff2 <- round(abs((mean(AMPEEst1) - mean(AMPEEst2)) / mean(AMPEEst2)*100),2)


## Load the Dataset
data <- read.dta(file= "south_african_heart_disease_data.dta") 
## Estimating effect ch ~ ldl
probit<-glm(data$chd~data$ldl,family=binomial(link="probit"))
## Model 2
probit2<-glm(data$chd~data$ldl+data$age,family=binomial(link="probit"))

## Check for Correlation
correlation_ldl_age <- cor.test(data$ldl, data$age)
resultCor <-as.numeric(correlation_ldl_age[4]) %>% round(2)

## Model 3
probit3<-glm(data$chd~poly(data$ldl,2, raw= TRUE),family=binomial(link="probit"))

## Extract regression coefficients:

beta.1 <- as.numeric(probit$coefficients)
beta.2 <- as.numeric(probit3$coefficients)

## Prepare design matrices:

ldl <- seq(1, 15, 0.01)
ldl.square <- ldl^2
int <- rep(1, length(ldl))

designmat.1 <- as.matrix(data.frame("intercept" = int, 
                                    "ldl" = ldl))

designmat.2 <- as.matrix(data.frame("intercept" = int, 
                                    "ldl" = ldl, 
                                    "ldl.square" = ldl.square))

## Calculate linear predictor:

xb1 <- designmat.1 %*% beta.1
xb2 <- designmat.2 %*% beta.2

## Calculate marginal probability effect for each probit model:

mpe1 <- dnorm(xb1)*beta.1[2]
mpe2 <- dnorm(xb2)*(beta.2[2] + 2*beta.2[3]*ldl)

## Plot them together:
png(filename="Plot4.png", height = 350)

plot(x = ldl, 
     y = mpe1, 
     t= "l", 
     col="red", 
     lwd=2, 
     xlim=c(1,15),
     ylim=c(min(mpe1, mpe2), max(mpe1, mpe2)),
     main = "Plot 4: MPE for ldl in two probit models",
     ylab = "MPE")

lines(x = ldl, y = mpe2, t= "l", col="blue", lwd=2)
abline(v=mean(data$ldl))

legend("right", legend = c("ldl", "ldl + ldl²"), fill=c("red", "blue"))
dev.off()

## Calculate MPE at means
mean_ldl= mean(data$ldl)

# calculate MPEs
mpe1_mean <- dnorm(beta.1[1]+ beta.1[2]*mean_ldl)*(beta.1[2])
mpe2_mean <- dnorm(beta.2[1]+ beta.2[2]*mean_ldl+beta.2[3]*mean_ldl^2)*(beta.2[2] + 2*beta.2[3]*mean_ldl)

# get model parameters
Probit_at_ldl_mean1 <- probitmfx(data$chd~data$ldl,data=data,atmean=TRUE)
Probit_at_ldl_mean2 <- probitmfx(data$chd~data$ldl+I(data$ldl^2),data=data,atmean=TRUE)