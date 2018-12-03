
# Assigment 1 -------------------------------------------------------------
library(tidyverse)

rm(list = ls())
## 1st Problem
## Simulate Probit Model with known Latent Variable

# set seed for comparibility
set.seed(101)
## set observation to 400
#*****************************************************************************
## Define Parameters
beta0 <- -30
beta1 <- 4
xNorm1 <- rnorm(400,10,2) %>% round(2)

# Error term
errorTerm <- rnorm(400,0,1)

# Latent Model for the joint multivariate Distribution
yLatent <- beta0 + beta1 * xNorm1 + errorTerm
## transform Latent Variable into Bernoulli Variable
## y <- yLatent %>% 
y <- yLatent %>% replace(yLatent<=0,0) %>% replace(yLatent>0,1)

# yLatent[yLatent<=0] <- 0
# yLatent[yLatent>0] <- 1

# Probit Model
probit<-glm(yLatent ~ xNorm1, family=binomial(link="probit"))
probit$coefficients
