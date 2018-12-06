### Assesment 1 , Part 2 

##Preparing the enviroment 
remove(list=ls())
library(foreign)
library(mfx)
library(margins)

## Instrictions: Load the dataset “south_african_heart_disease_data.dta” and estimate the effect of ldl-
##(bad)cholesterol (ldl) in blood on the probability of suffering from heart disease (chd equals 1 if one suffers from it). 

## Loading the data
data <- read.dta(file= "south_african_heart_disease_data.dta") 

## Estimating effect ch ~ ldl

probit<-glm(data$chd~data$ldl,family=binomial(link="probit"))
summary(probit)

##a) [0.5P] Can you learn anything from the estimated coefficients? Explain shortly. 

## In non-linear regression models, such as the probit model, coefficients cannot be interpreted as marginal effects.They could be interpreted with respect to 
## the latent variable but since it has no specific unit, we cannot use any meaningful scale.
##Therefore, the coefficients give only the signs of the partial effect of each x on the response probability, and its statistical sygnificance by giving 
## us the information whether we can reject null hypothesis at sufficiently small sygnificance level. In our model the estimated coefficient of ldl is 0.17 
## and is sygnifficant at < 1% level. Conclusion: The ldl-cholesterol increases the chances of suffering from the disease and is higly sygnifficnat.  

## b) [0.5P] Are the S.E. valid, or do you need to adjust them for heteroscedasticity? Explain

## The S.E of the estimated coefficient is the square root of its asymptotic variance. Therefore, the S.E will be valid, if the sample size is large. 
## In this case n = 462. We assume Asymptitic Normality of the Maximul Likelihood Estimator; the ML is asymptotically efficient, since it is asymptotically
## normally distributed with mean 0 and variance equal to inverse of the Fisher-Information Matrix. Since we assume Normality of the estimator, there is no need to test for 
## heteroscedastisity. In practice, this assumption not always holds.

## c) [2P] Re-estimate the model from a) but this time include age in addition to ldl. You see that the estimated coefficient of ldl changes. 
##  Explain why? Additionally, show that your explanation is supported by the data. 

probit2<-glm(data$chd~data$ldl+data$age,family=binomial(link="probit"))
summary(probit2)

## In both models, the coefficients are estimated under ceteris paribus assumtion, meaning that while estimating each of them we control for all other variables. 
## If both variable were uncorralated, the coefficient of ldl should not chage in the second model, controling for age should have no impact on ldl. 
## Therefore, seeing that the two coefficients for ild are diffrent we may suspect ldl correlation with age. Sense check: the older people get, the more likely it is 
## they will have higher level of ldl-cholesterol. Now, let's use data to prove it: 

correlation_ldl_age <- cor(data$ldl, data$age)
print(correlation_ldl_age)

##Bingo!

## d) Finally, estimate the model from a) but include ldl squared next to ldl as a control variable.
ldlsq<-data$ldl^2
probit3<-glm(data$chd~data$ldl+ldlsq,family=binomial(link="probit"))
summary(probit3)

## i. [1P] Based on the estimated coefficients from a) and d) draw the two resulting marginal probability effects of ldl as a function of 
## ldl for ldl ∈ [1; 15] next to each other.
   


## ii. [0.5P] Are any of the marginal probability effects linear in ldl? Explain why.

## As we are dealing with non-lineat model the Marginal Probability Effect is estimated for every observation
## sperepatly, it does not have linear effect. To obtain the "overall" average value for the whole sample we can calculate  
## the mean of all the MPEs (recommendated) or take the MPE of an average for the sample ovservarion. 


## iii. [1P] What is the advantage of the marginal probability effect based on the estimation in d) over the one based on a)? Explain shortly.

## Often the relationship between y and x is nonlinear. There are a variety of solutions. One solution is to add polynomial terms and the first one to look at is usually x2. But you should first look at a scatterplot of x and y; you should also look at the residuals from the linear model without the quadratic term. 


## iv. [1.5P] Calculate and properly interpret both marginal probability effects for the mean value of ldl in the sample (You do not need to 
## compute standard errors).
mp_1 <- probitmfx(data$chd~data$ldl, data=data)
mp_2 <- probitmfx(data$chd~data$ldl+ldlsq, data=data)
mp_1
mp_2


## v. [1P] Are any of the effects computed in iv), ceteris paribus effects? Explain shortly. 

## Yes.
