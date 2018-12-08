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

## In non-linear regression models, such as the probit model, coefficients cannot be interpreted as marginal effects. In probit model, the coefficients could 
## be interpreted with respect to the latent variable but since it has no specific unit, we cannot use any meaningful scale.
##Therefore, the coefficients give only the signs of the partial effect of each x on the response probability, and its statistical sygnificance by giving 
## us the information whether we can reject null hypothesis at sufficiently small sygnificance level. In our model the estimated coefficient of ldl is 0.17 
## and is sygnifficant at < 1% level. Conclusion: The ldl-cholesterol increases the chances of suffering from the disease and is higly sygnifficnat.  

## b) [0.5P] Are the S.E. valid, or do you need to adjust them for heteroscedasticity? Explain

##Under restrictive assumptions (i.e which includes homoscedasticity of error term), the Maximum Likelihood Estimator
## is consistent and asymptotically efficient. Therefore, we assume the assimptotic normality of coefficients;the SE are valid 
##and there is no need to test for heteroscedasticity.  


## c) [2P] Re-estimate the model from a) but this time include age in addition to ldl. You see that the estimated coefficient of ldl changes. 
##  Explain why? Additionally, show that your explanation is supported by the data. 

probit2<-glm(data$chd~data$ldl+data$age,family=binomial(link="probit"))
summary(probit2)
probit5 <-glm(data$chd~data$ldl*data$age,family=binomial(link="probit"))
summary(probit5)
## In both models, the coefficients are estimated under ceteris paribus assumtion, meaning that while estimating each of them we control for all other variables. 
## If both variable were uncorralated, the coefficient of ldl should not chage in the second model; controling for age should have no impact on ldl. 
## Therefore, seeing that the two coefficients for ldk are diffrent we may suspect ldl correlation with age. Sense check: the older people get, the more likely it is 
## they will have higher level of ldl-cholesterol. Now, let's use data to prove it: 

correlation_ldl_age <- cor(data$ldl, data$age)
print(correlation_ldl_age)

##Bingo! Hiwever this question is for two pints, maybe too easy answer....

## d) Finally, estimate the model from a) but include ldl squared next to ldl as a control variable.

probit3<-glm(data$chd~poly(data$ldl,2, raw= TRUE),family=binomial(link="probit"))
summary(probit3)

## i. [1P] Based on the estimated coefficients from a) and d) draw the two resulting marginal probability effects of ldl as a function of 
## ldl for ldl ∈ [1; 15] next to each other.

summary(data$ldl)
### Subsetting the data ( there are 2 values that are outside of [1;15] range ? 
## Let's assume no

## Running the probit for both models

xb1 <-probit$linear.predictors
xb2 <-probit3$linear.predictors

mpe_1 <-(dnorm(xb1)*(probit$coefficients[2]))
mpe_2 <-(dnorm(xb2)*(probit3$coefficients[2]*probit3$coefficients[3]*data$ldl))

plot1<-plot(data$ldl,mpe_1)
plot2<-plot(data$ldl,mpe_2)

### Kevin, could you just continue to make a joint graph? plot2 doesn't makes sense to me it seems that MPE never reaches 0.


## ii. [0.5P] Are any of the marginal probability effects linear in ldl? Explain why.

## As we are dealing with non-linear model the Marginal Probability Effect is estimated for every observation sperepatly, it does not have linear effect. In the 
##To obtain the "overall" average value for the whole sample we can calculate the mean of all the MPEs (recommendated) or take the MPE of an average for the sample ovservarion, however that average 
## observation may acctually not exist in reality.  


## iii. [1P] What is the advantage of the marginal probability effect based on the estimation in d) over the one based on a)? Explain shortly.

## Adding ldl^2 variable allows us to estimate wether the marginal effect of ldl increases or decreases as the ldl is higher, 
## as oppose to model a) where we assume independet of the ldl value effect. 
## In estimation d) the marginal effect of ldl^2  has a negative coefficient, which tells us that the popability function of P(Y=1|ldl) is concaved and 
## will reach its maximum when dP/dldl is equal to 0.
## Interpretation: Until certain level the increase of ldl will increase the chances of the getting the heart disease, however above this level the
## impact of increase of ldl on the probability of getting heart disease will be smaller. Sense check : Once the individual reaches certain level of ldl it is
## very likely they have a heart disease, therefore even further increase of ldl has a smaller impact on probability of getting heart disease. 
## This interpretation is not possible based on the estimation a), where we can only observed that increase in the ldl level increases the probability 
## of getting the heart disease. NOTE: The ldl^2 coefficient is only statistical at 10% level, which could be a limitation.

## Kevin, maybe we can inclustrate this point with the graph as on the slide 83? 

## iv. [1.5P] Calculate and properly interpret both marginal probability effects for the mean value of ldl in the sample (You do not need to 
## compute standard errors).

Probit_at_ldl_mean1<-probitmfx(data$chd~data$ldl,data=data,atmean=TRUE)
Probit_at_ldl_mean2 <- probitmfx(data$chd~poly(data$ldl,2, raw= TRUE),data=data,atmean=TRUE)
Probit_at_ldl_mean1
Probit_at_ldl_mean2
## In the model a), dF/dx of ldl at mean is equal to 0.0609. 
## Interpretation: At less then 1% sygniffiance level, the marginal probability in the sample suggest that for an individual with the exact mean characteristics
## increase of their ldl-cholesterol by one unit, increases their probability of getting the heart disease by 6.09 procentage points. 

## In model d), dF/dx of ldl at mean is equal to 0.1362 and of ldl^2 at mean is equal to - 0.006. 

ldl_ldlsq <- (-0.136/(2*(-0.006)))
ldl_ldlsq
##Interpretation: At 10% sygnifficance level, for an individual with the exact mean characteristics increase of their ldl-cholesterol by one unit, 
## increases their probability of getting the heart disease by 11.33 procentage points. 

## Kevin, I am not sure about the last point. 
## 1) Do we take the level of sygnificance of the square variable.? 
## 2) Do we change the interpretation of ldl based on the ldl^2 ? What I used is here: https://www.statalist.org/forums/forum/general-stata-discussion/general/1408413-significance-level-of-quadratic-term
## We could also calculate it manually (see slide 83) 

## If I use mpe_fuction2 at the mean, I get:  


## v. [1P] Are any of the effects computed in iv), ceteris paribus effects? Explain shortly. 

## No, none of the computed in iv) effects are ceteris paribus.In the model from d), once we fix the value of ldl at 
## the  mean the value of ldl^2 is also automatically fixed at (mean of ldl)^2. Therefore,in both models there is only one variable, which can cause change in the result
## and that is ldl.This means that once we fix the value of ldl at mean, the is no other variables that we need to control for
## - there is no chance for change in the outcome caused by another variable .

