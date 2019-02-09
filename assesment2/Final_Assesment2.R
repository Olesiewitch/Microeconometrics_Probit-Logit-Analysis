
## title: "Microeconometrics 2018/2019 - 2nd Assignment"
## date: "9th of Feb, 2018"
## Instructor: Tomas Repasky
##
## Group Members: 
##
##  - Johannes Wagner  
##    Humboldt-Universitat Berlin - ID 598797 - Msc Statistics  
##    wagnejoh@hu-berlin.de
##
## - Malgorzata Paulina Olesiewicz  
##   Humboldt-Universitat Berlin - ID 598939 - Msc Statistics  
##   malgorzata.paulina.olesiewicz@student.hu-berlin.de
##
## - Kevin Hoppe  
##   Humboldt-Universitat Berlin - ID 598247 -  Msc Statistics  
##   Kevin.Hoppe1@web.de


## [this .R-file is based on our .Rmd-file that we used to create the PDF-output.
##  R-Code is left as is and LaTeX-formatted explanations are commented out]

###################################################################

### Setting up the enviroment
library(foreign, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(margins, quietly = TRUE)
library(mfx, quietly = TRUE)
library(caret,quietly = TRUE)
# install.packages("InformationValue")
library(InformationValue)
library(xtable)
library(olsrr)
#install.packages("nnet")
library(nnet)
library(dplyr)
#install.packages("mlogit")
library(mlogit)

# clear environment
rm(list = ls())
# Define Directory
dir <- "~/R/git/Micro_Ass1"

setwd(dir)
### loading the data and looking up meaning of the variables
data <- read.dta("loanapp1.dta")
(variables=cbind(colnames(data),attr(data, "var.labels")))
# Model Selection ---------------------------------------------------------
# For simplicity reasons we can use a simple OLS model for the selection process
# remove variables which are directly related to rejection and have no explanatory power given our social theory: unver
# remove missings for the purpose of model selection
data_woNA <- data[complete.cases(data),] %>% dplyr::select(-c(approve, action, unver))
# names(data_woNA)
# regress the depend variable on all regressors
fit <- lm(reject~.,data=data_woNA)
## Perform forward stepwise selection to find most important variables
#ols_step_forward_p(fit, details = FALSE, penter = 0.001)
## Select the four significant variables which can be explained trough our theory:
# inson, pubrec, white, obrat



##Task 2 :Results interpretation

## Calculate logit model with MPEs at means
dataSelect <- dplyr::select(data, c(reject, sch, married, netw, inson, pubrec, white, obrat))
dataSelect <- dataSelect[complete.cases(dataSelect),] 
# logit_results <- logitmfx(reject~sch + married + inson + cons + pubrec + white + obrat, data = dataSelect, atmean=TRUE)
# logit_results
model=glm(reject~sch + married + netw + inson + pubrec + white + obrat, data =data, family = "binomial"(link="logit"))

## Calculate the probability for the average individual in our sample:
# Vector of Means
mean<-c(mean(dataSelect$sch),mean(dataSelect$married),mean(dataSelect$netw),mean(dataSelect$inson),mean(dataSelect$pubrec),mean(dataSelect$white),mean(dataSelect$obrat))
dim(mean)<-c(1,7)
names(mean)<-c("sch","married","netw","inson","pubrec","white","obrat");mean
# add Intercept
mean<-cbind(1,mean)
# Calculate P(Y=1|X) for Mean values in the sample
Prob1=1/(1+exp(-(model$coefficients%*%t(mean))))

##Task 3 :Sensitivity and Specificity

##1

predicted=predict(model,dataSelect, type="response")
sen_1=sensitivity(dataSelect$reject,predicted,threshold = 0.5)
spec_1=specificity(dataSelect$reject,predicted,threshold = 0.5)

##2
sen_2=sensitivity(dataSelect$reject,predicted,threshold = 0.3)
spec_2=specificity(dataSelect$reject,predicted,threshold = 0.3)

##3

cut= optimalCutoff(actuals = dataSelect$reject, predictedScores =predicted,optimiseFor="Both")
sen_3=sensitivity(dataSelect$reject,predicted,threshold =cut)
spec_3=specificity(dataSelect$reject,predicted,threshold =cut)

### Taks4:  Multinomial Logit

## Load data and extract the 7 variables that we chose 
## as a group:
##
## 1. Education (sch): [0,1]
## 2. marital status (married): [0,1]
## 3. Existence of private mortgage insurance (inson): [0,1]
## 4. Net worth (netw): [-7919, 28023]
## 5. Former bankruptcy (pubrec): [0,1]
## 6. Ethnic majority (white): [0,1]
## 7. Ratio of income and liabilities (obrat) : [0-95]
dat.red <- dplyr::select(data, reject, sch, married, inson, cons, pubrec, white, obrat, netw); as.tibble(dat.red)
## Coding categorial variables as dummy variables:
dummy.var <- function(var){
  
  if(min(unique(var), na.rm = T) == 0){
    var <- var + 1
  }
  
  len.var <- length(var)
  levels.var <- var %>% as.factor %>% levels %>% length
  levels.dummy <- levels.var - 1
  
  dummy <- matrix(NA,
                  nrow = len.var,
                  ncol = levels.dummy)
  for(j in 1:levels.var){
    
    match.loc <- which(var == unique(var)[j])
    
    for(i in match.loc){
      
      dummy.vec <- rep(0, levels.dummy)
      dummy.vec[unique(var)[j] - 1] <- 1
      dummy[i,] <- dummy.vec
    }
  }
  return(dummy)
}
sch <- dummy.var(dat.red$sch); head(sch)
married <- dummy.var(dat.red$married); head(married)
inson <- dummy.var(dat.red$inson); head(inson)
#cons <- dummy.var(dat.red$cons); head(cons)
pubrec <- dummy.var(dat.red$pubrec); head(pubrec)
white <- dummy.var(dat.red$white); head(white)
## Put together data set for regression models:
dat <- data.frame("reject" = dat.red$reject,
                  "educ" = sch,
                  "marry" = married,
                  "insur" = inson,
                  "netw" = dat.red$netw,
                  "bankr" = pubrec,
                  "white" = white,
                  "oblig" = dat.red$obrat)
dat.mnl <- mlogit.data(dat, shape = "wide", choice = "reject")
## Formulas for regression models:
form <- 'reject ~ educ + marry + insur + netw + bankr + white + oblig'
form.mnl <- mFormula(appr ~ educ + marry + insur + netw + bankr + white + oblig)
## Estimate a binomial logit model:
head(dat); binom <- glm(formula = form, data = dat, family = "binomial"(link = "logit"))
summary(binom)
## Estimate a multinomial logit model:
head(dat.mnl); multinom <- multinom(formula = form, data = dat)
summary(multinom)
## Comparison of coefficient estimates:
beta.diff <- cbind(coef(binom) - coef(multinom)); beta.diff





