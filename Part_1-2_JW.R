# Assigment Part 1 and 2 JW

### Setting up the enviroment
library(foreign, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(margins, quietly = TRUE)
library(mfx, quietly = TRUE)
library(caret,quietly = TRUE)
# install.packages("InformationValue")
library(InformationValue)

library(MASS)
library(olsrr)

# Define Directory
dir <- "~/R/git/Micro_Ass1"

setwd(dir)
# clear environment
rm(list = ls())
### loading the data and looking up meaning of the variables

data <- read.dta("loanapp1.dta")
(variables=cbind(colnames(data),attr(data, "var.labels")))
names(data)

# Model Selection ---------------------------------------------------------
# For simplicity reasons we can use a simple OLS model for the selection process
# remove variables which are directly related to rejection and have no social explanatory power: unver
# remove missings for the purpose of model selection
data_woNA <- data[complete.cases(data),] %>% dplyr::select(-c(approve, action, unver))
names(data_woNA)
# regress the depend variable on all regressors
fit <- lm(reject~.,data=data_woNA)

## Backwards regression
selectedMod <- step(fit, direction = "backward")
summary(selectedMod)

## Perform forward stepwise selection to find most important variables
ols_step_forward_p(fit, details = FALSE, penter = 0.001)
# "reject"   "=1 if action == 3" 
## => model with ten variables; not all make sense though
## Order: how much they contribute to adjRsquare
# 2. inson: private mortgage insurance approved: positive approval history : Yes
# 4. cons: credit history on consumer stuff: check for credit history: Yes
# 3. pubrec: =1 if filed bankruptcy": check for negative record: Yes
# 5. typur: type of purchaser of loan: We dont really have information on that
# 6. white: control for racial discrimination: Yes
# 7. obrat: other oblgs,  % total inc: Further Information on reliabilities (ratio): Yes

## For interpretability reasons we include the social characteristics of age and marriage and we have seven

## We check for Multicolleanarity
fit2 <- lm(reject~ sch + married + inson+cons+ pubrec+  white+ obrat ,data=data_woNA)
all_vifs <- car::vif(fit2)
print(all_vifs)

## Calculate Probit model and MPEs at means
dataSelect <- dplyr::select(data, c(reject, sch, married, inson, cons, pubrec, white, obrat))
model=glm(reject~sch + married + inson + cons + pubrec + white + obrat, data=dataSelect, family = "binomial"(link="logit"))
summary(model)
(probitMPE <- logitmfx(reject~sch + married + inson + cons + pubrec + white + obrat,data=dataSelect,atmean=TRUE))

# Be aware that income is included in obrat
