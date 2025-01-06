# In this script we test if an intervention had an effect on a group of nursing applicants, 
# while controlling for background factors using multiple linear regression.
# The outcome variable (dependent variable) for the intervention is the PNPI2 score 
# as calculated and analysed in https://github.com/eglerean/sulo
#
# The various independent variables are the background factors as listed in the manuscript
# and the binary variable had access to the intervention tool yes/no
#
# Preprocessing for the data is done with custom python scripts available at 
# https://github.com/eglerean/sulo
#
# Original data cannot be shared due to legal and ethical restrictions. 
# If you want to reuse this code, prepare the data in tidy format: each row is a subject and the columns
# can be seen in section 5 below. Some columns are categorical and those are mapped in section 4.


## load needed libraries
# uncomment the missing libraries if you need to install them
#install.packages("olsrr")
library("olsrr")
#install.packages("fBasics")
library("fBasics")
#install.packages("lmtest")
library("lmtest") #dwtest
#install.packages("robustbase")
library("robustbase")
#install.packages("sandwhich")
library("sandwich")
#install.packages("vtable")
library("vtable")
#install.packages("KableExtra")
library("kableExtra")
#install.packages("dplyr")
library("dplyr")
#install.packages("ggpubr")
library("ggpubr")

## custom functions 
# Function that does quality control of the linear model

testbgf <- function(testcols,alldata){
  for (i in 1:length(testcols)){
    if(class(alldata[[testcols[i]]]) == "factor"){
      table2 <- table(alldata[[testcols[i]]])
      pptot <- prop.table(table2)
      sumouttot <- paste(pptot[2])
      table3 <- table(alldata$digital_platform, alldata[[testcols[i]]])
      pp <- prop.table(table3, margin = 1)
      sumout <- paste(pp[,2])
      res <- chisq.test(alldata[[testcols[i]]], alldata$digital_platform, correct=FALSE)
    }else
    {
      tottemp <- summary(alldata[[testcols[i]]])  
      sumouttot <- paste(as.numeric(tottemp), collapse=' ')
      ctemp <- tapply(alldata[[testcols[i]]], alldata$digital_platform, summary)
      ctrl_gr <- paste(as.numeric(ctemp[[1]]), collapse=' ')
      inter_gr <- paste(as.numeric(ctemp[[2]]), collapse=' ')
      
      sumout <- paste(ctrl_gr, inter_gr, collapse=' ', sep = ",")
      fmla <- formula(paste(testcols[i], " ~ digital_platform"))
      res <- wilcox.test(fmla,
                         data = alldata,
                         exact = FALSE)
    }
    #summary(res) # uncomment to see the output of the residuals
    temp <- paste(sumout, collapse = ' ')
    out <- paste(testcols[i], sumouttot, temp, res$p.value)
    print(out)
  }
  
}

lmqc <- function(model) {
  ## quality control for the fitted model
  if(class(model) == "lm"){
    BP1=ols_test_breusch_pagan(model) # non significant is good
    print(BP1)
  }
  # test residuals are normally distributed
  JT=jarqueberaTest(model$resid) 
  print(JT)
  # test that residuals are independent (p>0.05 we are fine)
  DW=dwtest(model) #Test for independence of residuals
  print(DW)
  BP2=bptest(model)
  print(BP2)
  # QC plots
  par(mfrow=c(2,2)) 
  plot(model)
  # check for skewness
  # rstudent returns the standardized residuals, we hope to see a normal distrib
  if(class(model) == "lm"){
    par(mfrow=c(2,1))
    hist(rstudent(model))
  }
}

## Script begins here

## 1. loading data
# load the data
data_file = "bg_and_scores.csv"
fulldata = read.table(data_file,header=TRUE, sep=',')
# check data briefly
head(fulldata)

## 2. Removing outliers using the Hampel filter
lower_bound <- median(fulldata$totalscore) - 3 * mad(fulldata$totalscore, constant = 1)
upper_bound <- median(fulldata$totalscore) + 3 * mad(fulldata$totalscore, constant = 1)
# those from organisations 7 and 9 had access to the platform but did not comply to the protocol as they reported that they did not use the digital platform
not_compliant <- (fulldata$site == 7 & fulldata$digital_platform == -1) | (fulldata$site == 9 & fulldata$digital_platform == -1)

# the final set of users are those with total score between rainge lower_bound and upper_bound and who complied to the protocol
alldata <- subset(fulldata,fulldata$totalscore >= lower_bound & fulldata$totalscore <= upper_bound & not_compliant == FALSE)

# 3. Summary of the number of participants before and after data cleaning
# before
initial_N = nrow(fulldata) # total N at beginning
initial_NT =nrow(subset(fulldata,fulldata$site ==7 | fulldata$site == 9))
initial_NTn = sum(not_compliant) # total treatment non compliant
initial_NTy = initial_NT - initial_NTn # total treatment compliant
initial_NC = initial_N - initial_NT # total controls
# after screening
final_N = nrow(alldata) # total N at beginning
final_NT =nrow(subset(alldata,alldata$digital_platform == 1))
final_NC = final_N - final_NT

print(paste0("The number of subjects who gave permission to use their data for the experiment was: ",initial_N," (",initial_NC, " controls and ",initial_NT," who received access to the platform, in percentages: ", round(100*initial_NC/initial_N,2), " C and ", round(100*initial_NT/initial_N,2), " T )"))

print(paste0("Of all the ", initial_NT, " with access to the platform ", initial_NTy, "(", round(100*initial_NTy/initial_NT,2) ,"%) complied with the protocol while ",initial_NTn,"(", round(100*initial_NTn/initial_NT,2), "%) decided to not use the platform"))

print(paste0("After removing outliers and those who did not comply with the protocol the final number of subjects was ", final_N ," with ", final_NC, "(", round(100*final_NC/final_N,2), "%) controls and ", final_NT, "(", round(100*final_NT/final_N,2), "%) intervention"    ))


## 4. Preparing the variables + descriptive statistics
# fix mapping of categorical variables
cols <- c("site", "sex", "education", "has_acq_info", "from_nurses", "from_friends_studying_in_field", "from_other_friends", "from_internet", "from_entertainment", "from_magazines", "from_supervisor", "from_event_by_uoas", "past_nurse_training", "has_worked",  "digital_platform")
alldata[cols] <- lapply(alldata[cols], factor)
# descriptive statistics for all data with outliers removed, this output is not
# used in the paper:
sumtable(alldata)

# descriptive statistics for intervention group, outliers removed
sumdata <- subset(alldata, alldata$digital_platform == 1);
sumdata[cols] <- lapply(sumdata[cols], factor)
sumtable(sumdata)

# descriptive statistics for control group, outliers removed
sumdata <- subset(alldata, alldata$digital_platform == -1);
sumdata[cols] <- lapply(sumdata[cols], factor)
sumtable(sumdata)

## 5. Test if the two groups are different in respect of some background factors
testcols <- c("age", "sex", "education", "has_acq_info", "from_nurses", "from_friends_studying_in_field", "from_other_friends", "from_internet", "from_entertainment", "from_magazines", "from_supervisor", "from_event_by_uoas", "past_nurse_training", "has_worked", "perception")

# the output is what is used for the paper descriptive statistics
testbgf(testcols,alldata)

## test if those who complied with the protocol vs those who did not comply
treatment_sites <- (fulldata$site == 7 ) | (fulldata$site == 9 )
alltreatdata <- subset(fulldata,fulldata$totalscore >= lower_bound & fulldata$totalscore <= upper_bound & treatment_sites == TRUE)
#alltreatdata <- subset(fulldata, treatment_sites == TRUE)
cols <- c("site", "sex", "education", "has_acq_info", "from_nurses", "from_friends_studying_in_field", "from_other_friends", "from_internet", "from_entertainment", "from_magazines", "from_supervisor", "from_event_by_uoas", "past_nurse_training", "has_worked",  "digital_platform")
alltreatdata[cols] <- lapply(alltreatdata[cols], factor)
testbgf(testcols,alltreatdata)

## 6. Multiple linear regression 
model_lm <- lm(totalscore ~ age + sex + education + from_nurses + from_friends_studying_in_field + from_other_friends + from_internet + from_entertainment + from_magazines + from_supervisor + from_event_by_uoas +past_nurse_training + has_worked   + perception + digital_platform, data = alldata)
model_lm
summary(model_lm)
lmqc(model_lm)

## 7. Quality control: let's use the sandwhich library to correct the variance for heteroschedasticity
# these are the MLR results that we report in the manuscript
results <- coeftest(model_lm, vcov = vcovHC(model_lm))
results
#install.packages("yhat")
library("yhat")
effect.size(model_lm)
#library(pwr)
#pwr.t.test(n=1000,d=0.17,sig.level=.05,alternative="greater")
library("parameters")
# these we also add to the manuscript:
parameters::standardize_parameters(model_lm)

library(effectsize)
#install.packages("mctest")
library("mctest")
mctestout <- mctest(model_lm,type="i",corr=TRUE)
mctestout
# compute correlation matrix, from mctest source code
x=mctestout$x
n=length(x[,1])
sx<-scale(x)/sqrt(n-1)
corR<-t(sx)%*%sx
write.csv(corR,"cov_alldata.csv")


## 8. Matched analysis
library("MatchIt")
# testing initial matching imbalance
m.out0 <- matchit(digital_platform ~ age + sex + education + from_nurses + from_friends_studying_in_field + from_other_friends + from_internet + from_entertainment + from_magazines + from_supervisor + from_event_by_uoas +past_nurse_training + has_worked   + perception, data = alldata,
                  method = NULL, distance = "glm")
summary(m.out0)


# let's match with the significant from descriptive stats
m.out1 <- matchit(digital_platform ~ age  + education  + from_nurses + from_internet + from_entertainment +  from_supervisor + from_event_by_uoas +  has_worked , data = alldata, method = "nearest", distance = "glm")

summary(m.out1)
plot(m.out1, type = "density", interactive = FALSE, which.xs = c("education"))


## let's get the subset of subjects that are matched with background factors
matcheddata <- match.data(m.out1)
sumtable(matcheddata)


## test if the new data is equal
testbgf(testcols,matcheddata)

## perform multiple linear regression using only the matched subjects
model_lm_MD <- lm(totalscore ~ age + sex + education + from_nurses + from_friends_studying_in_field + from_other_friends + from_internet + from_entertainment + from_magazines + from_supervisor + from_event_by_uoas +past_nurse_training + has_worked   + perception + digital_platform, data = matcheddata)
model_lm_MD
summary(model_lm_MD)
lmqc(model_lm_MD)
mctestout_MD <- mctest(model_lm_MD,type="i",corr=TRUE)
mctestout_MD

# compute correlation matrix, from mctest source code
x=mctestout_MD$x
n=length(x[,1])
sx<-scale(x)/sqrt(n-1)
corR_MD<-t(sx)%*%sx
write.csv(corR_MD,"cov_matcheddata.csv")
# fix results to account for heteroskedasticity 
results <- coeftest(model_lm_MD, vcov = vcovHC(model_lm_MD))
# these are the ones we report
results
parameters::standardize_parameters(model_lm_MD)
