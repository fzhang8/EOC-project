library(scales)
library(class)
library(ggplot2)
library(glmnet)
library(nnet)
library(pls)
library(cvAUC)
library(randomForest)
library(e1071)
library(gbm)

load("./clean_data/clean_MS.RData")
source("./code/function_HBMBgenerate.R")
source("./code/function_kNN_fold.R")
source("./code/function_PLS_fold.R")
source("./code/function_bagging_fold.R")
source("./code/function_RandomForest_fold.R")
source("./code/function_SVM_fold.R")
source("./code/function_Boosting_fold.R")
source("./code/function_script.R")


####### record time ################
totaltime1 <- proc.time()["elapsed"]
####################################

############## Illustration ####################
#**** If it's NULL the corresponding method will not be run, except svmpolycost and svmradialcost,
#**** for which if they are null then the svmcost will be assigned to them, otherwise they will use 
#**** the specified values. All the values specified here will override the defaults that can be checked 
#**** in the "function_script.R" script.

#**** All the methods will be called from within the "function_script.R", and it is called from this script
#**** with different value settings. Each dataset is associated with one script like this including the 
#**** interested value settings.

#**** There should be three folders "clean_data","code" and "outputs" under the working directory. 
#**** The outcome will be named as (the first argument)_output.Rdata in the folder: "./outputs" where "." represents
#**** the current working directory. I put all the codes including this one into the ./code folder 
#**** and run under working directory with command: source("./code/MS_script.R").



#*** example 1 ***
#script("MS2",kNN_k=30,plscpnt=30,baggingntree=3000,rfntree=3000,svmcost=2^seq(-10,15,1),
#svmdegree=seq(2,5,1),svmgma=2^seq(-15,15,1),svmpolycost=NULL,svmradialcost=NULL,
#bstdepth=seq(1,5,1),bstshrnk=c(0.01,0.001),bstntree=c(1000,5000))

#*** example 2 ***
#script("NOESYcsvm5",kNN_k=NULL,plscpnt=NULL,baggingntree=NULL,rfntree=NULL,svmcost=2^seq(-15,40,1),
#svmdegree=seq(2,5,1),svmpolycost=2^seq(15,140,1),svmgma=2^seq(-30,30,1),svmradialcost=2^seq(-15,40,1),
#bstdepth=NULL,bstshrnk=NULL,bstntree=NULL)

#*** example 3 ***
#script("MSRF",kNN_k=NULL,plscpnt=NULL,baggingntree=NULL,rfntree=1000,svmcost=NULL,
#svmdegree=NULL,svmpolycost=NULL,svmgma=NULL,svmradialcost=NULL,
#bstdepth=NULL,bstshrnk=NULL,bstntree=NULL)



########## time used ###############
totaltime2 <- proc.time()["elapsed"]
cat("Total time used: ",(totaltime2 - totaltime1) %/% 3600,"h",round(((totaltime2 - totaltime1) %% 3600)/60,1),"m\n")
####################################
