
################## script function for all MS and NMR data ###############
set.seed(12345678)

script <- function(setname,kNN_k=15,plscpnt=15,baggingntree=3000,rfntree=3000,svmcost=2^seq(-10,15,1),
										svmdegree=seq(2,3,1),svmpolycost=NULL,svmgma=2^seq(-15,15,1),svmradialcost=NULL,
										bstdepth=seq(1,3,1),bstshrnk=c(0.01,0.001),bstntree=c(10,100)){

MB <- HBMBgenerate(YY,"MB") #** factor
HB <- HBMBgenerate(YY,"HB")	#** factor
if("Healthy" %in% levels(YY)){
	Class <- factor(YY,levels = c("Healthy","Benign","Cancer"))	#** factor "Healthy"-1,"Benign"-2,"Cancer"-3
}else{
	Class <- factor(YY,levels = c("Normal","Benign","Cancer"))	#** factor "Normal"-1,"Benign"-2,"Cancer"-3
}

folds <- 3
nobs <- nrow(XX)
p <- ncol(XX)

########## fold generate ##########
cvidx <- sample(1:nobs)  #** generate a permutation from 1 to Number of obs
fold1 <- cvidx[1:floor(nobs/folds)]  #** the first one third of the permutation will be fold 1
fold2 <- cvidx[(floor(nobs/folds)+1):floor(2*nobs/folds)] #** the second one third will be fold 2
fold3 <- cvidx[(floor(2*nobs/folds)+1):nobs] #** the last one third part will be fold 3

######### kNN ############
if(is.null(kNN_k)){kNNoutput <- NULL}else{

kNNoutput <- list()
cat("kNN","\n")
kvalues <- seq(1,kNN_k,2) ##** chosen k's 33
for(md in c("HB","MB","Class")){ ##** for loop to deal with all 3 models
	accuracy <- data.frame()
	roc <- list()
	missclsfd <- list()
	dataset <- cbind(get(md),scale(XX))
	for(kk in kvalues){ ##** for loop to deal with all potential k
		cat("k:",kk,"\n")
		for(i in 1:folds){
			tempout <- kNNcv(dataset,get(paste("fold",i,sep="")),kk)
			assign(paste("folderr_",i,sep=""), tempout$er) ##** test error of each fold
			assign(paste("foldpred_",i,sep=""), tempout$pred) ##** prediction of each fold
			assign(paste("foldtrue_",i,sep=""), tempout$truevalue) ##** true values of each fold
			assign(paste("foldmis_",i,sep=""), tempout$mis) ##** misclassified id for each fold
		}
		assign(paste("test.accuracy.",kk,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
		assign(paste("test.accuracy.fold1",kk,sep=""),1 - folderr_1)
		assign(paste("test.accuracy.fold2",kk,sep=""),1 - folderr_2)
		assign(paste("test.accuracy.fold3",kk,sep=""),1 - folderr_3)
		accuracy <- rbind(accuracy,data.frame(K=kk,Accuracy=get(paste("test.accuracy.",kk,sep="")),
																							Accuracy_fd1=get(paste("test.accuracy.fold1",kk,sep="")),
																							Accuracy_fd2=get(paste("test.accuracy.fold2",kk,sep="")),
																							Accuracy_fd3=get(paste("test.accuracy.fold3",kk,sep="")),
																						Accuracy_sd=sd(c(get(paste("test.accuracy.fold1",kk,sep="")),
																														get(paste("test.accuracy.fold2",kk,sep="")),
																														get(paste("test.accuracy.fold3",kk,sep=""))))))
		pred <- list(Foldpred_1=foldpred_1,Foldpred_2=foldpred_2,Foldpred_3=foldpred_3)
		True <- list(Foldtrue_1=foldtrue_1[,1],Foldtrue_2=foldtrue_2[,1],Foldtrue_3=foldtrue_3[,1])
		roc[[paste("k",kk,sep="")]] <- list(Prediction=pred,Truelabel=True)
		missclsfd[[paste("k",kk,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
	}
	assign(paste(md,"_ouput",sep=""),list(Accuracy=accuracy,ROC=roc,Missclassified=missclsfd)) ##** each model's results are included in one list
	kNNoutput[md] <- list(get(paste(md,"_ouput",sep=""))) ##** all 3 models' results are included in one list
}
}

############ PLS ###########
if(is.null(plscpnt)){pls_output <- NULL}else{

pls_output <- list()
cat("PLS","\n")
components <- seq(1,plscpnt,1) #30
for(md in c("HB","MB","Class")){ ##** for loop to deal with all 3 models
	accuracy <- data.frame()
	roc <- list()
	missclsfd <- list()
	dataset <- data.frame(as.numeric(get(md))-1,scale(XX))
	colnames(dataset)[1] <- md
	fmula <- as.formula(paste(md,"~."))
	for(cp in components){
		cat("Cp",cp,"\n")
		for(i in 1:folds){
			tempout <- plscv(fmula,dataset,get(paste("fold",i,sep="")),cp)
			assign(paste("folderr_",i,sep=""), tempout$er) 
			assign(paste("foldpred_",i,sep=""), tempout$pred) 
			assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
			assign(paste("foldmis_",i,sep=""), tempout$mis) 
		}
		assign(paste("test.accuracy.",cp,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
		assign(paste("test.accuracy.fold1",cp,sep=""),1 - folderr_1)
		assign(paste("test.accuracy.fold2",cp,sep=""),1 - folderr_2)
		assign(paste("test.accuracy.fold3",cp,sep=""),1 - folderr_3)
		accuracy <- rbind(accuracy,data.frame(NumberLV=cp,Accuracy=get(paste("test.accuracy.",cp,sep="")),
																							Accuracy_fd1=get(paste("test.accuracy.fold1",cp,sep="")),
																							Accuracy_fd2=get(paste("test.accuracy.fold2",cp,sep="")),
																							Accuracy_fd3=get(paste("test.accuracy.fold3",cp,sep="")),
																							Accuracy_sd=sd(c(get(paste("test.accuracy.fold1",cp,sep="")),
																														get(paste("test.accuracy.fold2",cp,sep="")),
																														get(paste("test.accuracy.fold3",cp,sep=""))))))
		pred <- list(Foldpred_1=foldpred_1[,,1],Foldpred_2=foldpred_2[,,1],Foldpred_3=foldpred_3[,,1])
		True <- list(Foldtrue_1=foldtrue_1,Foldtrue_2=foldtrue_2,Foldtrue_3=foldtrue_3)
		roc[[paste("Cp_",cp,sep="")]] <- list(Prediction=pred,Truelabel=True)
		missclsfd[[paste("Cp_",cp,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
	}
	assign(paste(md,"_ouput",sep=""),list(Accuracy=accuracy,ROC=roc,Missclassified=missclsfd)) 
	pls_output[md] <- list(get(paste(md,"_ouput",sep=""))) 
}
}

############ Bagging ################
if(is.null(baggingntree)){bagging_output <- NULL}else{

bagging_output <- list()
cat("bagging","\n")
for(md in c("HB","MB","Class")){ 
	dataset <- data.frame(get(md),XX)
	colnames(dataset)[1] <- md
	fmula <- as.formula(paste(md,"~."))
	for(i in 1:folds){
		tempout <- baggingcv(fmula,dataset,get(paste("fold",i,sep="")),baggingntree)
		assign(paste("folderr_",i,sep=""), tempout$er) 
		assign(paste("foldpred_",i,sep=""), tempout$pred) 
		assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
		assign(paste("foldmis_",i,sep=""), tempout$mis) 
	}
	accuracy <- 1 - (folderr_1 + folderr_2 + folderr_3)/folds
	accuracy.fold1 <- 1 - folderr_1
	accuracy.fold2 <- 1 - folderr_2
	accuracy.fold3 <- 1 - folderr_3
	pred <- list(Foldpred_1=foldpred_1,Foldpred_2=foldpred_2,Foldpred_3=foldpred_3)
	True <- list(Foldtrue_1=foldtrue_1,Foldtrue_2=foldtrue_2,Foldtrue_3=foldtrue_3)
	roc <- list(Prediction=pred,Truelabel=True)
	missclsfd <- c(foldmis_1,foldmis_2,foldmis_3)
	assign(paste(md,"_ouput",sep=""),list(Accuracy=data.frame(Accuracy=accuracy,Accuracy_fd1=accuracy.fold1,
	Accuracy_fd2=accuracy.fold2,Accuracy_fd3=accuracy.fold3,
	Accuracy_sd=sd(c(accuracy.fold1,accuracy.fold2,accuracy.fold3))),ROC=roc,Missclassified=missclsfd)) 
	bagging_output[md] <- list(get(paste(md,"_ouput",sep=""))) 
}
}

########### RandomForest #############
if(is.null(rfntree)){RF_output <- NULL}else{

RF_output <- list()
cat("RF","\n")
#variabletrys <- floor(sqrt(p)) ##** mtry sequence -10,10
for(md in c("HB","MB","Class")){ 
	accuracy <- data.frame()
	roc <- list()
	missclsfd <- list()
	dataset <- data.frame(get(md),XX)
	colnames(dataset)[1] <- md
	fmula <- as.formula(paste(md,"~."))
#	for(mmtry in variabletrys){
#		cat("mtry",mmtry,"\n")
		for(i in 1:folds){
			tempout <- RFcv(fmula,dataset,get(paste("fold",i,sep="")),rfntree) ##mmtry,
			assign(paste("folderr_",i,sep=""), tempout$er) 
			assign(paste("foldpred_",i,sep=""), tempout$pred) 
			assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
			assign(paste("foldmis_",i,sep=""), tempout$mis) 
			mmtry <- tempout$mtry
		}
		cat(md,"mtry",mmtry,"\n")
		assign(paste("test.accuracy.",mmtry,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
		assign(paste("test.accuracy.fold1",mmtry,sep=""),1 - folderr_1)
		assign(paste("test.accuracy.fold2",mmtry,sep=""),1 - folderr_2)
		assign(paste("test.accuracy.fold3",mmtry,sep=""),1 - folderr_3)
		
		accuracy <- rbind(accuracy,data.frame(SplitTrys=mmtry,Accuracy=get(paste("test.accuracy.",mmtry,sep="")),
																					Accuracy_fd1=get(paste("test.accuracy.fold1",mmtry,sep="")),
																					Accuracy_fd2=get(paste("test.accuracy.fold2",mmtry,sep="")),
																					Accuracy_fd3=get(paste("test.accuracy.fold3",mmtry,sep="")),
																					Accuracy_sd=sd(c(get(paste("test.accuracy.fold1",mmtry,sep="")),
																												get(paste("test.accuracy.fold2",mmtry,sep="")),
																												get(paste("test.accuracy.fold3",mmtry,sep=""))))))
		pred <- list(Foldpred_1=foldpred_1,Foldpred_2=foldpred_2,Foldpred_3=foldpred_3)
		True <- list(Foldtrue_1=foldtrue_1,Foldtrue_2=foldtrue_2,Foldtrue_3=foldtrue_3)
		roc[[paste("mtry_",mmtry,sep="")]] <- list(Prediction=pred,Truelabel=True)
		missclsfd[[paste("mtry_",mmtry,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
#	}
	assign(paste(md,"_ouput",sep=""),list(Accuracy=accuracy,ROC=roc,Missclassified=missclsfd)) 
	RF_output[md] <- list(get(paste(md,"_ouput",sep=""))) 
}
}

############### SVM ####################
if(is.null(svmcost)){svm_output <- NULL}else{

svm_output <- list()
cat("svm","\n")
costs <- svmcost #-10,15
if(is.null(svmpolycost)){svmpolycosts <- svmcost}else{svmpolycosts <- svmpolycost}
if(is.null(svmradialcost)){svmradialcosts <- svmcost}else{svmradialcosts <- svmradialcost}

degrees <- svmdegree
gammas <- svmgma #-15,15
for(md in c("HB","MB")){
	for(knel in c("linear","polynomial","radial")){
		cat(knel,"\n")
		accuracy <- data.frame()
		missclsfd <- list()
		dataset <- data.frame(get(md),XX)
		colnames(dataset)[1] <- md
		fmula <- as.formula(paste(md,"~."))
		
		
			if(knel == "linear"){
				for(cost in costs){
			
				for(i in 1:folds){
					tempout <- svmcv(fmula,dataset,get(paste("fold",i,sep="")),cost,knel)
					assign(paste("folderr_",i,sep=""), tempout$er) 
					assign(paste("foldpred_",i,sep=""), tempout$pred) 
					assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
					assign(paste("foldmis_",i,sep=""), tempout$mis) 
				}
				assign(paste("test.accuracy.",cost,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
				assign(paste("test.accuracy.fold1",cost,sep=""),1 - folderr_1)
				assign(paste("test.accuracy.fold2",cost,sep=""),1 - folderr_2)
				assign(paste("test.accuracy.fold3",cost,sep=""),1 - folderr_3)
				
				accuracy <- rbind(accuracy,data.frame(Cost=cost,Accuracy=get(paste("test.accuracy.",cost,sep="")),
																					Accuracy_fd1=get(paste("test.accuracy.fold1",cost,sep="")),
																					Accuracy_fd2=get(paste("test.accuracy.fold2",cost,sep="")),
																					Accuracy_fd3=get(paste("test.accuracy.fold3",cost,sep="")),
																					Accuracy_sd=sd(c(get(paste("test.accuracy.fold1",cost,sep="")),
																												get(paste("test.accuracy.fold2",cost,sep="")),
																												get(paste("test.accuracy.fold3",cost,sep=""))))))
				missclsfd[[paste("Cost_",cost,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
				cat(md,knel,"cost:",cost,"\n")
			
				}
			}
			
			
			if(knel == "polynomial"){
				for(cost in svmpolycosts){
				
				for(degree in degrees){
					for(i in 1:folds){
						tempout <- svmcv(fmula,dataset,get(paste("fold",i,sep="")),cost,knel,degree)
						assign(paste("folderr_",i,sep=""), tempout$er) 
						assign(paste("foldpred_",i,sep=""), tempout$pred) 
						assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
						assign(paste("foldmis_",i,sep=""), tempout$mis) 
					}
					assign(paste("test.accuracy.cost.",cost,"degree",degree,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
					assign(paste("test.accuracy.cost.fold1",cost,"degree",degree,sep=""),1 - folderr_1)
					assign(paste("test.accuracy.cost.fold2",cost,"degree",degree,sep=""),1 - folderr_2)
					assign(paste("test.accuracy.cost.fold3",cost,"degree",degree,sep=""),1 - folderr_3)
					accuracy <- rbind(accuracy,data.frame(Cost=cost,Degree=degree,Accuracy=get(paste("test.accuracy.cost.",cost,"degree",degree,sep="")),
												Accuracy_fd1=get(paste("test.accuracy.cost.fold1",cost,"degree",degree,sep="")),
												Accuracy_fd2=get(paste("test.accuracy.cost.fold2",cost,"degree",degree,sep="")),
												Accuracy_fd3=get(paste("test.accuracy.cost.fold3",cost,"degree",degree,sep="")),
												Accuracy_sd=sd(c(get(paste("test.accuracy.cost.fold1",cost,"degree",degree,sep="")),
																			get(paste("test.accuracy.cost.fold2",cost,"degree",degree,sep="")),
																		get(paste("test.accuracy.cost.fold3",cost,"degree",degree,sep=""))))))
					missclsfd[[paste("Cost_",cost,"degree_",degree,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
					cat(md,knel,"cost:",cost,"degree:",degree,"\n")
				}
				}
			}
			
			
			if(knel == "radial"){
				for(cost in svmradialcosts){
			
				for(gamma in gammas){
					for(i in 1:folds){
						tempout <- svmcv(fmula,dataset,get(paste("fold",i,sep="")),cost,knel,,gamma)
						assign(paste("folderr_",i,sep=""), tempout$er) 
						assign(paste("foldpred_",i,sep=""), tempout$pred) 
						assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
						assign(paste("foldmis_",i,sep=""), tempout$mis) 
					}
					assign(paste("test.accuracy.cost.",cost,"gma",gamma,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
					assign(paste("test.accuracy.cost.fold1",cost,"gma",gamma,sep=""),1 - folderr_1)
					assign(paste("test.accuracy.cost.fold2",cost,"gma",gamma,sep=""),1 - folderr_2)
					assign(paste("test.accuracy.cost.fold3",cost,"gma",gamma,sep=""),1 - folderr_3)
					accuracy <- rbind(accuracy,data.frame(Cost=cost,Gamma=gamma,Accuracy=get(paste("test.accuracy.cost.",cost,"gma",gamma,sep="")),
												Accuracy_fd1=get(paste("test.accuracy.cost.fold1",cost,"gma",gamma,sep="")),
												Accuracy_fd2=get(paste("test.accuracy.cost.fold2",cost,"gma",gamma,sep="")),
												Accuracy_fd3=get(paste("test.accuracy.cost.fold3",cost,"gma",gamma,sep="")),
												Accuracy_sd=sd(c(get(paste("test.accuracy.cost.fold1",cost,"gma",gamma,sep="")),
																				get(paste("test.accuracy.cost.fold2",cost,"gma",gamma,sep="")),
																				get(paste("test.accuracy.cost.fold3",cost,"gma",gamma,sep=""))))))
					missclsfd[[paste("Cost_",cost,"gma_",gamma,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
					cat(md,knel,"cost:",cost,"gma",gamma,"\n")
				}
			
			}
			}
			
			
		
		assign(paste(knel,"_output",sep=""),list(Accuracy=accuracy,Missclassified=missclsfd))
	}
	assign(paste(md,"_ouput",sep=""),list(Linear=linear_output,Polynomial=polynomial_output,Radial=radial_output))
	svm_output[md] <- list(get(paste(md,"_ouput",sep=""))) 
}
}

########### Boosting #############
if(is.null(bstdepth)){boosting_output <- NULL}else{

boosting_output <- list()
cat("Boosting","\n")
depth <- bstdepth
shrinkage <- bstshrnk
treenum <- bstntree #1000,5000
for(md in c("HB","MB","Class")){
	accuracy <- data.frame()
	missclsfd <- list()
	if(md=="Class"){dataset <- data.frame(get(md),XX)}else{dataset <- data.frame(as.numeric(get(md))-1,XX)}
	colnames(dataset)[1] <- md
	fmula <- as.formula(paste(md,"~."))
	for(ntree in treenum){
		for(lamda in shrinkage){
			for(d in depth){
			cat(md,ntree,lamda,d,"\n")
					for(i in 1:folds){
						time3 <- proc.time()["elapsed"]
						tempout <- boostingcv(fmula,dataset,get(paste("fold",i,sep="")),ntree,lamda,d)
						time4 <- proc.time()["elapsed"]
						cat("Time used for boosting of fold_",i,", ",(time4 - time3) %/% 3600,"h",round(((time4 - time3) %% 3600)/60,1),"m\n")
						assign(paste("folderr_",i,sep=""), tempout$er) 
						assign(paste("foldpred_",i,sep=""), tempout$pred) 
						assign(paste("foldtrue_",i,sep=""), tempout$truevalue) 
						assign(paste("foldmis_",i,sep=""), tempout$mis) 
					}
					assign(paste("test.accuracy.ntree.",ntree,"lamda",lamda,"d",d,sep=""),1 - (folderr_1 + folderr_2 + folderr_3)/folds)
					assign(paste("test.accuracy.ntree.fold1",ntree,"lamda",lamda,"d",d,sep=""),1 - folderr_1)
					assign(paste("test.accuracy.ntree.fold2",ntree,"lamda",lamda,"d",d,sep=""),1 - folderr_2)
					assign(paste("test.accuracy.ntree.fold3",ntree,"lamda",lamda,"d",d,sep=""),1 - folderr_3)
					accuracy <- rbind(accuracy,data.frame(TreeNum=ntree,Shrinkage=lamda,Depth=d,Accuracy=get(paste("test.accuracy.ntree.",ntree,"lamda",lamda,"d",d,sep="")),
									Accuracy_fd1=get(paste("test.accuracy.ntree.fold1",ntree,"lamda",lamda,"d",d,sep="")),
									Accuracy_fd2=get(paste("test.accuracy.ntree.fold2",ntree,"lamda",lamda,"d",d,sep="")),
									Accuracy_fd3=get(paste("test.accuracy.ntree.fold3",ntree,"lamda",lamda,"d",d,sep="")),
									Accuracy_sd=sd(c(get(paste("test.accuracy.ntree.fold1",ntree,"lamda",lamda,"d",d,sep="")),
															get(paste("test.accuracy.ntree.fold2",ntree,"lamda",lamda,"d",d,sep="")),
															get(paste("test.accuracy.ntree.fold3",ntree,"lamda",lamda,"d",d,sep=""))))))
					missclsfd[[paste("TreeNum.",ntree,sep="")]][[paste("Shrinkage.",lamda,sep="")]][[paste("Depth.",d,sep="")]] <- c(foldmis_1,foldmis_2,foldmis_3)
			}
		}
	}
	assign(paste(md,"_ouput",sep=""),list(Accuracy=accuracy,Missclassified=missclsfd))
	boosting_output[md] <- list(get(paste(md,"_ouput",sep=""))) 
}
}

########## save all ############
assign(paste(setname,"_output",sep=""),list(kNN=kNNoutput,PLS=pls_output,Bagging=bagging_output,
																RandomForest=RF_output,	SVM=svm_output,Boosting=boosting_output))
save(list=paste(setname,"_output",sep=""),file=paste("./outputs/",setname,"_output.Rdata",sep=""))

}