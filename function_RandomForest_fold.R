RFcv <- function(fomula,set,foldidx,rfntree=3000){ ##mmtry,
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	model.fit <- randomForest(fomula,data=trainset,ntree= rfntree)  #mtry = mmtry,
	predY <- predict(model.fit,newdata=testset)
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testset)[testY != predY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst,mtry=model.fit$mtry))
}
