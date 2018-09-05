baggingcv <- function(fomula,set,foldidx,baggingntree=3000){
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	model.fit <- randomForest(fomula,data=trainset,mtry = ncol(set)-1,ntree= baggingntree)
	predY <- predict(model.fit,newdata=testset)
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testset)[testY != predY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))
}
