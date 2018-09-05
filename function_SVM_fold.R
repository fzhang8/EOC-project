svmcv <- function(fomula,set,foldidx,cost,knel,degree=0,gamma=0){
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	if(knel == "linear"){
		model.fit <- svm(fomula,data=trainset,kernel=knel,cost=cost,scale=FALSE)
	}
	if(knel == "polynomial"){
		model.fit <- svm(fomula,data=trainset,kernel=knel,cost=cost,scale=FALSE,degree=degree)
	}
	if(knel == "radial"){
		model.fit <- svm(fomula,data=trainset,kernel=knel,cost=cost,scale=FALSE,gamma=gamma)
	}
	
	predY <- predict(model.fit,newdata=testset)
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testset)[testY != predY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))
}
