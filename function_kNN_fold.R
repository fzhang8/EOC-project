kNNcv <- function(set,foldidx,kk){
	trainX <- as.matrix(set[-foldidx,-1])
	trainY <- as.matrix(set[-foldidx,1])
	testX <- as.matrix(set[foldidx,-1])
	testY <- as.matrix(set[foldidx,1])
	predY <- knn(trainX,testX,trainY,k=kk)
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testY)[predY != testY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))
}