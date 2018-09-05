logiscv <- function(fomula,set,foldidx){
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	model.fit <- multinom(fomula,data=trainset)
	predY <- predict(model.fit,newdata=testset)
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testset)[predY != testY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))
}