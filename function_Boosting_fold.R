boostingcv <- function(fomula,set,foldidx,ntree,lamda,d){
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	model.fit <- gbm(fomula,data=trainset,n.trees=ntree,interaction.depth=d,shrinkage=lamda)
	predY <- predict(model.fit,newdata=testset,n.trees=ntree,type="response")
	if(names(set)[1] == "Class"){
		predY <- predY[,,1]
		pred <- c()
		for(i in 1:nrow(predY)){
			pred[i] <- ifelse(names(which.max(predY[i,]))=="Benign",2,ifelse(names(which.max(predY[i,]))=="Healthy",1,3))
		}
		predY <- pred
		Md.test.er <- mean(predY != as.numeric(testY))
		missclslst <- rownames(testset)[as.numeric(testY) != predY]
	}
	else{ #** HB MB
		predY[predY >= 0.5] <- 1
		predY[predY < 0.5] <- 0
		Md.test.er <- mean(predY != testY)
		missclslst <- rownames(testset)[testY != predY]
	}
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))
}
