plscv <- function(fomula,set,foldidx,ncp){
	trainset <- set[-foldidx,]
	testset <- set[foldidx,]
	testY <- testset[,1]
	model.fit <- plsr(fomula,data=trainset,scale=FALSE,ncomp=ncp)
	predY <- predict(model.fit,newdata=testset,ncomp=ncp)
	if(length(unique(testY)) == 2){ #** HB or MB
		predY[predY >=0.5] <- 1
		predY[predY <0.5] <- 0
	}
	else{ #** Class
		predY[predY >=1.5] <- 2
		predY[predY >=0.5 & predY <1.5] <- 1
		predY[predY <0.5] <- 0
	}
	Md.test.er <- mean(predY != testY)
	missclslst <- rownames(testset)[testY != predY]
	return(list(er=Md.test.er, pred = predY,truevalue= testY,mis=missclslst))

}
