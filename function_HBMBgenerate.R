HBMBgenerate <- function(YY,md){
	if(md == "HB"){
		HB <- as.character(YY)
		HB[HB != "Cancer"] <- "Negtive"
		HB[HB == "Cancer"] <- "Positive"
		HB <- as.factor(HB)
		return(HB)
	}
	else{ #** MB
		if("Normal" %in% levels(YY)){ #** "Benign" "Cancer" "Normal" coding
			MB <- as.character(YY)
			MB[MB != "Normal"] <- "Positive"
			MB[MB == "Normal"] <- "Negtive"
			MB <- as.factor(MB)
			return(MB)
		}
		else{ #** "Benign"  "Cancer"  "Healthy" coding
			MB <- as.character(YY)
			MB[MB != "Healthy"] <- "Positive"
			MB[MB == "Healthy"] <- "Negtive"
			MB <- as.factor(MB)
			return(MB)
		}
	}
}
