HWE_test <- function(n2, n1, n0){
	# if(any(c(n2, n1, n0) <= 5)){
		# warning('at least one genotype count is less than 5')
		# return(NA)
	# }
	n <- n2+n1+n0
	p2 <- n2/n
	p <- (n2*2+n1)/n/2
	delta <- p2-p^2
	return(n*delta^2/p^2/(1-p)^2)
}


sdMAF_test <- function(gF, gM, snp_type, test_type){
	if(!is.numeric(gF) | length(gF) != 3){
		stop('Female genotype counts should be a numeric vector of length 3')
	}
	if(!is.numeric(gM)){
		stop('Male genotype counts should be a numeric vector')
	}
	if(snp_type == 'NPR' & length(gM) !=2){
		stop('Male genotype counts should be 2 for a NPR SNP')
	}
	if(snp_type != 'NPR' & length(gM) !=3){
		stop('Male genotype counts should be 3 for a PAR SNP')
	}
	fAA <- gF[1]; fAB <- gF[2]; fBB <- gF[3]
	nF <- sum(gF); nM <- sum(gM)
	pF <- (fAA*2+fAB)/(2*nF)
	if(snp_type == 'NPR'){
		mA <- gM[1]; mB <- gM[2]
		p <- (fAA*2 + fAB + mA)/(2*nF+nM)
		pM <- mA/nM
		}else{
			mAA <- gM[1]; mAB <- gM[2]; mBB <- gM[3]
			p <- (fAA*2 + fAB + mAA*2 + mAB)/(2*nF+2*nM)
			pM <- (mAA*2 + mAB)/(2*nM)
		}
	if(any(c(p, pF, pM) < 0.05 | c(p, pF, pM) > 0.95)){
		warning('MAF is less than 0.05')
		return(NA)
	}
	if(snp_type=='NPR'){
		score_test <- (pF-pM)^2/(1/nF/2+1/nM)/(p*(1-p))
		wald_test <- (pF-pM)^2/(pF*(1-pF)/nF/2 + pM*(1-pM)/nM)
	}else{
		score_test <- (pF-pM)^2/((1/nF/2 + 1/nM/2)*p*(1-p))
		wald_test <- (pF-pM)^2/(pF*(1-pF)/nF/2 + pM*(1-pM)/nM/2)
	}
	if(test_type=='score'){
		return(pchisq(score_test, df=1, lower.tail=F))
	}else if(test_type=='wald'){
		return(pchisq(wald_test, df=1, lower.tail=F))
	}else{
		stop('wrong test type')
	}
}


Pearson_Xchr <- function(gF, gM, snp_type, sdMAF, use_male, ECDF=NULL, size=10^7){
	if(!is.numeric(gF) | length(gF) != 3){
		stop('Female genotype counts should be a numeric vector of length 3')
	}
	if(!is.numeric(gM)){
		stop('Male genotype counts should be a numeric vector')
	}
	if(snp_type == 'NPR' & length(gM) !=2){
		stop('Male genotype counts should be 2 for a NPR SNP')
	}
	if(snp_type != 'NPR' & length(gM) !=3){
		stop('Male genotype counts should be 3 for a PAR SNP')
	}
	fAA <- gF[1]; fAB <- gF[2]; fBB <- gF[3]
	nF <- sum(gF); nM <- sum(gM)
	pF <- (fAA*2+fAB)/(2*nF)
	if(snp_type == 'NPR'){
		mA <- gM[1]; mB <- gM[2]
		p <- (fAA*2 + fAB + mA)/(2*nF+nM)
		pM <- mA/nM
		}else{
			mAA <- gM[1]; mAB <- gM[2]; mBB <- gM[3]
			p <- (fAA*2 + fAB + mAA*2 + mAB)/(2*nF+2*nM)
			pM <- (mAA*2 + mAB)/(2*nM)
		}
	if(any(c(p, pF, pM) < 0.05 | c(p, pF, pM) > 0.95)){
		warning('MAF is less than 0.05')
		return(NA)
	}
	if(snp_type=='NPR'){
		if(use_male){
			if(sdMAF){
				stop('At a NPR SNP, male sample cannot be used to test HWE if assuming sdMAF exists')
			}else{
				stat <- (mA-nM*p)^2/(nM*p) + (mB-nM*(1-p))^2/(nM*(1-p)) + (fAA-nF*p^2)^2/(nF*p^2) + (fAB-nF*2*p*(1-p))^2/(nF*2*p*(1-p)) + (fBB-nF*(1-p)^2)^2/(nF*(1-p)^2)
				return(pchisq(stat, df=2, lower.tail=F))
			}
		}else{
			if(sdMAF){
				stat <- (fAA-nF*pF^2)^2/(nF*pF^2) + (fAB-nF*2*pF*(1-pF))^2/(nF*2*pF*(1-pF)) + (fBB-nF*(1-pF)^2)^2/(nF*(1-pF)^2)
				return(pchisq(stat, df=1, lower.tail=F))
			}else{
				stat <- (fAA-nF*p^2)^2/(nF*p^2) + (fAB-nF*2*p*(1-p))^2/(nF*2*p*(1-p)) + (fBB-nF*(1-p)^2)^2/(nF*(1-p)^2)
				if(is.null(ECDF)){
					sex_ratio <- (2*nF+nM)/(nM*2)
					thres <- rchisq(n=size, df=1) + rgamma(n=size, shape=1/2, rate=sex_ratio)
					ECDF <- ecdf(thres)
				}
				return(1-ECDF(stat))
			}
		}
	}else{
		if(sdMAF){
			stat <- HWE_test(mAA, mAB, mBB) + HWE_test(fAA, fAB, fBB)
			return(pchisq(stat, df=2, lower.tail=F))
		}else{
			stat <- HWE_test(n2=(mAA+fAA), n1=(mAB+fAB), n0=(mBB+fBB))
			return(pchisq(stat, df=1, lower.tail=F))
		}
	}
}

RA_Xchr <- function(gF, gM, snp_type, sdMAF, use_male, size=10^7){
	if(!is.numeric(gF) | length(gF) != 3){
		stop('Female genotype counts should be a numeric vector of length 3')
	}
	if(!is.numeric(gM)){
		stop('Male genotype counts should be a numeric vector')
	}
	if(snp_type == 'NPR' & length(gM) !=2){
		stop('Male genotype counts should be 2 for a NPR SNP')
	}
	if(snp_type != 'NPR' & length(gM) !=3){
		stop('Male genotype counts should be 3 for a PAR SNP')
	}
	fAA <- gF[1]; fAB <- gF[2]; fBB <- gF[3]
	nF <- sum(gF); nM <- sum(gM)
	pF <- (fAA*2+fAB)/(2*nF)
	pAA <- fAA/nF
	if(snp_type == 'NPR'){
		mA <- gM[1]; mB <- gM[2]
		p <- (fAA*2 + fAB + mA)/(2*nF+nM)
		pM <- mA/nM
		}else{
			mAA <- gM[1]; mAB <- gM[2]; mBB <- gM[3]
			p <- (fAA*2 + fAB + mAA*2 + mAB)/(2*nF+2*nM)
			pM <- (mAA*2 + mAB)/(2*nM)
		}
	if(any(c(p, pF, pM) < 0.05 | c(p, pF, pM) > 0.95)){
		warning('MAF is less than 0.05')
		return(NA)
	}
	if(snp_type=='NPR'){
		if(use_male){
			if(sdMAF){
				stop('At a NPR SNP, male sample cannot be used to test HWE if assuming sdMAF exists')
			}else{
				stat <- (pF-pM)^2/(1/nF/2 + 1/nM)/p/(1-p) + nF*(pAA - pF^2 + (p-pF)^2)^2/p^2/(1-p)^2
				return(pchisq(stat, df=2, lower.tail=F))
			}
		}else{
			if(sdMAF){
				delta_f <- pAA - pF^2
				stat <- nF*delta_f^2/pF^2/(1-pF)^2
				return(pchisq(stat, df=1, lower.tail=F))
			}else{
				delta_f <- pAA - pF^2
				stat <- nF*(delta_f + (p-pF)^2)^2/p^2/(1-p)^2
				return(pchisq(stat, df=1, lower.tail=F))
			}
		}
	}else{
		if(sdMAF){
			stat <- HWE_test(mAA, mAB, mBB) + HWE_test(fAA, fAB, fBB)
			return(pchisq(stat, df=2, lower.tail=F))
		}else{
			stat <- HWE_test(n2=(mAA+fAA), n1=(mAB+fAB), n0=(mBB+fBB))
			return(pchisq(stat, df=1, lower.tail=F))
		}
	}
}
