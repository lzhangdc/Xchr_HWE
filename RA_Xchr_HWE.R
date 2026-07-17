HWE_test <- function(n2, n1, n0){
	n <- n2+n1+n0
	p2 <- n2/n
	p <- (n2*2+n1)/n/2
	delta <- p2-p^2
	return(n*delta^2/p^2/(1-p)^2)
}

RA_Xchr <- function(gF, gM, snp_type, sdMAF, joint_test){
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
	if(snp_type == 'PAR' & joint_test == T){
		stop('Cannot jointly test HWE and no sdMAF for a PAR SNP')
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
		if(joint_test){
			## 2df joint test of HWD and sdMAF
			stat <- (pF-pM)^2/(1/nF/2 + 1/nM)/p/(1-p) + nF*(pAA - pF^2 + (p-pF)^2)^2/p^2/(1-p)^2
			return(pchisq(stat, df=2, lower.tail=F))
		}else{
			if(sdMAF){
				## 1 df of HWE assuming sdMAF
				delta_f <- pAA - pF^2
				stat <- nF*delta_f^2/pF^2/(1-pF)^2
				return(pchisq(stat, df=1, lower.tail=F))
			}else{
				## 1 df of HWE assuming no sdMAF
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
