#####################################################
#  Spatially variable Sex-Specific Selection and the
#  maintenance of SA polymorphism in simultaneous
#  hermaphrodites
#
#  Necessary functions for simulations and analyses
#
#  Author: Colin Olito
#
#  NOTES:  
#          

###############
# Dependencies



##########################################
##########################################
# 2-Patch Levene Model for Hermaphrodites

#' General lambda for q = 0
#'
#' @title General lambda for q = 0
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function.
#' @param hm     Dominance coefficient for male sex function
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
lambda0  <-  function(C, delta, hf, hm, sf, sm) {
	(2*(-2 + sm + hm*sm + hf*(sf - sf*sm)) + C*(2 + sf + sm - 4*hm*sm - sf*sm + 
    	4*(-1 + hf*sf)*(-1 + sm)*delta) + (C^2)*((-1 + 2*hm)*sm + 
    		2*(-1 + sm)*delta - (-1 + 2*hf)*sf*(-1 + sm)*(-1 + 2*delta))) / 
				(2*(-2 + C)*(-1 + sm)*(-1 + C*delta))
}
#' General lambda for q = 1
#'
#' @title General lambda for q = 0
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function.
#' @param hm     Dominance coefficient for male sex function
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
lambda1  <-  function(C, delta, hf, hm, sf, sm) {
	((-1 + C)*(-1 + sf)*(2 - 2*hm*sm + C*(-1 + (-1 + 2*hm)*sm)) - 
		(2 - C + 2*(-1 + C)*hf*sf)*(-1 + C*(-1 + 2*delta))) / 
			(2*(2 - C)*(-1 + sf)*(-1 + C*delta))
}

#' 1-Patch Invasion conditions for q = 0
#'
#' @title General lambda for q = 0
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function.
#' @param hm     Dominance coefficient for male sex function
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
InvB  <-  function(C, delta, hf, hm, sm) {
	((C - 1)*(2 - C + 2*(C - 1)*hm)*sm) / 
		((2*(C - 1)*hf - C)*(sm - 1)*(-1 + C*(2*delta - 1)))
}

#' 1-Patch Invasion conditions for q = 1
#'
#' @title General lambda for q = 1
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function.
#' @param hm     Dominance coefficient for male sex function
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
InvA  <-  function(C, delta, hf, hm, sm) {
	((C - 1)*(2*(C - 1)*hm - C)*sm) / 
		(2 - 2*hf + 2*hm*sm + (C^2)*(-1 - sm + 2*hm*sm + hf*(2 - 4*delta) + 
			2*delta) +  C*(1 + sm - 4*hm*sm + 4*(hf - 1)*delta))
}


#' Single-locus SA equilibrium allele frequencies (Additive SA fitness effects)
#'
#' @title Single-locus SA equilibrium allele frequencies
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
qHatAdd  <-  function(C, delta, sf, sm) {
	qHat  <-  (sf - sm + sf*sm + C*(sf + sm - sf*sm - 2*sf*delta)) /(2*(sf*sm - C*sf*sm*delta))
	if(qHat < 0 | is.nan(qHat) | is.na(qHat))
		qHat  <-  0
	if(qHat > 1)
		qHat  <-  1
	qHat
}

#' Single-locus SA equilibrium allele frequencies (partially recessive SA fitness effects; dominance reversals)
#'
#' @title Single-locus SA equilibrium allele frequencies
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param sf     Selection coefficient for female sex function
#' @param sm     Selection coefficient for male sex function
#' @export
qHatDomRev  <-  function(C, delta, sf, sm, h) {
	qHat  <-  ((-1 + C)*sm*(-2*h + C*(-1 + 2*h + delta)) + sf*(2 - 2*h + C*(-1 + 2*h - delta))*(-1 + C*(-1 + 2*delta))) / (2*(-1 + C)*(-1 + 2*h)*((-1 + C)*sm + sf*(-1 + C*(-1 + 2*delta))))
	if(qHat < 0 | is.nan(qHat) | is.na(qHat))
		qHat  <-  0
	if(qHat > 1)
		qHat  <-  1
	qHat
}



#' Convenience function for calculating covariances for general 
#' invasion conditions under additive fitness
#'
#' @title Calculate covariances for general invasion conditions
#' @param svals	Matrix of seleciton values for female and male fitness
#' @export 
matCov  <-  function(svals) {
	i    <- length(svals)/2
	sfs  <-  svals[1:i]
	sms  <-  svals[(i+1):(2*i)]
	cov(sfs,sms)
}


#' General sex-specific multipatch invasion conditions based on additive fitness
#'
#' @title General lambda for q = 0
#' @param q      Frequency of a allele (determines which boundary 
#' 				 to do stability analysis for)
#' @param k      Number of patches (integer >= 1)
#' @param C      Population selfing rate
#' @param sf     matrix of selection coefficients for female sex function across k patches
#' @param sm     matrix of selection coefficients for male sex function across k patches
#' @export
generalAddInv  <-  function(q, k, C, sfs, sms) {
	
	if(q != 0 & q != 1) {
		stop('q must equal 0 or 1')
	}
	if(k < 1 | k != ncol(sfs) | k != ncol(sms)) {
		stop('k must equal ncol(sfs) = ncol(sms)')
	}

	# storage structures
	kTot  <-  ncol(sfs)
	n     <-  nrow(sfs)
	Wf    <-  array(0, dim=c(n, 3, ncol(sfs)))
	Wm    <-  array(0, dim=c(n, 3, ncol(sfs)))

	# Calculate fitness through each sex function (additive fitness)
	for(i in 1:kTot) {
		x  <-  Wf.fit(hf = 1/2, sf=sfs[,i])[-1]
		Wf[,1,i]  <-  1
		Wf[,2,i]  <-  x[1:n]
		Wf[,3,i]  <-  x[(n+1):(2*n)]
		x  <-  Wm.fit(hm = 1/2, sm=sms[,i])
		Wm[,1,i]  <-  x[1:n]
		Wm[,2,i]  <-  x[(n+1):(2*n)]
		Wm[,3,i]  <-  1
	}

	# Calculate selection coefficients
	sfMat     <-  matrix(0,ncol=kTot, nrow=n)
	smMat     <-  matrix(0,ncol=kTot, nrow=n)
	Lambdas   <-  matrix(0,ncol=kTot, nrow=n)

	# Loop over number of patches
	for(i in 1:kTot) {
		if(q == 0) {
			sfMat[,i]  <-  (Wf[,2,i] - Wf[,1,i]) / Wf[,1,i]
			smMat[,i]  <-  (Wm[,2,i] - Wm[,1,i]) / Wm[,1,i]
		}
		if(q == 1) {
			sfMat[,i]  <-  (Wf[,2,i] - Wf[,3,i]) / Wf[,3,i]
			smMat[,i]  <-  (Wm[,2,i] - Wm[,3,i]) / Wm[,3,i]
		}
	
		# If there is only one patch
		if(i == 1) {
			sfBar        <-  sfMat[,i]
			smBar        <-  smMat[,i]
			Lambdas[,i]  <-  ((1 + C)/(2 - C))*sfBar + ((1 - C)/(2 - C))*smBar
		}

		# For multiple patches
		if(i > 1) {
			sfBar  <-  rowMeans(sfMat[,1:i])
			smBar  <-  rowMeans(smMat[,1:i])
			sfVar  <-  apply(sfMat[,1:i], MARGIN=1, var)
			smVar  <-  apply(smMat[,1:i], MARGIN=1, var)
			CoVar  <-  apply(cbind(sfMat[,1:i],smMat[,1:i]), MARGIN=1, FUN = matCov)
	
			Lambdas[,i]  <-  ((1 + C)/(2 - C))*sfBar + ((1 + C)^2/(2 - C)^2)*sfVar +
							((1 - C)/(2 - C))*smBar + ((1 - C)^2/(2 - C)^2)*smVar +
							(2*(1 - C^2)/(2 - C)^2)*CoVar
		}
	}

	Lambdas[Lambdas > 0]   <-  1
	Lambdas[Lambdas <= 0]  <-  0
	Lambdas
}



#' SA invasion conditions, optimized for matrix input
#'
#' @title General lambda for q = 0
#' @param q      Frequency of a allele (determines which boundary 
#' 				 to do stability analysis for)
#' @param k      Number of patches (integer >= 1)
#' @param C      Population selfing rate
#' @param sf     matrix of selection coefficients for female sex function across k patches
#' @param sm     matrix of selection coefficients for male sex function across k patches
#' @export
SAAddInv  <-  function(q, k, C, sfs, sms) {
	
	if(q != 0 & q != 1) {
		stop('q must equal 0 or 1')
	}
	if(k < 1 | k != ncol(sfs) | k != ncol(sms)) {
		stop('k must equal ncol(sfs) = ncol(sms)')
	}

	# storage structures
	kTot  <-  ncol(sfs)
	n     <-  nrow(sfs)

	# Calculate selection coefficients
	Lambdas        <-  matrix(0,ncol=kTot, nrow=n)
	kPatchLambdas  <-  matrix(0,ncol=kTot, nrow=n)

	# Loop over number of patches
	for(i in 1:kTot) {
		if(q == 0) {
			Lambdas[,i]        <-  lambda0(C=C, delta=0, hf=0.5, hm=0.5, sf=sfs[,i], sm=sms[,i])
		}
		if(q == 1) {
			Lambdas[,i]        <-  lambda1(C=C, delta=0, hf=0.5, hm=0.5, sf=sfs[,i], sm=sms[,i])
		}
		if(i == 1) {
			kPatchLambdas[,i]  <-  Lambdas[,i]
		}
		if(i > 1) {
			kPatchLambdas[,i]  <-  rowMeans(Lambdas[,1:i])
		}	
	}
	kPatchLambdas[kPatchLambdas <= 1]  <-  0
	kPatchLambdas[kPatchLambdas > 1]   <-  1
	kPatchLambdas
}


###############################################################
###############################################################
#' Simulation to compare proportion of parameter space where
#' polymorphism is maintained in the 2-patch model compared
#' with the 1-patch expectation 
#'
#' @title 2-Patch SA in hermaphrodites Simulation
#' @param n      Sample size (number of selection coefficients 
#' 				 randomly drawn from uniform distribution)
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function. Assume equal across patches.
#' @param hm     Dominance coefficient for male sex function. Assume equal across patches.
#' @param sMax   Maximum selection coefficient in Patch 1 (determines range
#' 				  of selection coefficient parameter space to be explored)
#' 				  (we assume that we always explore a square parameter space 
#' 				  (i.e., sMax is the same for males and females, and equal 
#' 				  across patches))
#' @export
simMultiPatch  <-  function(n, C, delta, hf, hm, sMax) {

	# Draw random seleciton coefficients for up to 5 patches
	sf1  <-  runif(n, max = sMax)
	sm1  <-  runif(n, max = sMax)
	sf2  <-  runif(n, max = sMax)
	sm2  <-  runif(n, max = sMax)
	sf3  <-  runif(n, max = sMax)
	sm3  <-  runif(n, max = sMax)
	sf4  <-  runif(n, max = sMax)
	sm4  <-  runif(n, max = sMax)
	sf5  <-  runif(n, max = sMax)
	sm5  <-  runif(n, max = sMax)

	# 1-Patch invasion criteria for boundaries of q = 0 and q = 1
	lambda0_1patch  <-  lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1)
	lambda1_1patch  <-  lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1)

	# 2-Patch invasion criteria for boundaries of q = 0 and q = 1
	lambda0_2patch  <-  (1/2)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2))
	lambda1_2patch  <-  (1/2)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2))

	# 3-Patch invasion criteria for boundaries of q = 0 and q = 1
	lambda0_3patch  <-  (1/3)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3))
	lambda1_3patch  <-  (1/3)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3))

	# 4-Patch invasion criteria for boundaries of q = 0 and q = 1
	lambda0_4patch  <-  (1/4)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf4, sm = sm4))
	lambda1_4patch  <-  (1/4)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf4, sm = sm4))

	# 5-Patch invasion criteria for boundaries of q = 0 and q = 1
	lambda0_5patch  <-  (1/5)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf4, sm = sm4) + 
							   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf5, sm = sm5))
	lambda1_5patch  <-  (1/5)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf1, sm = sm1) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf2, sm = sm2) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf3, sm = sm3) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf4, sm = sm4) + 
							   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf5, sm = sm5))


	# Calculate proportion of parameter space 
	# where polymorphism is predicted to be 
	# maintained by balancing selection
	poly1Patch  <-  sum(1 < lambda0_1patch & 1 < lambda1_1patch)/n
	poly2Patch  <-  sum(1 < lambda0_2patch & 1 < lambda1_2patch)/n
	poly3Patch  <-  sum(1 < lambda0_3patch & 1 < lambda1_3patch)/n
	poly4Patch  <-  sum(1 < lambda0_4patch & 1 < lambda1_4patch)/n
	poly5Patch  <-  sum(1 < lambda0_5patch & 1 < lambda1_5patch)/n

	# Save and return results
	res  <-  list(
				  "poly1"  =  poly1Patch,
				  "poly2"  =  poly2Patch,
				  "poly3"  =  poly3Patch,
				  "poly4"  =  poly4Patch,
				  "poly5"  =  poly5Patch
				  )
	return(res)
}




#' Simulation to compare proportion of parameter space where
#' polymorphism is maintained in the 2-patch model compared
#' with the 1-patch expectation across a gradient of selection
#' coefficients (wrapper function for sim2Patch())
#'
#' @title 2-Patch SA in hermaphrodites Simulation
#' @param n      Sample size (number of selection coefficients 
#' 				 randomly drawn from uniform distribution)
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function. Assume equal across patches.
#' @param hm     Dominance coefficient for male sex function. Assume equal across patches.
#' @param sMax   Maximum selection coefficient in Patch 1 (determines range
#' 				  of selection coefficient parameter space to be explored)
#' 				  (we assume that we always explore a square parameter space 
#' 				  (i.e., sMax is the same for males and females, and equal 
#' 				  across patches))
#' @param resolution resolution for selection coefficient gradient
#' @export
simMultiPatchSgrad  <-  function(n, C, delta, hf, hm, sMax = 1, resolution = 0.01) {

	# Initialize storage structures
	Poly1  <-  c()
	Poly2  <-  c()
	Poly3  <-  c()
	Poly4  <-  c()
	Poly5  <-  c()

	# Create vector of sMaxes (we assume that we always explore a 
	# square parameter space (i.e., sMax is the same for males and females))
	sMaxes  <- seq(from = resolution, to = sMax, by = resolution)
	
	# loop over sMaxes
	for (i in 1:length(sMaxes)) {
		res  <-  simMultiPatch(n = n, C = C, delta = delta, hf = hf, hm = hm, sMax = sMaxes[i]) 
		Poly1[i]  <-  res$poly1
		Poly2[i]  <-  res$poly2
		Poly3[i]  <-  res$poly3
		Poly4[i]  <-  res$poly4
		Poly5[i]  <-  res$poly5
	}

	# Save results as data frame
	data  <-  data.frame(
						 "sMax"      =  sMaxes,
						 "poly1"     =  Poly1,
						 "poly2"     =  Poly2,
						 "poly3"     =  Poly3,
						 "poly4"     =  Poly4,
						 "poly5"     =  Poly5,
						 "diffPoly12"  =  Poly2 - Poly1,
						 "diffPoly13"  =  Poly3 - Poly1,
						 "diffPoly14"  =  Poly4 - Poly1,
						 "diffPoly15"  =  Poly5 - Poly1
						 )
	
	# Export data
	filename  <-  paste("./output/data/simMultiPatchSgrad", "_C", C, "_delta", delta, "_hf", hf, "_hm", hm, "_sMax", sMax, ".csv", sep="")
	write.csv(data, file=filename, row.names = FALSE)
}



###############################################################
###############################################################
#' Simulation to compare proportion of parameter space where
#' polymorphism is maintained in the 2-patch model compared
#' with the 1-patch expectation 
#'
#' @title 2-Patch SA in hermaphrodites Simulation
#' @param n      Sample size (number of selection coefficients 
#' 				 randomly drawn from uniform distribution)
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function. Assume equal across patches.
#' @param hm     Dominance coefficient for male sex function. Assume equal across patches.
#' @param sMax   Maximum selection coefficient in Patch 1 (determines range
#' 				  of selection coefficient parameter space to be explored)
#' 				  (we assume that we always explore a square parameter space 
#' 				  (i.e., sMax is the same for males and females, and equal 
#' 				  across patches))
#' @export
simMultiPatchGeneralSA  <-  function(n, C, k, delta = 0, hf = 1/2, hm = 1/2, sMax) {

	# Draw random selection coefficients for up to 5 patches
	sfMat  <-  matrix(runif(k*n, max = sMax), ncol=k, nrow=n)
	smMat  <-  matrix(runif(k*n, max = sMax), ncol=k, nrow=n)

	# Perform invasion analysis using general invasion conditions
	generalInv_q0  <-  generalAddInv(q=0, k=k, C=C, sfs=sfMat, sms=smMat)
	generalInv_q1  <-  generalAddInv(q=1, k=k, C=C, sfs=sfMat, sms=smMat)

	# Perform invasion analysis using SA invasion conditions
	SAInv_q0  <-  SAAddInv(q=0, k=k, C=C, sfs=sfMat, sms=smMat)
	SAInv_q1  <-  SAAddInv(q=1, k=k, C=C, sfs=sfMat, sms=smMat)

	# Calculate proportion of parameter space 
	# where polymorphism is predicted to be 
	# maintained by balancing selection
	polyGen  <-  colSums(generalInv_q0 == 1 & generalInv_q1 == 1)/n
	polySA   <-  colSums(SAInv_q0 == 1 & SAInv_q1 == 1)/n

	# Save and return results
	res  <-  list(
				  "polyGen"  =  polyGen,
				  "polySA"   =  polySA
				  )
	return(res)
}



#' Simulation to compare general multi-patch invasion 
#' conditions with those from the SA model across a 
#' gradient of selection coefficients (wrapper function for simMultiPatch())
#'
#' @title Multi-Patch SA in hermaphrodites 
#' @param n      Sample size (number of selection coefficients 
#' 				 randomly drawn from uniform distribution)
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function. Assume equal across patches.
#' @param hm     Dominance coefficient for male sex function. Assume equal across patches.
#' @param sMax   Maximum selection coefficient in Patch 1 (determines range
#' 				  of selection coefficient parameter space to be explored)
#' 				  (we assume that we always explore a square parameter space 
#' 				  (i.e., sMax is the same for males and females, and equal 
#' 				  across patches))
#' @param resolution resolution for selection coefficient gradient
#' @export
simMultiPatchSgradCompareGeneralSAInv  <-  function(n, k, C, delta=0, hf = 1/2, hm = 1/2, 
													sMax = 1, resolution = 0.01) {

	# Create vector of sMaxes (we assume that we always explore a 
	# square parameter space (i.e., sMax is the same for males and females))
	sMaxes  <- seq(from = resolution, to = sMax, by = resolution)
	
	# Initialize storage structures
	PolyGen  <-  matrix(0,ncol=k, nrow=length(sMaxes))
	PolySA   <-  matrix(0,ncol=k, nrow=length(sMaxes))

	# loop over sMaxes
	for (i in 1:length(sMaxes)) {
		res  <-  simMultiPatchGeneralSA(n = n, k=k, C = C, delta = delta, hf = hf, hm = hm, sMax = sMaxes[i]) 
		PolyGen[i,]  <-  res$polyGen
		PolySA[i,]   <-  res$polySA

		# Print Progress
		if(i %% 5 == 0) {
			print(paste("Progress: ", i," / ", length(sMaxes)))
		}
	}

	# make colnames
	PolyGenNames  <-  c()
	PolySANames   <-  c()
	for(i in 1:k) {
		PolyGenNames[i]  <-  paste0("PolyGen",i)
		PolySANames[i]   <-  paste0("PolySA",i)
	}
	colnames(PolyGen)  <-  PolyGenNames	
	colnames(PolySA)   <-  PolySANames	

	# Save results as data frame
	data  <-  data.frame(cbind(sMaxes, PolyGen, PolySA))

	# Export data
	filename  <-  paste("./output/data/simMultiPatchSgrad_General_SA", "_k", k, "_C", C, "_delta", delta, "_hf", hf, "_hm", hm, "_sMax", sMax, ".csv", sep="")
	write.csv(data, file=filename, row.names = FALSE)
}


##################################################################
##################################################################
#  Functions for genotypic recursions & deterministic simulations

#############################################
## Average fitness through each sex function

Wf.av  <-  function(Fij, Wf, ...){
   (Fij[1]*Wf[1]) + (Fij[2]*Wf[2]) + (Fij[3]*Wf[3])
}
Wm.av  <-  function(Fij, Wm, ...){
   (Fij[1]*Wm[1]) + (Fij[2]*Wm[2]) + (Fij[3]*Wm[3])
}
W.tot  <-  function(C, delta, ...){
   1 - (C *delta)
}


#########################################
## Genotypic frequency recursions

FAA.pr  <-  function(Fij, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm, ...) {
	((1 - C)*(((Fij[1]*Wf[1]/Wf.av(Fij,Wf)) * (Fij[1]*Wm[1]/Wm.av(Fij,Wm)))   +
			  ((Fij[1]*Wf[1]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/2 +
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[1]*Wm[1]/Wm.av(Fij,Wm)))/2 +
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/4) +
	 	C*(1 - delta)*(((Fij[1]*Wf[1])/Wf.av(Fij,Wf)) + 
	 		 		   ((Fij[2]*Wf[2])/(4*Wf.av(Fij,Wf))))) / W.tot(C = C, delta = delta)
}
FAa.pr  <-  function(Fij, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm, ...) {
	((1 - C)*(((Fij[1]*Wf[1]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/2 +
			  ((Fij[1]*Wf[1]/Wf.av(Fij,Wf)) * (Fij[3]*Wm[3]/Wm.av(Fij,Wm)))   +
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[1]*Wm[1]/Wm.av(Fij,Wm)))/2 +
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/2 + 
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[3]*Wm[3]/Wm.av(Fij,Wm)))/2 + 
			  ((Fij[3]*Wf[3]/Wf.av(Fij,Wf)) * (Fij[1]*Wm[1]/Wm.av(Fij,Wm)))   +
			  ((Fij[3]*Wf[3]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/2) +
	 	C*(1 - delta)*(((Fij[2]*Wf[2]))/(2*Wf.av(Fij,Wf)))) / W.tot(C = C, delta = delta)
}
Faa.pr  <-  function(Fij, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm, ...) {
	((1 - C)*(((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/4 + 
			  ((Fij[2]*Wf[2]/Wf.av(Fij,Wf)) * (Fij[3]*Wm[3]/Wm.av(Fij,Wm)))/2 + 
			  ((Fij[3]*Wf[3]/Wf.av(Fij,Wf)) * (Fij[2]*Wm[2]/Wm.av(Fij,Wm)))/2 +
			  ((Fij[3]*Wf[3]/Wf.av(Fij,Wf)) * (Fij[3]*Wm[3]/Wm.av(Fij,Wm)))) +
	 	C*(1 - delta)*(((Fij[2]*Wf[2])/(4*Wf.av(Fij,Wf))) + 
	 				   ((Fij[3]*Wf[3]))/Wf.av(Fij,Wf))) / W.tot(C = C, delta = delta)
}


########################
##  Simulation function
########################

#' Forward deterministic simulation of genotypic recursions for
#' 1-locus spatially variable sex-specific selection in 
#' hermaphrodites
#'
#' @title Forward deterministic simulation of genotypic recursions
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen   =  5000,
#'				   C     =  0,
#'				   delta =  0,
#'				   sm    =  0.1,
#'				   sf    =  0.1,
#'				   hm    =  0.5,
#'				   hf    =  0.5
#'				   )
#' @param Fij.init A vector of initial genotypic frequencies (must have length = 3).
#' 				   c(0.99,0.01) for invasion of aabb into population 'fixed' for AABB.
#' 				   c(0.01,0,0.99) for invasion of AABB into population 'fixed' for aabb.
#' @return Returns a list with timeseries for each genotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, Fij.init, threshold = 1e-6) 
recursionFwdSim  <-  function(gen = 5000, C = 0, delta =  0, sm = 0.1, sf = 0.1, hm = 0.5, hf = 0.5, threshold = 1e-7, ...) {

	##  Warnings
	if(any(c(C, delta, sm, sf, hm, hf) < 0) | any(c(C, delta, sm, sf, hm, hf) > 1))
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(hf  !=  hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(hf != 0.5 & hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	##  Fitness vectors
	Wf  <-  c(1       , (1 - hf*sf), (1 - sf))
	Wm  <-  c((1 - sm), (1 - hm*sm),        1)

	##  Initilize data storage structures
	Fij.gen  <-  matrix(0, ncol=3, nrow=gen)

	##  Initial frequencies
	if(hf == hm & hf == 0.5) {
		qHat  <-  qHatAdd(C = C, delta = delta, sf = sf, sm = sm)
		if(qHat <= 0.5)
			Fij.init    <-  c(0.999,0,0.001)
		if(qHat >= 0.5)
			Fij.init    <-  c(0.001,0,0.999)
	}
	if(hf == hm & hf == 0.25) {
		qHat  <-  qHatDomRev(C = C, delta = delta, sf = sf, sm = sm, h = hf)
		if(qHat <= 0.5)
			Fij.init    <-  c(0.999,0,0.001)
		if(qHat >= 0.5)
			Fij.init    <-  c(0.001,0,0.999)
	}

	##  Generation Loop
		# initialize
		Fij.gen[1,1]   <-  as.numeric(rounded(FAA.pr(Fij = Fij.init, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fij.gen[1,2]   <-  as.numeric(rounded(FAa.pr(Fij = Fij.init, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fij.gen[1,3]   <-  as.numeric(rounded(Faa.pr(Fij = Fij.init, Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))


	# Start simulation
	i       <-  2
	diffs   <-  rep(1,3)

	while (i < gen & any(abs(diffs) >= threshold)) {
		Fij.gen[i,1]   <-  as.numeric(rounded(FAA.pr(Fij = Fij.gen[i-1,], Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fij.gen[i,2]   <-  as.numeric(rounded(FAa.pr(Fij = Fij.gen[i-1,], Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fij.gen[i,3]   <-  as.numeric(rounded(Faa.pr(Fij = Fij.gen[i-1,], Wf = Wf, Wm = Wm, C = C, delta = delta, hf = hf, hm = hm), precision=12))
		diffs  <-  Fij.gen[i,] - Fij.gen[i-1,]
		i      <-  i+1
	}

	if (i > gen) {
		MaxGen  <-  1
#		print('Warning: maximum runtime reached. Results may not represent equilibrium frequencies')
	}
	else MaxGen  <-  0

	##  Is equilibrium polymorphic?
	if (any(Fij.gen[i-1,] > 0.9999)) 
		 simPoly  <-  0
	else simPoly  <-  1

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies
	lambda0_1patch  <-  lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = sf, sm = sm)
	lambda1_1patch  <-  lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = sf, sm = sm)

	if (lambda0_1patch > 1 & lambda1_1patch > 1 )
		 eigPoly  <-  1
	else eigPoly  <-  0

	##  Does simulation result agree with Eigenvalues?
	if (simPoly == eigPoly)
		 agree  <-  1
	else agree  <-  0

	##  Output list
	par.list  <-  list(
				   	   gen   =  gen,
					   C     =  C,
					   delta =  delta,
					   sm    =  sm,
					   sf    =  sf,
					   hm    =  hm,
					   hf    =  hf
					   )
	res  <-  list(
				  "par.list" =  par.list,
				  "Fij.gen"  =  Fij.gen[1:i-1,],
				  "EQ.freq"  =  Fij.gen[i-1,],
				  "MaxGen"   =  MaxGen,
				  "l.A"      =  lambda0_1patch,
				  "l.a"      =  lambda1_1patch,
				  "simPoly"  =  simPoly,
				  "eigPoly"  =  eigPoly,
				  "agree"    =  agree
 				 )
	return(res)
}




##  Functions for fitness vectors used in kPatchRecursionFwdSim
Wf.fit  <-  function(hf, sf){
	c(1, (1 - hf*sf), (1 - sf))
}
Wm.fit  <-  function(hm, sm){
	c((1 - sm), (1 - hm*sm), 1)
}


#' Forward deterministic simulation of genotypic recursions for
#' 1-locus spatially variable sex-specific selection in 
#' hermaphrodites
#'
#' @title Forward deterministic simulation of genotypic recursions
#' @param par.list 	A list with desired parameter values for the simulation with structure:
#' @param gen      	Maximum simulation runtime
#' @param k 		Number of patches
#' @param C 		Population selfing rate
#' @param delta 	Population inbreeding depression level
#' @param s.vals	A k x 2 matrix with selection coefficients for each patch
#' @param hm 		Dominance for male-function SA trait
#' @param hf 		Dominance for female-funciton SA trait
#' @param threshold Maximum difference in genotypic frequencies determining when equilibrium is reached
#' @return Returns a list with timeseries for each genotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, Fij.init, threshold = 1e-6)
kPatchRecursionFwdSim  <-  function(gen = 5000, k = 5, C = 0, delta =  0, s.vals, hm = 0.5, hf = 0.5, threshold = 1e-6, ...) {

	##  Warnings
	if(any(c(C, delta, s.vals, hm, hf) < 0) | any(c(C, delta, s.vals, hm, hf) > 1))
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(hf  !=  hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(hf != 0.5 & hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(nrow(s.vals)  !=  k | nrow(s.vals)  !=  k)
		stop('s.vals must have k rows')


	##  Initilize data storage structures
	Fijk.init   <-  array(0, dim=c(1, 3, k))
	Fijk.gen    <-  array(0, dim=c(gen, 3, k))
	Fij.gen     <-  matrix(0, nrow=gen, ncol=3)
	qHats       <-  c(0, length=k)
	Wf.patches  <-  matrix(0,nrow=k, ncol=3)
	Wm.patches  <-  matrix(0,nrow=k, ncol=3)

	##  Calculate initial frequencies for each patch
	for(p in 1:k) {
		if(hf == hm & hf == 0.5) {
			qHat  <-  qHatAdd(C = C, delta = delta, sf = s.vals[p,1], sm = s.vals[p,2])
			if(qHat <= 0.5)
				Fijk.init[,,p]    <-  c(0.999,0,0.001)
			if(qHat >= 0.5)
				Fijk.init[,,p]    <-  c(0.001,0,0.999)
		}
		if(hf == hm & hf == 0.25) {
			qHat  <-  qHatDomRev(C = C, delta = delta, sf = s.vals[p,1], sm = s.vals[p,2], h = hf)
			if(qHat <= 0.5)
				Fijk.init[,,p]    <-  c(0.999,0,0.001)
			if(qHat >= 0.5)
				Fijk.init[,,p]    <-  c(0.001,0,0.999)
		}
	
		# Calculate fitness expressions for each patch
		Wf.patches[p,]  <-  Wf.fit(hf = hf, sf = s.vals[p,1])
		Wm.patches[p,]  <-  Wm.fit(hm = hm, sm = s.vals[p,2])

		# initialize Fijk.gen with  
		Fijk.gen[1,1,p]   <-  as.numeric(rounded(FAA.pr(Fij = Fijk.init[,,p], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fijk.gen[1,2,p]   <-  as.numeric(rounded(FAa.pr(Fij = Fijk.init[,,p], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
		Fijk.gen[1,3,p]   <-  as.numeric(rounded(Faa.pr(Fij = Fijk.init[,,p], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
	}

	# Overall genotypic frequencies after migration in generation 1
	Fij.gen[1,]  <-  apply(Fijk.gen[1,,], MARGIN=1, mean)

	##  Generation Loop
	# Start simulation
	i      <-  2
	diffs  <-  matrix(1, nrow=3, ncol=k)

	while (i <= gen & any(abs(diffs) >= threshold)) {

		for(p in 1:k) {
			Fijk.gen[i,1,p]   <-  as.numeric(rounded(FAA.pr(Fij = Fij.gen[i-1,], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
			Fijk.gen[i,2,p]   <-  as.numeric(rounded(FAa.pr(Fij = Fij.gen[i-1,], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
			Fijk.gen[i,3,p]   <-  as.numeric(rounded(Faa.pr(Fij = Fij.gen[i-1,], Wf = Wf.patches[p,], Wm = Wm.patches[p,], C = C, delta = delta, hf = hf, hm = hm), precision=12))
		}
		Fij.gen[i,]  <-  apply(Fijk.gen[i,,], MARGIN=1, mean)
		diffs  <-  Fij.gen[i,] - Fij.gen[i-1,]
		i      <-  i+1
	}

	##  Is equilibrium polymorphic?
	if (any(Fij.gen[i-1,] > 0.9999)) 
		 simPoly  <-  0
	else simPoly  <-  1

	if (i > gen) {
		MaxGen  <-  1
#		print('Warning: maximum runtime reached. Results may not represent equilibrium frequencies')
	}
	else MaxGen  <-  0

	##  Output list
	par.list  <-  list(
				   	   k     =  k,
				   	   gen   =  gen,
					   C     =  C,
					   delta =  delta,
					   s.vals =  s.vals,
					   hm    =  hm,
					   hf    =  hf
					   )
	res  <-  list(
				  "par.list" =  par.list,
				  "MaxGen"   =  MaxGen,
				  "Fij.gen"  =  Fij.gen[1:i-1,],
				  "EQ.freq"  =  Fij.gen[i-1,],
				  "simPoly"  =  simPoly
 				 )
	return(res)
}

#' Simulation loop wrapping forward deterministic simulations 
#' of genotypic recursions 
#' #'
#' @title Forward deterministic simulation of genotypic recursions.
#' @param n number of randomly generated values for sf & sm. 
#' Determines resolution with which parameter space is explored.
#' @param gen Maximum number of generations for each simulation (as in par.list).
#' @param C The fixed selfing rate (as in par.list).
#' @param hf Dominance through female expression (as in par.list).
#' @param hm Dominance through male expression (as in par.list).
#' @param r.vals Values of recombination rate to explore(as in par.list).
#' @param threshold Threshold difference between genotypic frequencies before simulation cuts off.
#' @return Returns a data frame with parameter values, a variable describing whether 
#' the final state of the simulation was polymorphic polymorphism, whether evaluating the eigenvalues
#' predicts polymorphism, and whether these two methods agree with one another.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop(n = 10000, gen = 5000, C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)
recursionFwdSimLoop  <-  function(n = 10000, gen = 10000, sRange = c(0,1), C = 0, delta = 0, hf = 0.5, hm = 0.5, threshold = 1e-7) {

	## Warnings
	if(any(c(C,hf,hm) < 0) | any(c(C,hf,hm) > 1))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-7)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will effect whether the simulations agree with the analytic results')

	#  initialize selection coeficients and storage structures
	s.vals   <-  matrix(runif(2*n, min=sRange[1], max=sRange[2]), ncol=2)
	simPoly  <-  c()
	eigPoly  <-  c()
	agree    <-  c()


	print('Running Deterministic Recursion Simulations')
	pb   <-  txtProgressBar(min=0, max=nrow(s.vals), style=3)
	setTxtProgressBar(pb, 0)

	##  Simulation Loop over values of r, sm, sf for fixed selfing rate (C)
		for (i in 1:nrow(s.vals)) {
				
			res         <-  recursionFwdSim(gen = gen, C = C, delta =  delta, sf = s.vals[i,1], sm = s.vals[i,2], hm = hm, hf = hf, threshold = threshold)
			simPoly[i]  <-  res$simPoly
			eigPoly[i]  <-  res$eigPoly
			agree[i]    <-  res$agree
		setTxtProgressBar(pb,i)
	}

	#  Compile results as data.frame
	results.df  <-  data.frame("hf"      = rep(hf, length(s.vals)),
							   "hm"      = rep(hm, length(s.vals)),
							   "C"       = rep(C,  length(s.vals)),
							   "delta"   = rep(delta,  length(s.vals)),
							   "sf"      = s.vals[,1],
							   "sm"      = s.vals[,2],
							   "simPoly" = simPoly,
							   "eigPoly" = eigPoly,
							   "agree"   = agree
							   )

	#  Write results.df to .txt file
	filename  <-  paste("./output/data/determFwdSimLoop", "_C", C, "_delta", delta, "_h", hf, "_sMax",sRange[2],"_n", n, ".txt", sep="")
	write.table(results.df, file=filename, col.names = TRUE, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}



#' Multi-patch deterministic simulation
#'
#' @title 2-Patch SA in hermaphrodites Simulation
#' @param n      Sample size (number of selection coefficients 
#' 				 randomly drawn from uniform distribution)
#' @param C      Population selfing rate
#' @param delta  Inbreeding depression
#' @param hf     Dominance coefficient for female sex function. Assume equal across patches.
#' @param hm     Dominance coefficient for male sex function. Assume equal across patches.
#' @param sMax   Maximum selection coefficient in Patch 1 (determines range
#' 				  of selection coefficient parameter space to be explored)
#' 				  (we assume that we always explore a square parameter space 
#' 				  (i.e., sMax is the same for males and females, and equal 
#' 				  across patches))
#' @param resolution resolution for selection coefficient gradient
#' @export
detSimMultiPatchSgrad  <-  function(n = 10000, gen = 5000, C=0, delta=0, hf=0.5, hm=0.5, threshold = 1e-7, sMax = 1, resolution = 0.1) {

	# Create vector of sMaxes (we assume that we always explore a 
	# square parameter space (i.e., sMax is the same for males and females))
	sMaxes  <- c(0.05,seq(from = resolution, to = sMax, by = resolution))
#	sMaxes  <-  c(0.05,0.2,0.4)
	nMaxed  <-  matrix(0, nrow=length(sMaxes), ncol=5)

	# Initialize storage structures
	pSimPoly1  <-  rep(0, length = length(sMaxes))
	pSimPoly2  <-  rep(0, length = length(sMaxes))
	pSimPoly3  <-  rep(0, length = length(sMaxes))
	pSimPoly4  <-  rep(0, length = length(sMaxes))
	pSimPoly5  <-  rep(0, length = length(sMaxes))
	pEigPoly1  <-  rep(0, length = length(sMaxes))
	pEigPoly2  <-  rep(0, length = length(sMaxes))
	pEigPoly3  <-  rep(0, length = length(sMaxes))
	pEigPoly4  <-  rep(0, length = length(sMaxes))
	pEigPoly5  <-  rep(0, length = length(sMaxes))

	# loop over sMaxes
	for (i in 1:length(sMaxes)) {
	print(paste0('Running Deterministic Recursion Simulations for sMax = ',sMaxes[i]))

		# Draw random selection coefficients for up to 5 patches
		s1  <-  matrix(runif(2*n, max = sMaxes[i]), nrow=n, ncol=2)
		s2  <-  matrix(runif(2*n, max = sMaxes[i]), nrow=n, ncol=2)
		s3  <-  matrix(runif(2*n, max = sMaxes[i]), nrow=n, ncol=2)
		s4  <-  matrix(runif(2*n, max = sMaxes[i]), nrow=n, ncol=2)
		s5  <-  matrix(runif(2*n, max = sMaxes[i]), nrow=n, ncol=2)

		simPoly1  <-  rep(0, times=n)
		simPoly2  <-  rep(0, times=n)
		simPoly3  <-  rep(0, times=n)
		simPoly4  <-  rep(0, times=n)
		simPoly5  <-  rep(0, times=n)
		eigPoly1  <-  rep(0, times=n)
		eigPoly2  <-  rep(0, times=n)
		eigPoly3  <-  rep(0, times=n)
		eigPoly4  <-  rep(0, times=n)
		eigPoly5  <-  rep(0, times=n)
		
		# Loop over randomly drawn pairs of selection coefficients for female & male function
		print('sampling sf x sm selection parameter space')
		pb   <-  txtProgressBar(min=0, max=n, style=3)
		setTxtProgressBar(pb, 0)
		maxedOut  <-  matrix(0, nrow=n, ncol=5)
		for(j in 1:n) {

			s.vals.n  <-  rbind(s1[j,],s2[j,],s3[j,],s4[j,],s5[j,])

			# Determine if polymorphic usng deterministic fwd simulations of k-patch recursions 
			res1  <-  recursionFwdSim(gen = gen, C = C, delta =  delta, sf = s.vals.n[1,1], sm = s.vals.n[1,2], hm = hm, hf = hf, threshold = threshold)
			res2  <-  kPatchRecursionFwdSim(gen = gen, k = 2, C = C, delta = delta, s.vals = s.vals.n[1:2,], hm = hm, hf = hf, threshold = threshold)
			res3  <-  kPatchRecursionFwdSim(gen = gen, k = 3, C = C, delta = delta, s.vals = s.vals.n[1:3,], hm = hm, hf = hf, threshold = threshold)
			res4  <-  kPatchRecursionFwdSim(gen = gen, k = 4, C = C, delta = delta, s.vals = s.vals.n[1:4,], hm = hm, hf = hf, threshold = threshold)
			res5  <-  kPatchRecursionFwdSim(gen = gen, k = 5, C = C, delta = delta, s.vals = s.vals.n, hm = hm, hf = hf, threshold = threshold)

			simPoly1[j]  <-  res1$simPoly
			simPoly2[j]  <-  res2$simPoly
			simPoly3[j]  <-  res3$simPoly
			simPoly4[j]  <-  res4$simPoly
			simPoly5[j]  <-  res5$simPoly

			maxedOut[j,1]  <-  res1$MaxGen
			maxedOut[j,2]  <-  res2$MaxGen
			maxedOut[j,3]  <-  res3$MaxGen
			maxedOut[j,4]  <-  res4$MaxGen
			maxedOut[j,5]  <-  res5$MaxGen

			# Determine if polymorphic from eigenvaluse evaluated at boundary equilibria
			# 1-Patch invasion criteria for boundaries of q = 0 and q = 1
			lambda0_1patch  <-  lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2])
			lambda1_1patch  <-  lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2])
	
			# 2-Patch invasion criteria for boundaries of q = 0 and q = 1
			lambda0_2patch  <-  (1/2)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]))
			lambda1_2patch  <-  (1/2)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]))
	
			# 3-Patch invasion criteria for boundaries of q = 0 and q = 1
			lambda0_3patch  <-  (1/3)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]))
			lambda1_3patch  <-  (1/3)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]))
	
			# 4-Patch invasion criteria for boundaries of q = 0 and q = 1
			lambda0_4patch  <-  (1/4)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[4,1], sm = s.vals.n[4,2]))
			lambda1_4patch  <-  (1/4)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[4,1], sm = s.vals.n[4,2]))
	
			# 5-Patch invasion criteria for boundaries of q = 0 and q = 1
			lambda0_5patch  <-  (1/5)*(lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[4,1], sm = s.vals.n[4,2]) + 
									   lambda0(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[5,1], sm = s.vals.n[5,2]))
			lambda1_5patch  <-  (1/5)*(lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[1,1], sm = s.vals.n[1,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[2,1], sm = s.vals.n[2,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[3,1], sm = s.vals.n[3,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[4,1], sm = s.vals.n[4,2]) + 
									   lambda1(C = C, delta = delta, hf = hf, hm = hm, sf = s.vals.n[5,1], sm = s.vals.n[5,2]))

			eigPoly1[j]  <-  sum(1 < lambda0_1patch & 1 < lambda1_1patch)
			eigPoly2[j]  <-  sum(1 < lambda0_2patch & 1 < lambda1_2patch)
			eigPoly3[j]  <-  sum(1 < lambda0_3patch & 1 < lambda1_3patch)
			eigPoly4[j]  <-  sum(1 < lambda0_4patch & 1 < lambda1_4patch)
			eigPoly5[j]  <-  sum(1 < lambda0_5patch & 1 < lambda1_5patch)
		
			setTxtProgressBar(pb, j)
		}

		# record multipatch polymorphism
		pSimPoly1[i]  <-  sum(simPoly1) / n
		pSimPoly2[i]  <-  sum(simPoly2) / n
		pSimPoly3[i]  <-  sum(simPoly3) / n
		pSimPoly4[i]  <-  sum(simPoly4) / n
		pSimPoly5[i]  <-  sum(simPoly5) / n
		pEigPoly1[i]  <-  sum(eigPoly1) / n
		pEigPoly2[i]  <-  sum(eigPoly2) / n
		pEigPoly3[i]  <-  sum(eigPoly3) / n
		pEigPoly4[i]  <-  sum(eigPoly4) / n
		pEigPoly5[i]  <-  sum(eigPoly5) / n

		nMaxed[i,]  <-  colSums(maxedOut)
	}


	# Save results as data frame
	data  <-  data.frame(
						 "sMax"       =  sMaxes,
						 "pSimPoly1"  =  pSimPoly1,
						 "pSimPoly2"  =  pSimPoly2,
						 "pSimPoly3"  =  pSimPoly3,
						 "pSimPoly4"  =  pSimPoly4,
						 "pSimPoly5"  =  pSimPoly5,
						 "pEigPoly1"  =  pEigPoly1,
						 "pEigPoly2"  =  pEigPoly2,
						 "pEigPoly3"  =  pEigPoly3,
						 "pEigPoly4"  =  pEigPoly4,
						 "pEigPoly5"  =  pEigPoly5,
						 "nMaxed1"    =  nMaxed[,1],
						 "nMaxed2"    =  nMaxed[,2],
						 "nMaxed3"    =  nMaxed[,3],
						 "nMaxed4"    =  nMaxed[,4],
						 "nMaxed5"    =  nMaxed[,5]
						 )
	
	# Export data
	filename  <-  paste("./output/data/determSimMultiPatchSgrad_wkSel", "_C", C, "_delta", delta, "_hf", hf, "_hm", hm, "_sMax", sMax, ".csv", sep="")
	write.csv(data, file=filename, row.names = FALSE)
}