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
#lambda0  <-  function(C=0, delta=0, hf=1/2, hm=1/2, sf=0.01, sm=0.01) {
#	(4 + 2*hf*sf*(sm - 1) - 2*(1 + hm)*sm + C*(-2 - sf - sm + 4*hm*sm + sf*sm - 
#    4*(-1 + hf*sf)*(-1 + sm)*delta) + (C^2)*(sm - 2*hm*sm + 2*delta - 
#    2*sm*delta + (-1 + 2*hf)*sf*(-1 + sm)*(-1 + 2*delta))) / (2*(-2 + C)*(-1 + sm))
#}
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
#lambda1  <-  function(C=0, delta=0, hf=1/2, hm=1/2, sf=0.01, sm=0.01) {
#	((-1 + C)*(-1 + sf)*(2 - 2*hm*sm + C*(-1 + (-1 + 2*hm)*sm)) - (2 - 
#    C + 2*(-1 + C)*hf*sf)*(-1 + C*(-1 + 2*delta))) / (2*(-2 + C)*(-1 + sf))
#}
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
#	onePatchLB  <-  Inv0SinglePatch(C = C, delta = delta, hf = hf, hm = hm, sm = sm1)
#	onePatchUB  <-  Inv1SinglePatch(C = C, delta = delta, hf = hf, hm = hm, sm = sm1)
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

	# 4-Patch invasion criteria for boundaries of q = 0 and q = 1
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




#' Simulation to icompare proportion of parameter space where
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




