#############################################################
#  1-locus SA with partial selfing and inbreeding depression
#
#  R code to run deterministic forward simulations
#  of genotypic frequency recursions
#
#  Appendix XXX:
#  
#  Article title goes here...
#  
#  Authors: Colin Olito, Jessica K. Abbott, Crispin Y. Jordan
#
#  NOTES:  
#          
#	  Deterministic simulations explore essentially the same 
#	  parameter conditions as presented in Fig.1 of the main
#	  text. The results demonstrate that the weak-selection 
#	  approximations used to derive the eigenvalues work well
#	  when compared against the deterministic genotypic recursions,
#	  which make no such assumption.
#	  
#		-- hf = hm = (1/2, 1/4)
#		-- C  =  (0, 0.25, 0.5, 0.75)
#		-- delta  =  (0, 0.5)



rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-figures.R')



######################################
#  Additive effects (hf = hm = 0.5)

##############
# 1: C = 0
	# delta = 0
	Add1  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0, delta = 0, hf = 0.5, hm = 0.5, threshold = 1e-10)
	# delta = 1/2
	Add2  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0, delta = 0.5, hf = 0.5, hm = 0.5, threshold = 1e-10)

##############
# 2: C = 0.25
	# delta = 0
	Add3  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.25, delta = 0, hf = 0.5, hm = 0.5, threshold = 1e-10)
	# delta = 1/2
	Add4  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.25, delta = 0.5, hf = 0.5, hm = 0.5, threshold = 1e-10)

##############
# 3: C = 0.5
	# delta = 0
	Add5  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.5, delta = 0, hf = 0.5, hm = 0.5, threshold = 1e-10)
	# delta = 1/2
	Add6  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.5, delta = 0.5, hf = 0.5, hm = 0.5, threshold = 1e-10)

##############
# 4: C = 0.75
	# delta = 0
	Add7  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.75, delta = 0, hf = 0.5, hm = 0.5, threshold = 1e-10)
	# delta = 1/2
	Add8  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.75, delta = 0.5, hf = 0.5, hm = 0.5, threshold = 1e-10)



######################################
#  Dominance Reversals (hf = hm < 0.5)

##############
# 1: C = 0
	# delta = 0
	DR1  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0, delta = 0, hf = 0.25, hm = 0.25, threshold = 1e-10)
	# delta = 1/2
	DR2  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0, delta = 0.5, hf = 0.25, hm = 0.25, threshold = 1e-10)

##############
# 2: C = 0.25
	# delta = 0
	DR3  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.25, delta = 0, hf = 0.25, hm = 0.25, threshold = 1e-10)
	# delta = 1/2
	DR4  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.25, delta = 0.5, hf = 0.25, hm = 0.25, threshold = 1e-10)

##############
# 3: C = 0.5
	# delta = 0
	DR5  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.5, delta = 0, hf = 0.25, hm = 0.25, threshold = 1e-10)
	# delta = 1/2
	DR6  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.5, delta = 0.5, hf = 0.25, hm = 0.25, threshold = 1e-10)

##############
# 4: C = 0.75
	# delta = 0
	DR7  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.75, delta = 0, hf = 0.25, hm = 0.25, threshold = 1e-10)
	# delta = 1/2
	DR8  <-  recursionFwdSimLoop(n = 10000, gen = 10000, sRange = c(0,1), C = 0.75, delta = 0.5, hf = 0.25, hm = 0.25, threshold = 1e-10)
