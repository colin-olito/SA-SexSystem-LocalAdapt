#########################################################################
#	2-patch spatially variable selection for simultaneous hermaphrodites
#	Run simulations to estimate proportion of parameter space shwere SA
#	polymorphism is maintained by balancing selection across a gradient
#	of selection coefficients ranging from 0 to sMax
#
#  R code for W-F forward simulations. Generates output data
#  as .csv files saved to ./output/data/simResults.
#
#
#  Author: Colin Olito
#
#  NOTES:  
#          

rm(list=ls())
#####################
##  Dependencies

source('R/functions-analyses.R')


# Set constant parameters for simulations
n           <-  100000
sMax        <-  1
resolution  <-  0.01

######################
##  Run Simulations

# Additive Fitness effects
	h      <-  1/2
	# No Inbreeding depression (delta = 0)
		delta  <-  0
		
		# C = 0
		simMultiPatchSgrad(n = n, C = 0,   delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/4
		simMultiPatchSgrad(n = n, C = 1/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/2
		simMultiPatchSgrad(n = n, C = 1/2, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 3/4
		simMultiPatchSgrad(n = n, C = 3/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)

	# Moderate inbreeding depression (delta = 0.4)
		delta  <-  0.4
		
		# C = 0
		simMultiPatchSgrad(n = n, C = 0,   delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/4
		simMultiPatchSgrad(n = n, C = 1/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/2
		simMultiPatchSgrad(n = n, C = 1/2, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 3/4
		simMultiPatchSgrad(n = n, C = 3/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)

# Dominance Reversal
	h      <-  1/4
	# No Inbreeding depression (delta = 0)
		delta  <-  0
		
		# C = 0
		simMultiPatchSgrad(n = n, C = 0,   delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/4
		simMultiPatchSgrad(n = n, C = 1/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 1/2
		simMultiPatchSgrad(n = n, C = 1/2, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)
		# C = 3/4
		simMultiPatchSgrad(n = n, C = 3/4, delta = delta, hf = h, hm = h, sMax=1, resolution = 0.01)

