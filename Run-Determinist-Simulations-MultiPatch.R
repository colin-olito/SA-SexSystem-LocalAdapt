#########################################################################
#	Multi-patch spatially variable selection for simultaneous hermaphrodites
#	Run simulations to estimate proportion of parameter space shwere SA
#	polymorphism is maintained by balancing selection across a gradient
#	of selection coefficients ranging from 0 to sMax
#
#	Generates output data as .csv files saved to ./output/data/simResults.
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
source('R/functions-figures.R')


# Set constant parameters for simulations
n           <-  10000
sMax        <-  1
resolution  <-  0.2

######################
##  Run Simulations

# Additive Fitness effects
	h      <-  1/2
	# No Inbreeding depression (delta = 0)
		delta  <-  0
		
		# C = 0
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 0,   delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/2
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/2, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 3/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 3/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)

	# Moderate inbreeding depression (delta = 0.5)
		delta  <-  1/2
		
		# C = 0
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 0,   delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/2
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/2, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 3/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 3/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)

# Dominance Reversal
	h      <-  1/4
	# No Inbreeding depression (delta = 0)
		delta  <-  0
		
		# C = 0
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 0,   delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/2
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/2, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 3/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 3/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)

	# Moderate inbreeding depression (delta = 0.5)
		delta  <-  0.5
		
		# C = 0
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 0,   delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 1/2
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 1/2, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
		# C = 3/4
		detSimMultiPatchSgrad(n = n, gen = 200000, C = 3/4, delta = delta, hf = h, hm = h, sMax = sMax, resolution = resolution, threshold = 1e-10)
