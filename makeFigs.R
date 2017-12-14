#  Functions to generate figures for: 
#    
#  Title: 	The interaction between sexual antagonism 
#			and local adaptation in species without
#			separate sexes
#
#			A Contributed paper for Phil. Trans.
#			Roy. Soc. Theme Issue put together by 
#		   	the ESEB Special Topics Network: linking 
#			local adaptation with the evolution of 
#			sex-differences.
#
#
#  Author: Colin Olito, Jessica K. Abbott, and Crispin Y. Jordan
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())
###############
# Dependencies
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

source('R/functions-analyses.R')
source('R/functions-figures.R')


########################
# Figures for the paper
########################

# toPdf(Fig.1(), figPath(name='Fig1.pdf'), width=5, height=7.75)
# embed_fonts(figPath(name='Fig1.pdf'))

# toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig2.pdf'))

# toPdf(Fig3Alt(), 
#             figPath(name='Fig3Alt.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig3Alt.pdf'))



########################
# Supplementary Figures
########################






######################
# Exploratory Figures
######################

# Compare SA polymorphism between 1 and 2 patch models
toPdf(proportionPoly2Patch(h = 1/2, delta = 0), figPath(name='propPoly2Patch_Additive_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPoly2Patch_Additive_delta0.pdf'))

toPdf(proportionPoly2Patch(h = 1/4, delta = 0), figPath(name='propPoly2Patch_DomRev_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPoly2Patch_DomRev_delta0.pdf'))


# Difference in SA polymorphism between 1 and 2 patch models
toPdf(diffPoly2Patch(h = 1/2, delta = 0), figPath(name='diffPoly2Patch_Additive_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='diffPoly2Patch_Additive_delta0.pdf'))

toPdf(diffPoly2Patch(h = 1/4, delta = 0), figPath(name='diffPoly2Patch_DomRev_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='diffPoly2Patch_DomRev_delta0.pdf'))