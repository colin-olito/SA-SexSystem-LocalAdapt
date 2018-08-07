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

toPdf(proportionPolyMultiPatch(h = 1/2, delta = 0), figPath(name='Fig1.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig1.pdf'))

toPdf(proportionPolyMultiPatchHDelta(h = c(0.5,0.25), delta=c(0,0)), figPath(name='Fig1Alt.pdf'), width=10, height=5)
embed_fonts(figPath(name='Fig1Alt.pdf'))

########################
# Supplementary Figures
########################

toPdf(proportionPolyMultiPatchHDelta(h = c(0.5,0.25), delta=c(0.5,0.5)), figPath(name='FigA1-alt.pdf'), width=10, height=5)
embed_fonts(figPath(name='FigA1-alt.pdf'))

toPdf(invConditionsSA(), figPath(name='deltaEffect.pdf'), width=10, height=5)
embed_fonts(figPath(name='deltaEffect.pdf'))

toPdf(compareSimEig1PatchFig(h = c(0.5,0.25), delta=0), figPath(name='compareEigSim1Patch.pdf'), width=10, height=5)
embed_fonts(figPath(name='compareEigSim1Patch.pdf'))

toPdf(compareSimEig1PatchFig(h = c(0.5,0.25), delta=0.5), figPath(name='compareEigSim1PatchDelta.pdf'), width=10, height=5)
embed_fonts(figPath(name='compareEigSim1PatchDelta.pdf'))

toPdf(compareEigLinesSimPolyMultiPatchHDelta(h = c(0.5,0.25), delta=c(0,0)), figPath(name='compareEigLinesSimMultiPatch.pdf'), width=10, height=5)
embed_fonts(figPath(name='compareEigLinesSimMultiPatch.pdf'))

toPdf(compareEigLinesSimPolyMultiPatchHDelta(h = c(0.5,0.25), delta=c(0.5,0.5)), figPath(name='compareEigLinesSimMultiPatchDelta.pdf'), width=10, height=5)
embed_fonts(figPath(name='compareEigLinesSimMultiPatchDelta.pdf'))

toPdf(PolyMultiPatchCompareGeneralSA(k = 5,sMax=1), figPath(name='compareGeneralSA-Polymorphism.pdf'), width=7, height=7)
embed_fonts(figPath(name='compareGeneralSA-Polymorphism.pdf'))

##################
# Joint Dist Figs
##################
toPdf(PolyMultiPatchJointDist(k=5, h=0.5, delta=0, sMax=1), figPath(name='Polymorphism-JointDist-h0.5-delta0.pdf'), width=9, height=18)
embed_fonts(figPath(name='Polymorphism-JointDist-h0.5-delta0.pdf'))

toPdf(PolyMultiPatchJointDist(k=5, h=0.25, delta=0, sMax=1), figPath(name='Polymorphism-JointDist-h0.25-delta0.pdf'), width=9, height=18)
embed_fonts(figPath(name='Polymorphism-JointDist-h0.25-delta0.pdf'))

toPdf(PolyMultiPatchJointDist(k=5, h=0.5, delta=0.5, sMax=1), figPath(name='Polymorphism-JointDist-h0.5-delta0.5.pdf'), width=9, height=18)
embed_fonts(figPath(name='Polymorphism-JointDist-h0.5-delta0.5.pdf'))

toPdf(PolyMultiPatchJointDist(k=5, h=0.25, delta=0.5, sMax=1), figPath(name='Polymorphism-JointDist-h0.25-delta0.5.pdf'), width=9, height=18)
embed_fonts(figPath(name='Polymorphism-JointDist-h0.25-delta0.5.pdf'))

######################
# Exploratory Figures
######################
source('R/functions-figures.R')
# Compare SA polymorphism between 1 and 2 patch models
toPdf(proportionPolyMultiPatch(h = 1/2, delta = 0), figPath(name='propPolyMultiPatch_Additive_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPolyMultiPatch_Additive_delta0.pdf'))

toPdf(proportionPolyMultiPatch(h = 1/2, delta = 0.5), figPath(name='propPolyMultiPatch_Additive_delta0.5.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPolyMultiPatch_Additive_delta0.5.pdf'))

toPdf(proportionPolyMultiPatch(h = 1/4, delta = 0), figPath(name='propPolyMultiPatch_DomRev_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPolyMultiPatch_DomRev_delta0.pdf'))

toPdf(proportionPolyMultiPatch(h = 1/4, delta = 1/2), figPath(name='propPolyMultiPatch_DomRev_delta0.5.pdf'), width=7, height=7)
embed_fonts(figPath(name='propPolyMultiPatch_DomRev_delta0.5.pdf'))


# Difference in SA polymorphism between 1 and 2 patch models
toPdf(diffPolyMultiPatch(h = 1/2, delta = 0), figPath(name='diffPolyMultiPatch_Additive_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='diffPolyMultiPatch_Additive_delta0.pdf'))

toPdf(diffPolyMultiPatch(h = 1/4, delta = 0), figPath(name='diffPolyMultiPatch_DomRev_delta0.pdf'), width=7, height=7)
embed_fonts(figPath(name='diffPolyMultiPatch_DomRev_delta0.pdf'))
