# The interaction between sexual antagonism and local adaptation in species without separate sexes


## Overview

This is a GitHub repository for the development of a collaborative research project begun at the 2017 ESEB Special Topics Network: *linking local adaptation with the evolution of sex-differences*. A major goal of the workshop was to put together a proposal for a Theme Issue of Phil Trans. Roy. Soc. B on sex-specific local adaptation. This project is intended to result in a contributed paper to this project. Here you can find all of the necessary code to reproduce the analyses and the pre-publication manuscript. The *Mathematica* `.nb` files where the models are developed are not tracked by git, but can be downloaded as an online appendix from the publisher website above (or contact me and I will be happy to provide if you don't have access to the appendices). Aside from the `.nb` files, all necessary code for creating the figures can be found in the `./R/functions-*.R` files. 


## Abstract for the proposed Phil. Trans. Theme Issue

Due to basic reproductive differences, males and females rarely share the same optimal phenotype. The resulting patterns of sex-specific selection may also vary considerably among environments (e.g., due to environmental heterogeneity). Both local and sex-specific selection are thought to be common, but research that addresses their interaction has traditionally focused on species with separate sexes. Here, we consider this issue in other sexual systems, focusing on hermaphrodites. Local adaptation appears widespread among hermaphrodites, despite some theoretical predictions that self-fertilization may erode local adaptation. However, theory also suggests that opposing selection through male and female sex-functions can maintain sexually antagonistic (SA) genetic variation for fitness in hermaphrodites. Here, we review the literature to explore similarities and differences in the consequences of SA selection on local adaptation for hermaphrodites versus separate-sexed species. We then develop a simple theoretical model to illustrate how self-fertilization interacts with spatial variation in sex-specific selection to determine local adaptation. Our results provide insight into counterintuitive empirical patterns of local adaptation in selfing vs. outcrossing species, and into our expectations for sex-specific local adaptation over a diverse spectrum of sexual systems.


## Citing information

Citing information for the resulting paper will be provided when it is made [available through the publisher](http://XXXXX). You can also contact me directly if you would like a reprint. 


## Reproducing the manuscript

The easiest way to reproduce the manuscript is to simply clone the repo, run `createFigs.R`, and then compile the manuscript file `doc/SA-SexSystem-LocalAdapt.tex` using whatever default LaTeX editor/engine you have. 


## Contact & bug reporting

Please report any bugs, problems, or issues by opening an issue on the XvAutosomeInversions github [issues page](https://github.com/colin-olito/SA-SexSystem-LocalAdapt/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito at monash dot edu.