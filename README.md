# DESCRIPTION
GroupStruct contains tools for analyzing and visualizing morphological data. Functions include body-size corrections of morphological characters to correct for ontogenetic variation in body size and various clustering/dimension reduction analyses to evince group structure.  
Please post questions and bug reports to groupstruct@googlegroups.com

# INSTALLATION
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chankinonn/GroupStruct")
```

# INPUT
Input data should be in csv format and configured as follow:

1st colum = Species/OTU identifier\
2nd colum = Standard body size measurement (e.g. snout-vent-length, SVL)\
3rd colum onwards = other morphometric characters.

Example:

OTU | SVL | HW | HL | SNL 
--- | --- | --- | --- | ---
A | 23.5 | 8.8 | 12.5 | 5.6
A | 24.0 | 8.5 | 12.0 | 5.2
A | 24.3 | 8.2 | 12.2 | 5.5
B | 22.9 | 7.2 | 10.3 | 4.8
B | 22.4 | 7.7 | 10.2 | 4.6
B | 22.4 | 7.6 | 10.4 | 5.0

## Allometric adjustment of morphological characters for multipecies and multipopulatin datasets
The `allom()` function uses the following allometric growth equation to adjust for ontogenetic variation (Thorpe, 1975):

<div align="center">X<sub>adj</sub> = log(X)-b[log(BL)-log(BL<sub>mean</sub>)]</div> 

\
where *X*<sub>adj</sub> = Adjusted value for character *X*; *X* = raw/unadjusted value for character *X*; *b* = regression coefficient (slope) for log(*X*) against log(BL); BL = measurement of body length/size; BL<sub>mean</sub> = grand mean of BL. This method removes all the information related to size, scales all individuals to the same size, and adjusts their shape to the new size according to allometry. 

A multispecies dataset comprises more than one species where each unique identifier in the first column represents a different species. If there are several populations/localities within a species, these should be grouped under a common identifier (pooled groups *sensu* Thorpe (1975)). A multipopulation dataset comprises a single species where each unique identifier represents a different population/locality (common within-group pooling). The only difference in calculations between multispecies and multipopulation datasets is in BLmean. For multispecies datasets, a separate BLmean will be calculated for each species, whereas in a multipopulation dataset, a single BLmean is calculated across all populations to represent the average size of the species.

Each OTU should be represented by more than two individuals for the adjustment to work and missing data must not be included. This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 

Example:\
Multispecies dataset:\
`m <- read.csv("foo.csv")`\
`allom(m, type = "species")`

Multipopulation dataset:\
`m <- read.csv("foo.csv")`\
`allom(m, type = "population")`

## Residuals
The function `resid()` calculates residuals of a linear regression of each character against body size

## PCA
The `ez_pca()` function peforms PCA using `prcomp()` with scaling and outputs `ggplot` graphs and a summary table. No transformations are done on the data so all transformations/adjustments (e.g. log transformation) should be done prior to using this function.

# CITATION
If you use this package, please cite it as:

Chan, K.O. & Grismer, L. L. (2021). Correcting for Body Size Variation in Morphometric Analysis. *bioRxiv*, 2021.05.17.444580; doi: https://doi.org/10.1101/2021.05.17.444580 

# References
Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to intraspecific variation in the Ringed Snake *Natrix natrix* (L.). *Biological Journal of the Linnaean Society*, 7: 27-43

Thorpe, R. S. (1976) Biometric analysis of geographic variation and racial affinities. *Biological Reviews*, 51: 407-452

Thorpe, R. S. (1983) A review of the numerical methods for recognising and analysing racial differentiation, in: Felsenstein, J. (Ed.), *Numerical Taxonomy*. Springer Verlag, Berlin, pp. 404-423

Reist, J. D. (1985) An empirical evaluation of several univariate methods that adjust for size variation in morphometric data. *Canadian Journal of Zoology*, 63: 1429-1439

Lleonart, J., Salat, J., and Torres, G. J. (2000) Removing Allometric Effects of Body Size in Morphological Analysis. *Journal of Theoretical Biology*, 205: 85-93

