# DESCRIPTION
GroupStruct contains tools for analyzing and visualizing morphological and DNA sequence data for the purpose of detecting group structure. Functions include body-size corrections of morphological characters to correct for ontogenetic variation in body size and various clustering/dimension reduction analyses to evince group structure.  
Please post questions and bug reports to groupstruct@googlegroups.com

# INSTALLATION
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chankinonn/GroupStruct")
```

# FUNCTIONS
Input data should be in csv format and arranged as follow:

1st colum = Species/OTU identifier\
2nd colum = Standard body size measurement (e.g. snout-vent-length, SVL)\
3rd colum onwards = other morphometric characters.

Example:

Population | SVL | HW | HL | SNL 
--- | --- | --- | --- | ---
A | 23.5 | 8.8 | 12.5 | 5.6
A | 24.0 | 8.5 | 12.0 | 5.2
A | 24.3 | 8.2 | 12.2 | 5.5
B | 22.9 | 7.2 | 10.3 | 4.8
B | 22.4 | 7.7 | 10.2 | 4.6
B | 22.4 | 7.6 | 10.4 | 5.0

Missing data should not be included


## Allometric adjustment of morphological characters for datasets containing multiple species (`allom()`) and multiple populations of a single species (`allom_pop()`)
The `allom()` function uses the following allometric growth equation to adjust for ontogenetic variation (Thorpe, 1975):\
*X*<sub>adj</sub> = log(*X*)-*b*[log(SVL)-log(SVL<sub>mean</sub>)]\
where *X*<sub>adj</sub> = Adjusted value for character *X*; *X* = raw/unadjusted value for character *X*; *b* = regression coefficient (slope) for log(*X*) against log(SVL); SVL = measured SVL; SVL<sub>mean</sub> = grand mean SVL for the species. This method removes all the information related to size, scales all individuals to the same size, and adjusts their shape to the new size according to allometry. 

This function assumes the each unique identifier in the first colum of the dataset represents a distinct species. If there are several populations/localities within a species, these should be grouped under a common name (pooled groups *sensu* Thorpe (1975)). Each species should be represented by more than two individuals for the adjustment to work and missing data must not be included. This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 

To adjust for intraspecific variation (e.g. geographic variation among populations of the same species), use the `allom_pop()` function. In this case, the input dataset consists of a single species with the first column representing population names. The slope (*b*) is thus calculated using common within-groups pooling and SVL<sub>mean</sub> is the grand mean calculated by averaging over all populations. 

Example:\
`m <- read.csv("foo.csv")`\
`allom(m)`

## Residuals
The function `resid()` calculates residuals of a linear regression of each character against body size

## PCA
The `GS_pca()` function is a wrapper that uses native `prcomp()` function with scaling and outputs `ggplot` graphs and a summary table

# CITATION
If you use this package, please cite it as:

Kin Onn Chan (2021). GroupStruct (R package). https://github.com/chankinonn/GroupStruct

# References
Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to intraspecific variation in the Ringed Snake *Natrix natrix* (L.). *Biological Journal of the Linnaean Society*, 7: 27-43

Thorpe, R. S. (1976) Biometric analysis of geographic variation and racial affinities. *Biological Reviews*, 51: 407-452

Thorpe, R. S. (1983) A review of the numerical methods for recognising and analysing racial differentiation, in: Felsenstein, J. (Ed.), *Numerical Taxonomy*. Springer Verlag, Berlin, pp. 404-423

Reist, J. D. (1985) An empirical evaluation of several univariate methods that adjust for size variation in morphometric data. *Canadian Journal of Zoology*, 63: 1429-1439

Lleonart, J., Salat, J., and Torres, G. J. (2000) Removing Allometric Effects of Body Size in Morphological Analysis. *Journal of Theoretical Biology*, 205: 85-93

