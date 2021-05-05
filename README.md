# DESCRIPTION
GroupStruct contains tools for analyzing and visualizing morphological and DNA sequence data for the purpose of detecting group structure. Functions include allometric adjustment of morphological characters to correct for ontogenetic variation in body size and various clustering/dimension reduction analyses to evince group structure.  
Please post questions and bug reports to groupstruct@googlegroups.com

# INSTALLATION
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chankinonn/GroupStruct")
```

# FUNCTIONS
Input data should be in csv format and arranged as follow:

1st colum = Population identifier\
2nd colum = Body size measurement (e.g. snout-vent-length, SVL)\
3rd colum onwards = other morphometric characters.

Example:

Population | SVL | HW | HL | SNL 
--- | --- | --- | --- | ---
A | 23.5 | 8.8 | 12.5 | 5.6
A | 24.0 | 8.5 | 12.0 | 5.2
B | 22.9 | 7.2 | 10.3 | 4.8
B | 22.4 | 7.7 | 10.2 | 4.6

Singletons, juvenile measurements, and missing data should not be included


## Allometric adjustment of morphological characters
The `body_adjust()` function uses an allometric growth equation to adjust for ontogenetic variation among **populations** within a species, i.e. **intraspecific variation**. Hence, this function should be applied separately for different species. The function uses the following allometric equation: *X*<sub>adj</sub> = log(*X*)-*b*[log(SVL)-log(SVL<sub>mean</sub>)], where *X*<sub>adj</sub> = Adjusted value for character *X*; *X* = raw/unadjusted value for character *X*; *b* = unstandardized regression coefficient for log(*X*) against log(SVL); SVL = measured SVL; SVL<sub>mean</sub> = mean SVL for all individuals.

Note that for *b*, three types of coefficients can be estimated according to group structure: pooled groups combines individuals from different localities to produce a single compound locality; common within-groups pool individuals from the same locality (each locality represents a separate pool); and individual within-groups considers each individual separately. Pooled groups should only be used if the slopes of within group pools have been shown to be equal. This function implements common within-groups pooling to account for potential geographically-based variation. The grouping category can be adjusted according to other types of group-based variation. This equation removes all information related to size, scales individuals to the same size, and adjusts their shape to the new size.  

This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 

Example:\
`m <- read.csv("foo.csv")`\
`body_adjust(m)`

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

