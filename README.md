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

First colum = Population identifier; second colum = Body size measurement such as snout-vent-length (SVL); third colum onwards = other morphometric characters.

Example:

Population | SVL | HW | HL | SNL 
--- | --- | --- | --- | ---
A | 23.5 | 8.8 | 12.5 | 5.6
A | 24.0 | 8.5 | 12.0 | 5.2
B | 22.9 | 7.2 | 10.3 | 4.8
B | 22.4 | 7.7 | 10.2 | 4.6

Singletons, juvenile measurements, and missing data should not be included


## Allometric adjustment of morphological characters
The `body_adjust()` function uses an allometric growth equation to adjust for ontogenetic variation among **populations** within a species, i.e. **intraspecific variation**. Hence, this function should be applied separately for different species. The function uses the following allometric equation: *X*<sub>adj</sub> = log(*X*)-*B*[log(SVL)-log(SVL<sub>mean</sub>)], where *X*<sub>adj</sub> = Adjusted value for character *X*; *X* = raw/unadjusted value for character *X*; *B* = unstandardized regression coefficient for *X* against SVL; SVL = measured SVL; SVL<sub>mean</sub> = mean SVL for all populations.

Note that each population will have its own *B* and SVL<sub>mean</sub> is the average SVL of all **populations** included in the dataset.

This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 

Example:
`m <- read.csv("foo.csv")
body_adjust(m)`

## PCA
The `GS_pca()` function is a wrapper that uses native `prcomp()` function with scaling and outputs `ggplot` graphs and a summary table

# CITATION
If you use this package, please cite it as:

Kin Onn Chan (2021). GroupStruct (R package). https://github.com/chankinonn/GroupStruct

# References
Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to intraspecific variation in the Ringed Snake *Natrix natrix* (L.). *Biological Journal of the Linnaean Society*, 7: 27-43

Thorpe, R. S. (1976) Biometric analysis of geographic variation and racial affinities. *Biological Reviews*, 51: 407-452

Thorpe, R. S. (1983) A review of the numerical methods for recognising and analysing racial differentiation, in: Felsenstein, J. (Ed.), *Numerical Taxonomy*. Springer Verlag, Berlin, pp. 404-423

Lleonart, J., Salat, J., and Torres, G. J. (2000) Removing Allometric Effects of Body Size in Morphological Analysis. *Journal of Theoretical Biology*, 205: 85-93

