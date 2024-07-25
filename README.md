# DESCRIPTION
GroupStruct contains tools for analyzing and visualizing morphological and genetic data to evince group structure. Functions include body-size corrections of morphological characters to correct for ontogenetic variation in body size, various clustering/dimension reduction analyses, and comparing genetic distances (p-distance) of species based on a user-specific threshold  
Please post questions and bug reports to groupstruct@googlegroups.com

# INSTALLATION
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("chankinonn/GroupStruct")
```

# Allometric adjustment of morphological characters for multispecies and multipopulation datasets
The `allom()` function uses the following allometric growth equation to adjust for ontogenetic variation (Thorpe, 1975):

<div align="center">X<sub>adj</sub> = log(X)-b[log(BL)-log(BL<sub>mean</sub>)]</div> 

\
where *X*<sub>adj</sub> = Adjusted value for character *X*; *X* = raw/unadjusted value for character *X*; *b* = regression coefficient (slope) for log(*X*) against log(BL); BL = measurement of body length/size; BL<sub>mean</sub> = grand mean of BL. This method removes all the information related to size, scales all individuals to the same size, and adjusts their shape to the new size according to allometry. 

A multispecies dataset comprises more than one species where each unique identifier in the first column represents a different species. If there are several populations/localities within a species, these should be grouped under a common identifier (pooled groups *sensu* Thorpe (1975)). A multipopulation dataset comprises a single species where each unique identifier represents a different population/locality. 

There are two ways to calculate the slope depending on how groups are structured: (i) Pooled groups combine different localities/populations to form a single compound locality to boost sample size. However, this should only be done if populations are not significantly heterogeneous (e.g., no geographic variation). (ii) Common within-group pooling considers each population separately and a different slope is calculated for each population. For type = population1 in a multipopulation dataset, a separate slope is calculated for each population (common within-group pooling). If type = population2, all populations are pooled (pooled group) and a single slope is calculated based on combined populations. For multispecies datasets (type = species), populations are pooled under a single species identifier and slopes are calculated based on pooled groups. If population-level variation is present within the multispecies dataset, each species should be analyzed as separate multipopulation datasets using type = population1.


The other difference in calculation between multispecies and multipopulation datasets is in BLmean. For multispecies datasets, a separate BLmean will be calculated for each species, whereas in a multipopulation dataset, a single BLmean is calculated across all populations to represent the average size of the species.

Each OTU should be represented by more than two individuals for the adjustment to work and missing data must not be included. This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 

## INPUT
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

Example:\
Multispecies datasets:\
`m <- read.csv("foo.csv")`\
`allom(m, type = "species")`

Multipopulation datasets (common within-group):\
`m <- read.csv("foo.csv")`\
`allom(m, type = "population1")`

Multipopulation datasets (pooled group):\
`m <- read.csv("foo.csv")`\
`allom(m, type = "population2")`

# PCA
The `ez_pca()` function peforms PCA using `prcomp()` with scaling and outputs `ggplot` graphs and a summary table. No transformations are done on the data so all transformations/adjustments (e.g. log transformation) should be done prior to using this function.

# P-distance comparisons
The `pdist_cutoff()` function calculates p-distances using a FASTA-formatted alignment and identifies pairs of species that exceeds a user-specified threshold. For example, this function can be used to identify combinations of species pairs that exceed a p-distance threshold of 3% (users can specify whatever threshold they want). For this function to work, multiple sequences need to be aligned as a single FASTA-formatted alignment file and labeled as genus_species. Other unique identifiers can be appended after the species name using additoinal underscores. E.g. Bufo_bufo_MN4578. The output table will be written to the working directory and a list of species that meet the cutoff will be printed.

Example:\
Set path to fasta alignment\
`myfasta <- "path/to/fasta/fasta_alignment.fas"`

Set cutoff value\
`mycutoff <- 0.05`

Set your output file name\
`raw_output <- "raw_output.csv"`
`aggregated_output <- "aggregated_output.csv"`

Get p-distances\
`results <- pdist_cutoff(myfasta, mycutoff, raw_output, aggregated_output)`

# Plotting P-dsitance comparisons
The `plot_pdist()` function takes the results of the `pdist_cutoff()` and plots boxplots of a chosen focal species compared against a list of selected species. E.g., if you wish you compare the p-distance of Species A (focal_species) with Species B, C, and D (comparison_species):\

Specify the focal species  and the list of comparison species\
`focal_species <- "Micryletta_subaraji"`
`comparison_species <- c("Micryletta_sumatrana", "Micryletta_dissimulans", "Micryletta_inornata", "Micryletta_nigromaculata")`

Boxplots for the specified species combinations\
`plot_pdist(results, focal_species, comparison_species)`


# CITATION
If you use this package, please cite it as:

Chan, K.O. & Grismer, L. L. (2021). A standardized and statistically defensible framework for quantitative morphological analyses in taxonomic studies. *Zootaxa*, 5023: 293-300 

# References
Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to intraspecific variation in the Ringed Snake *Natrix natrix* (L.). *Biological Journal of the Linnaean Society*, 7: 27-43

Thorpe, R. S. (1976) Biometric analysis of geographic variation and racial affinities. *Biological Reviews*, 51: 407-452

Thorpe, R. S. (1983) A review of the numerical methods for recognising and analysing racial differentiation, in: Felsenstein, J. (Ed.), *Numerical Taxonomy*. Springer Verlag, Berlin, pp. 404-423

Reist, J. D. (1985) An empirical evaluation of several univariate methods that adjust for size variation in morphometric data. *Canadian Journal of Zoology*, 63: 1429-1439

Lleonart, J., Salat, J., and Torres, G. J. (2000) Removing Allometric Effects of Body Size in Morphological Analysis. *Journal of Theoretical Biology*, 205: 85-93

