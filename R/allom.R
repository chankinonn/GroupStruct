#' Allometric adjustment of morphometric characters to correct for
#' ontogenetic variation in multispecies or multipopulation datasets
#'
#' A multispecies dataset should comprise more than 1 species where each unique identifier in the first column represents
#' a different species. If there are several populations/localities within a species, these should be grouped under a common identifier (pooled groups).
#' A multipopulation dataset comprises a single species where each unique identifier represents
#' a different population/locality. Datasets should be in csv format with the 1st colum = OTU identifier; 2nd column = body size measurement; 3rd column onwards = measurements of other morphometric traits
#'
#' This function uses the following allometric growth equation to adjust for ontogenetic variation:
#' Xadj = log10(X)-B[log10(BL)-log10(BLmean)], where Xadj = Adjusted value for character X;
#' X = raw/unadjusted value for character X; B = pooled regression coefficient (slope) of log10(X) against log10(BL);
#' BL = observed value of the standard measure of body length (e.g. snout-vent-length); BLmean = grand mean for BL.
#'
#' There are two ways to calculate the slope depending on how groups are structured: (i) Pooled groups combine different localities/populations to form a single
#' compound locality to boost sample size. However, this should only be done if populations are not significantly heterogeneous (e.g., no geographic variation).
#' (ii) Common within-group pooling considers each population separately and a different slope is calculated for each population.
#' For type = population1 in a multipopulation dataset, a separate slope is
#' calculated for each population (common within-group pooling). If type = population2, all populations are pooled (pooled group)
#' and a single slope is calculated based on combined populations. For multispecies datasets (type = species), populations
#' are pooled under a single species identifier and slopes are calculated based on pooled groups. If population-level variation is present within the multispecies dataset,
#' each species should be analyzed as separate multipopulation datasets using type = population1.
#' The other difference in calculations between multispecies and
#' multipopulation datasets is in BLmean. For multispecies datasets, a separate BLmean will be calculated for each species,
#' whereas in a multipopulation dataset, a single BLmean is calculated to represent the standard size of the species averaged across all populations.
#'
#'
#' Each OTU should be represented by more than two individuals for the adjustment to work and missing data must not be included.
#' This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism.
#'
#' @param data csv file with species identifiers in the 1st column, followed by body length measurement (e.g. snout-vent-length) in the 2nd column. Other morphometric measurements are contained in the 3rd column onwards. Singletons (only one sample per species) and missing data are not allowed
#' @param type Enter "species" for multispecies and "population" for multipopulation datasets.
#' @return Returns log-transformed and body-size adjusted data in a table called outfile.csv
#' @import dplyr
#' @export
#' @examples
#' For multispecies datasets:
#' m <- read.csv("foo.csv") ## foo.csv = raw measurements
#' allom(m, type="species")
#'
#' For multipopulation datasets (slope calculated based on common within-group pooling):
#' allom(m, type="population1")
#'
#' For multipopulation datasets (slope calculated based pooled groups):
#' allom(m, type="population2")
#'
#' @references
#' Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to
#' intraspecific variation in the Ringed Snake Natrix natrix (L.). Biological Journal of the Linnaean Society, 7: 27-43

allom <- function(data, type){
  # Reorder input data alphabetically by first column
  data <- data[order(data[,1]), ]

  species_list <- as.matrix(unique(data[,1]))  # Create list of unique species
  finalmatrix <- vector("list", length(species_list))  # Preallocate list for results

  # Loop through subsets of species and perform calculations
  for (i in seq_along(species_list)) {
    species_subsets <- filter(data, data[,1] == species_list[i])

    # Define functions for allometric adjustment
    allo <- function(x) {
      y <- species_subsets[,2]
      temp <- lm(log10(x) ~ log10(y), species_subsets)  # Linear regression
      beta <- as.numeric(temp$coefficients[2])  # Extract slope
      log10(x) - beta * (log10(y) - log10(mean(y)))  # Apply equation
    }

    allo_pop1 <- function(x) {
      y <- species_subsets[,2]
      z <- mean(data[,2])
      temp <- lm(log10(x) ~ log10(y), species_subsets)
      beta <- as.numeric(temp$coefficients[2])
      log10(x) - beta * (log10(y) - log10(z))
    }

    allo_pop2 <- function(x) {
      y <- data[,2]
      z <- mean(data[,2])
      temp <- lm(log10(x) ~ log10(y), data)
      beta <- as.numeric(temp$coefficients[2])
      log10(x) - beta * (log10(y) - log10(z))
    }

    # Apply the appropriate function based on 'type'
    if (type == "species") {
      finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo)
    } else if (type == "population1") {
      finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo_pop1)
    } else if (type == "population2") {
      finalmatrix[[i]] <- apply(data[,3:ncol(data)], 2, allo_pop2)
    }
  }

  # Combine output, log-transform, and write to table
  logsize <- log10(data[,2])
  all_combined <- cbind(logsize, do.call(rbind, finalmatrix))  # Bind results
  final_adjusted <- data.frame(cbind(data[,1], all_combined))  # Bind species column

  # Preserve original column names
  colnames(final_adjusted) <- colnames(data)

  # Write output
  write.csv(final_adjusted, "allom_outfile.csv", row.names = FALSE)
  print(final_adjusted)
}
