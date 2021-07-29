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
#' In a multipopulation dataset, there are two ways of calculating the slope. For type = population1, a separate slope is
#' calculated for each population (common within-group pooling). If type = population2, all populations are pooled (pooled group)
#' and a single slope is calculated based on combined populations. The latter should be used if the slopes of each population have
#' been ascertained to be equal (Thorpe 1976). The other difference in calculations between multispecies and
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

  species_list <- as.matrix(unique(data[1]))  ## Create list of unique species
  finalmatrix <- NULL  ### Create empty matrix and start loop for calculation

  ### Loop through subsets of species and perform calculations
  for (i in 1:nrow(species_list)) {

    ### Create species subsets
    species_subsets <- filter(data, data[,1]==species_list[i])

    ### Define fx for allometric equation as defined above
    allo <- function(x){
      y <- species_subsets[,2]
      temp <- lm(log10(x)~log10(y), species_subsets)  ### Perform linear regression
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log10(x)-temp*(log10(y)-log10(mean(y)))  ### equation for interspecies
    }

    allo_pop1 <- function(x){
      y <- species_subsets[,2]
      z <- mean(data[,2])
      temp <- lm(log10(x)~log10(y), species_subsets)  ### Perform linear regression, slope=common within groups
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log10(x)-temp*(log10(y)-log10(z))  ### Plug into equation
    }

    allo_pop2 <- function(x){
      y <- data[,2]
      z <- mean(data[,2])
      temp <- lm(log10(x)~log10(y), data)  ### Perform linear regression, slope=pooled groups
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log10(x)-temp*(log10(y)-log10(z))  ### Plug into equation
    }

    ### Loop through each species subset and apply fx from column 3 onwards; store results in "final"
    if(type == "species"){

      finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo)

    }

    if(type =="population1"){
      finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo_pop1)

    }

    if(type =="population2"){
      finalmatrix[[i]] <- apply(data[,3:ncol(data)], 2, allo_pop2)

    }
  }

  ### Combine output, log-transform, and write to table
  logsize <- log10(data[2])
  all_combined <- cbind(logsize, do.call(rbind, finalmatrix)) ### Combine results from loop and bind the column "SVL"
  final_adjusted <- data.frame(cbind(data[1], all_combined)) ### Bind the column "Species" to the final dataset
  write.csv(final_adjusted, "allom_outfile.csv", row.names = FALSE)
  print(final_adjusted)
}
