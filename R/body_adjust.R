#' Allometric adjustment of morphometric characters to correct for
#' ontogenetic variation in a multispecies dataset
#'
#' This function uses the following allometric growth equation to adjust for ontogenetic variation:
#' Xadj = log10(X)-B[log10(SVL)-log10(SVLmean)], where Xadj = Adjusted value for character X;
#' X = raw/unadjusted value for character X; B = pooled regression coefficient (slope) of log10(X) against log10(SVL);
#' SVL = measured SVL; SVLmean = grand mean SVL for a particular species
#'
#' Following Thorpe(1975), pooling refers to the grouping of multiple localities (of the same species)
#' to produce a compound locality with a higher sample size. This produces a more accurate slope compared to
#' slopes derived from subpopulations with low sample sizes. However, if population-level sampling is adequate
#' and the slopes of different populations are not equal (e.g. geographic variation), the common within-group
#' slope should be used (separate slope for each population)
#'
#' This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism.
#'
#' @param data csv file with species identifiers in the 1st column, followed by body length measurement (e.g. snout-vent-length) in the 2nd column. Other morphometric measurements are contained in the 3rd column onwards. Singletons (only one sample per species) and missing data are not allowed
#' @return Returns log-transformed and body-size adjusted data in a table called outfile.csv
#' @import dplyr
#' @export
#' @examples
#' m <- read.csv("foo.csv") ## foo.csv = raw measurements
#' body_adjust(m)

body_adjust <- function(data){

  species_list <- as.matrix(unique(data[1]))  ## Create list of unique species
  finalmatrix <- NULL  ### Create empty matrix and start loop for calculation

  ### Loop through subsets of species and perform calculations
  for (i in 1:nrow(species_list)) {

    ### Create species subsets
    species_subsets <- filter(data, data[,1]==species_list[i])

    ### Define fx for allometric equation as defined above
    allo <- function(x){
      y <- species_subsets$SVL
      temp <- lm(log10(x)~log10(y), species_subsets)  ### Perform linear regression
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log10(x)-temp*(log10(y)-log10(mean(y)))  ### Plug into equation
    }

    ### Loop through each species subset and apply fx from column 3 onwards; store results in "final"
    finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo)
  }

  ### Combine output, log-transform, and write to table
  all_combined <- cbind(SVL=log10(data$SVL), do.call(rbind, finalmatrix)) ### Combine results from loop and bind the column "SVL"
  final_adjusted <- data.frame(cbind(Species=data[1], all_combined)) ### Bind the column "Species" to the final dataset
  write.csv(final_adjusted, "outfile.csv", row.names = FALSE)
  print(final_adjusted)
}

