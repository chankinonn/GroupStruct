#' Allometric adjustment of morphometric characters to correct for intraspecific
#' ontogenetic variation
#'
#' This function uses an allometric growth equation to adjust for ontogenetic variation
#' among populations within a species, i.e. intraspecific variation. Hence, this function
#' should be applied separately for different species. The function uses the following allometric
#' equation: Xadj = log(X)-B[log(SVL)-log(SVLmean)], where Xadj
#' = Adjusted value for character X; X = raw/unadjusted value for character X; B =
#' unstandardized regression coefficient for X against SVL; SVL = measured SVL; SVLmean =
#' mean SVL for all populations.
#'
#' Note that each population will have its own B and SVLmean is the average SVL of all populations
#' included in the dataset.
#'
#' This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism.
#'
#' @param data csv file with population identifiers in the 1st column, followed by body legnth measurement (e.g. snout-vent-length) in the 2nd column. Other morphometric measurements are contained in the 3rd column onwards. Singletons (only one sample per population) and missing data are not allowed
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
      temp <- lm(x~y, species_subsets)  ### Perform linear regression
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log(x)-temp*(log(y)-log(mean(y)))  ### Plug into equation
    }

    ### Loop through each species subset and apply fx from column 3 onwards; store results in "final"
    finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, allo)
  }

  ### Combine output, log-transform, and write to table
  all_combined <- cbind(SVL=log(data$SVL), do.call(rbind, finalmatrix)) ### Combine results from loop and bind the column "SVL"
  final_adjusted <- data.frame(cbind(Species=data[1], all_combined)) ### Bind the column "Species" to the final dataset
  write.csv(final_adjusted, "outfile.csv", row.names = FALSE)
  print(final_adjusted)
}

