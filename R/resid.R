#' Calculate residuals of a linear regression of each character against body size
#'
#'
#' @param data csv file with species identifiers in the 1st column, followed by body length measurement (e.g. snout-vent-length) in the 2nd column. Other morphometric measurements are contained in the 3rd column onwards. Singletons (only one sample per species) and missing data are not allowed
#' @return Returns residuals (slope) for each character
#' @import dplyr
#' @export
#' @examples
#' m <- read.csv("foo.csv") ## foo.csv = raw measurements
#' resid(m)


resid <- function(data){

  species_list <- as.matrix(unique(data[1]))  ## Create list of unique species
  finalmatrix <- NULL  ### Create empty matrix and start loop for calculation

  ### Loop through subsets of species and perform calculations
  for (i in 1:nrow(species_list)) {

    ### Create species subsets
    species_subsets <- filter(data, data[,1]==species_list[i])

    ### Define fx for allometric equation as defined above
    residuals <- function(x){
      y <- species_subsets$SVL
      beta <- lm(log10(x)~log10(y), species_subsets)  ### Perform linear regression
      beta <- as.numeric(beta$coefficients[2]) ### Extract beta coefficient
      #adj <- log10(x)-temp*(log10(y)-log10(mean(y)))  ### Plug into equation
    }

    ### Loop through each species subset and apply fx from column 3 onwards; store results in "final"
    finalmatrix[[i]] <- apply(species_subsets[,3:ncol(species_subsets)], 2, residuals)
  }
  ### Combine output, log-transform, and write to table
  all_combined <- cbind(SVL=log10(data$SVL), do.call(rbind, finalmatrix)) ### Combine results from loop and bind the column "SVL"
  final_adjusted <- data.frame(cbind(Species=data[1], all_combined)) ### Bind the column "Species" to the final dataset
  write.csv(final_adjusted, "allom_outfile.csv", row.names = FALSE)
  print(final_adjusted)
}

