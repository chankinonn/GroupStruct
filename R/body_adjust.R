#' Allometric adjustment of morphometric characters to correct for intraspecific
#' ontogenetic variation.
#'
#' This function assumes that each Operational Taxonomic
#' Unit represents a distinct species. The input file is data frame in csv format. The firt column contains OTU identifiers, followed by
#' morphological characters. No other data should be included in the data frame.
#' Singleton OTU's and juvenile measurements should be exlucded.
#'
#'
#' @param data Path to the input file
#' @return Returns log-transformed and adjusted data
#' @import dplyr
#' @export

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

}

