#' Allometric adjustment of morphometric characters to correct for
#' ontogenetic variation in a multipopulation (single species) dataset
#'
#' This function uses the following allometric growth equation to adjust for ontogenetic variation:
#' Xadj = log10(X)-B[log10(SVL)-log10(SVLmean)], where Xadj = Adjusted value for character X;
#' X = raw/unadjusted value for character X; B = pooled regression coefficient (slope) of log10(X) against log10(SVL);
#' SVL = measured SVL; SVLmean = grand mean SVL for a particular species
#'
#' Following Thorpe(1975), this function uses common within-group pooling to obtain a separate slope for each population.
#' SVLmean is the grand mean calculated by averaging over all populations. 
#' Each species should be represented by more than two individuals for the adjustment to work and missing data must not be included. 
#' This adjustment should also be performed separately on different sexes to account for possible sexual dimorphism. 
#
#'
#' @param data csv file with population identifiers in the 1st column, followed by body length measurement (e.g. snout-vent-length) in the 2nd column. Other morphometric measurements are contained in the 3rd column onwards. Singletons (only one sample per species) and missing data are not allowed
#' @return Returns log-transformed and body-size adjusted data in a table called outfile.csv
#' @import dplyr
#' @export
#' @examples
#' m <- read.csv("foo.csv") ## foo.csv = raw measurements
#' allom_pop(m)
#' @references
#' Thorpe, R. S. (1975) Quantitative handling of characters useful in snake systematics with particular reference to
#' intraspecific variation in the Ringed Snake Natrix natrix (L.). Biological Journal of the Linnaean Society, 7: 27-43

allom_pop <- function(data){
  
  pop_list <- as.matrix(unique(data[1]))  ## Create list of unique species
  finalmatrix <- NULL  ### Create empty matrix and start loop for calculation
  
  ### Loop through subsets of species and perform calculations
  for (i in 1:nrow(pop_list)) {
    
    ### Create species subsets
    pop_subsets <- filter(data, data[,1]==pop_list[i])
    
    ### Define fx for allometric equation as defined above
    allo_pop <- function(x){
      y <- pop_subsets$SVL
      z <- mean(data$SVL)
      temp <- lm(log10(x)~log10(y), pop_subsets)  ### Perform linear regression
      temp <- as.numeric(temp$coefficients[2]) ### Extract beta coefficient
      adj <- log10(x)-temp*(log10(y)-log10(z))  ### Plug into equation
    }
    
    ### Loop through each species subset and apply fx from column 3 onwards; store results in "final"
    finalmatrix[[i]] <- apply(pop_subsets[,3:ncol(pop_subsets)], 2, allo_pop)
  }
  
  ### Combine output, log-transform, and write to table
  all_combined <- cbind(SVL=log10(data$SVL), do.call(rbind, finalmatrix)) ### Combine results from loop and bind the column "SVL"
  final_adjusted <- data.frame(cbind(Species=data[1], all_combined)) ### Bind the column "Species" to the final dataset
  write.csv(final_adjusted, "allom_pop_outfile.csv", row.names = FALSE)
  print(final_adjusted)
}

