#' Create boxplots comparing p-distances of a focal species vs. a list of other species
#'
#' Uses output from the pdist_cutoff function
#'

#' @param x Output from the pdist_cutoff function
#' @param focal_species Focal species that you want to make comparisons with
#' @param comparison_species A list of other species to be compared with the focal species
#' @return Writes two CSV-formatted files. One contains raw p-distances for every individual pairwise comparison. A second table aggregates p-distances by species showing min, max, and mean. Species that meet the cut-off will be printed
#' @import ape tidyr dplyr
#' @export
#' @examples
#' Set path to fasta alignment
#' myfasta <- "path/to/fasta/fasta_alignment.fas"
#'
#' Set cutoff value
#' mycutoff <- 0.05
#'
#' Set your output file name
#' raw_output <- "raw_output.csv"
#' aggregated_output <- "aggregated_output.csv"
#'
#' Get p-distances
#' results <- pdist_cutoff(myfasta, mycutoff, raw_output, aggregated_output)
#'
#' Specify the focal species  and the list of comparison species
#' focal_species <- "Micryletta_subaraji"
#' comparison_species <- c("Micryletta_sumatrana", "Micryletta_dissimulans", "Micryletta_inornata", "Micryletta_nigromaculata")
#'
#' Boxplots for the specified species combinations
#' plot_pdist(results, focal_species, comparison_species)

plot_pdist <- function(raw_distance_df, focal_species, comparison_species) {
  # Filter the results to include only rows where the specified species is involved in the specified combinations
  plot_data <- raw_distance_df %>%
    filter((Species1 == focal_species & Species2 %in% comparison_species) |
             (Species2 == focal_species & Species1 %in% comparison_species)) %>%
    mutate(Comparison_Species = ifelse(Species1 == focal_species, Species2, Species1))

  # Plot the boxplot
  ggplot(plot_data, aes(x = Comparison_Species, y = Freq)) +
    geom_boxplot() +
    labs(title = paste("P-Distances for", focal_species, "and specified comparisons"),
         x = "Comparison Species",
         y = "P-Distance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
