#' Create boxplots comparing p-distances of a focal species vs. a list of other species
#'
#' Uses output from the pdist_cutoff function
#'

#' @param x Output from the pdist_cutoff function
#' @param focal_species Focal species that you want to make comparisons with
#' @param comparison_species A list of other species to be compared with the focal species
#' @param output_table Name of output table
#' @param colors Optional. List of manual colors for each comparison
#' @return A boxplot comparing the focal species with the supplied list of comparison species. Also writes a pairwise comparison table with min-max p-distances
#' @import ape tidyr dplyr ggplot2
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
#' Specify the focal species, a list of other species for comparison, name for output table, and optional plotting colors
#' focal_species <- "Species_A"
#' comparison_species <- c("Species_B", "Species_C", "Species_D")
#' output_table <- "comparison_table.csv"
#' manual_colors <- c("Species_B" = "red", "Species_C" = "blue", "Species_D" = "green")
#'
#' Boxplots for the specified species combinations
#' pdist_compare(results, focal_species, comparison_species, colors = manual_colors, output_table)

pdist_compare <- function(raw_distance_df, focal_species, comparison_species, colors = NULL, output_table) {
  # Filter the results to include only rows where the specified species is involved in the specified combinations
  plot_data <- raw_distance_df %>%
    filter((Species1 == focal_species & Species2 %in% comparison_species) |
             (Species2 == focal_species & Species1 %in% comparison_species)) %>%
    mutate(Comparison_Species = ifelse(Species1 == focal_species, Species2, Species1))

  # Plot the boxplot
  plot <- ggplot(plot_data, aes(x = Comparison_Species, y = Freq, fill = Comparison_Species)) +
    geom_boxplot() +
    labs(title = paste("P-Distances for", focal_species, "and specified comparisons"),
         x = "Comparison Species",
         y = "P-Distance") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (!is.null(colors)) {
    plot <- plot + scale_fill_manual(values = colors)
  }
  print(plot)

  # Create a pairwise comparison table with min-max values
  species_list <- unique(c(focal_species, comparison_species))
  comparison_matrix <- matrix("", nrow = length(species_list), ncol = length(species_list),
                              dimnames = list(species_list, species_list))

  for (i in seq_along(species_list)) {
    for (j in seq_along(species_list)) {
      if (i != j) {
        comparison_values <- raw_distance_df$Freq[
          (raw_distance_df$Species1 == species_list[i] & raw_distance_df$Species2 == species_list[j]) |
            (raw_distance_df$Species1 == species_list[j] & raw_distance_df$Species2 == species_list[i])
        ]
        if (length(comparison_values) > 0) {
          min_value <- round(min(comparison_values), 2)
          max_value <- round(max(comparison_values), 2)
          comparison_matrix[i, j] <- paste0(min_value, "-", max_value)
        }
      }
    }
  }

  comparison_df <- as.data.frame(comparison_matrix)

  # Save the pairwise comparison table to a CSV file
  write.csv(comparison_df, file = output_table, row.names = TRUE)
  cat("Pairwise comparison table saved to", output_table, "\n")
}
