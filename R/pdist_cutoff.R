#' Calculate p-distances and identify species pairs that are above a user-specified threshold
#'
#' Uses the dist_dna function in ape to calculate p-distances
#'

#' @param fasta_file Alignned FASTA sequences. Make sure sequences are labeled as genus_species. You can add other unique identifiers after the species name using a second underscore. E.g. Bufo_bufo_MN09667
#' @param cutoff Set your desired cut-off. E.g. 0.03 for a 3 percent cut-off
#' @param raw_output_file Name of output file for raw p-distances
#' @param aggregated_output_file Name of output file for p-distances aggregated by species
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
#' Run example
#' pdist_cutoff(myfasta, mycutoff, raw_output, aggregated_output)

pdist_cutoff <- function(fasta_file, cutoff, raw_output_file, aggregated_output_file) {
  # Read the aligned sequences from the FASTA file
  sequences <- read.dna(fasta_file, format = "fasta")

  # Calculate the p-distances using the complete deletion method
  p_distances <- dist.dna(sequences, model = "raw", pairwise.deletion = FALSE)

  # Convert the distance matrix to a data frame for easier manipulation
  distance_df <- as.data.frame(as.table(as.matrix(p_distances)))

  # Convert Var1 and Var2 to character vectors
  distance_df$Var1 <- as.character(distance_df$Var1)
  distance_df$Var2 <- as.character(distance_df$Var2)

  # Extract species names from sequence names
  distance_df <- distance_df %>%
    mutate(Species1 = sapply(strsplit(Var1, "_"), function(x) paste(x[1:2], collapse = "_")),
           Species2 = sapply(strsplit(Var2, "_"), function(x) paste(x[1:2], collapse = "_")))

  # Save the original p-distance results to a CSV file
  write.csv(distance_df, file = raw_output_file, row.names = FALSE)
  cat("Raw p-distance results saved to", raw_output_file, "\n")

  # Filter pairs where the value exceeds the cutoff
  filtered_pairs <- distance_df %>%
    filter(Freq > cutoff) %>%
    select(Species1, Species2, Freq)

  # Aggregate results by species pairs and summarize the distance ranges
  aggregated_results <- filtered_pairs %>%
    group_by(Species1, Species2) %>%
    summarize(Min_Distance = min(Freq),
              Max_Distance = max(Freq),
              Median_Distance = median(Freq)) %>%
    arrange(Species1, Species2)

  # Save the aggregated results to a CSV file
  write.csv(aggregated_results, file = aggregated_output_file, row.names = FALSE)
  cat("Aggregated results saved to", aggregated_output_file, "\n")

  # Extract unique species names from the first two columns
  unique_species <- unique(c(aggregated_results$Species1, aggregated_results$Species2))

  # Print the unique species names
  cat("Species that meet the cutoff:\n")
  print(unique_species)

  return(distance_df)  # Return the raw p-distance data frame
}
