% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdist_compare.R
\name{pdist_compare}
\alias{pdist_compare}
\title{Create boxplots comparing p-distances of a focal species vs. a list of other species}
\usage{
pdist_compare(
  raw_distance_df,
  focal_species,
  comparison_species,
  colors = NULL,
  output_table
)
}
\arguments{
\item{raw_distance_df}{Results from the pdist_cutoff function}

\item{focal_species}{Focal species that you want to make comparisons with}

\item{comparison_species}{A list of other species to be compared with the focal species}

\item{colors}{Optional. List of manual colors for each comparison}

\item{output_table}{Name of output table}
}
\value{
A boxplot comparing the focal species with the supplied list of comparison species. Also writes a pairwise comparison table with min-max p-distances
}
\description{
Uses output from the pdist_cutoff function
}
\examples{
Set path to fasta alignment
myfasta <- "path/to/fasta/fasta_alignment.fas"

Set cutoff value
mycutoff <- 0.05

Set your output file name
raw_output <- "raw_output.csv"
aggregated_output <- "aggregated_output.csv"

Get p-distances
results <- pdist_cutoff(myfasta, mycutoff, raw_output, aggregated_output)

Specify the focal species, a list of other species to compare, name of output table, and optional plotting colors
focal_species <- "Species_A"
comparison_species <- c("Species_B", "Species_C", "Species_D")
output_table <- "comparison_table.csv"
manual_colors <- c("Species_B" = "red", "Species_C" = "blue", "Species_D" = "green")

Plot
pdist_compare(results, focal_species, comparison_species, colors = manual_colors, output_table)
}
