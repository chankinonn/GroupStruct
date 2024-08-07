% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdist_cutoff.R
\name{pdist_cutoff}
\alias{pdist_cutoff}
\title{Calculate p-distances and identify species pairs that are above a user-specified threshold}
\usage{
pdist_cutoff(fasta_file, cutoff, raw_output_file, aggregated_output_file)
}
\arguments{
\item{fasta_file}{Alignned FASTA sequences. Make sure sequences are labeled as genus_species. You can add other unique identifiers after the species name using a second underscore. E.g. Bufo_bufo_MN09667}

\item{cutoff}{Set your desired cut-off. E.g. 0.03 for a 3 percent cut-off}

\item{raw_output_file}{Name of output file for raw p-distances}

\item{aggregated_output_file}{Name of output file for p-distances aggregated by species}
}
\value{
Writes two CSV-formatted files. One contains raw p-distances for every individual pairwise comparison. A second table aggregates p-distances by species showing min, max, and mean. Species that meet the cut-off will be printed
}
\description{
Uses the dist_dna function in ape to calculate p-distances
}
\examples{
Set path to fasta alignment
myfasta <- "path/to/fasta/fasta_alignment.fas"

Set cutoff value
mycutoff <- 0.05

Set your output file name
raw_output <- "raw_output.csv"
aggregated_output <- "aggregated_output.csv"

Run example
pdist_cutoff(myfasta, mycutoff, raw_output, aggregated_output)
}
