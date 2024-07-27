#' Perform allometric size correction, PCA, and ANOVA/TUKEY analyses on morphological data
#'
#' Wrapper function that takes raw morphological data, performs allometric size correction using `allom()`, followed by PCA and ANOVA-Tukey test on size-corrected data.
#'

#' @param input_file Raw morphological data in CSV format. See `allom()` documention on how to format input data
#' @param type Type of allometric correction to use. Choose from "species", "population1" or "population2". See `allom()` documention
#' @param colors Optional vector of colors for the PCA plot. Number of colors need to equal number of groups/species
#' @param draw_ellipses Boolean. Whether to draw confidence ellipse on PCA plot
#' @param confidence_level Level of confidence interval to draw ellipse
#' @return PCA plot of PC1 and PC2 scores. Also writes two tables. One containing a summary of the PCA analysis and another containing the results of the ANOVA/Tukey analysis
#' @import ggplot2 dplyr multcomp broom
#' @export
#' @examples
#' Path to raw morphological data
#' mydata <- "path/to/morpho_data.csv"
#'
#' morpho_struct(mydata)

morpho_struct <- function(input_file, type = "species", colors = NULL, draw_ellipses = TRUE, confidence_level = 0.95) {
  # Read the input file as a CSV
  input_data <- read.csv(input_file)

  # Extract the species column
  species <- input_data[, 1]

  # Execute the allom() function from GroupStruct package
  allom_result <- allom(input_data, type)

  # Check if the result contains the data required for PCA
  if (!is.null(allom_result) && nrow(allom_result) > 0) {
    # Perform PCA using prcomp (excluding the first column which is species)
    pca_result <- prcomp(allom_result[, -1], scale. = TRUE)

    # Calculate the percentage of variance explained
    pca_variance <- summary(pca_result)$importance[2, ]
    x_variance <- round(pca_variance[1] * 100, 2)
    y_variance <- round(pca_variance[2] * 100, 2)

    # Create a data frame for ggplot2
    pca_data <- data.frame(
      Sample = species,
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      Group = as.factor(species)
    )

    # Generate the PCA plot using ggplot2
    pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3) +
      labs(x = paste("PC 1 (", x_variance, "%)", sep = ""),
           y = paste("PC 2 (", y_variance, "%)", sep = "")) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "gray", linetype = "dotted"),
        axis.line = element_line(color = "black"),
        panel.border = element_blank()
      )

    # Draw ellipses if requested
    if (draw_ellipses) {
      pca_plot <- pca_plot + stat_ellipse(aes(group = Group), level = confidence_level)
    }

    # Apply manual colors if provided
    if (!is.null(colors)) {
      pca_plot <- pca_plot + scale_color_manual(values = colors)
    }

    # Print the PCA plot
    print(pca_plot)

    # Create an output table with loadings, summary stats, and eigenvalues
    loadings <- pca_result$rotation
    temp <- summary(pca_result)
    sum <- temp$importance
    eigen <- pca_result$sdev^2
    output_table <- rbind(sum, eigen, loadings)

    # Write the output table to a CSV file
    write.csv(output_table, file = "output_pca.csv")

    # Perform ANOVA and Tukey's HSD test
    anova_results <- lapply(names(allom_result)[-1], function(var) {
      model <- aov(as.formula(paste(var, "~ species")), data = allom_result)
      tukey <- TukeyHSD(model)
      tidy(tukey) %>%
        mutate(Variable = var)  # Add variable name to each row
    })

    # Combine the results into a single data frame
    anova_summary <- bind_rows(anova_results)

    # Write all columns of the Tukey HSD results to a CSV file
    write.csv(anova_summary, file = "output_anova.csv")

    # Return the PCA results, plot, and ANOVA summary for further analysis if needed
    return(list(pca_result = pca_result, pca_plot = pca_plot, output_table = output_table, anova_summary = anova_summary))
  } else {
    stop("The result from allom() does not contain the required data for PCA.")
  }
}
