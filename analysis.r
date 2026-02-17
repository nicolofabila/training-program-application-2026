# ---------------------------------------------------------

# Melbourne Bioinformatics Training Program

# This exercise to assess your familiarity with R and git. Please follow
# the instructions on the README page and link to your repo in your application.
# If you do not link to your repo, your application will be automatically denied.

# Leave all code you used in this R script with comments as appropriate.
# Let us know if you have any questions!


# You can use the resources available on our training website for help:
# Intro to R: https://mbite.org/intro-to-r
# Version Control with Git: https://mbite.org/intro-to-git/

# ----------------------------------------------------------

# Load libraries -------------------
# You may use base R or tidyverse for this exercise

library(tidyverse)

# Load data here ----------------------
# Load each file with a meaningful variable name.

expression_data <- read.csv("data/GSE60450_GeneLevel_Normalized(CPM.and.TMM)_data.csv", 
                            check.names = FALSE)
metadata <- read.csv("data/GSE60450_filtered_metadata.csv", 
                     check.names = FALSE)



# Inspect the data -------------------------

# What are the dimensions of each data set? (How many rows/columns in each?)
# Keep the code here for each file.

## Expression data
dim(expression_data)

## Metadata
dim(metadata)


# Prepare/combine the data for plotting ------------------------
# How can you combine this data into one data.frame?

# Transform expression data from wide to long format

colnames(expression_data)[1] <- "gene_id"

# Pivot expression data to long format using pivot_longer()
# gene_id and gene_symbol as identifier columns, pivot all samples
expression_long <- expression_data %>%
  pivot_longer(cols = starts_with("GSM"), 
               names_to = "sample_id", 
               values_to = "expression") %>%
  select(gene_id, gene_symbol, sample_id, expression)

# Merge with metadata using sample_id (first column of metadata)
colnames(metadata)[1] <- "sample_id"

# Use left join to Combine expression and metadata
combined_data <- expression_long %>%
  left_join(metadata, by = "sample_id")



# Plot the data --------------------------
## Plot the expression by cell type
## Can use boxplot() or geom_boxplot() in ggplot2


expression_plot <- ggplot(combined_data, 
                         aes(x = immunophenotype, 
                             y = expression, 
                             fill = immunophenotype)) +
  geom_boxplot() +
  scale_y_log10() +  # log scale for better visualization of expression data
  labs(title = "Gene Expression by Cell Type",
       x = "Cell Type",
       y = expression(Expression~(log[10]~CPM/TMM)),
       fill = "Cell Type") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "none")



## Save the plot
### Show code for saving the plot with ggsave() or a similar function

ggsave("results/expression_by_celltype_boxplot.png", 
       plot = expression_plot,
       width = 7,
       height = 5,
       dpi = 300)