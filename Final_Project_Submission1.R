
library(tidyr)
library(tidyverse)
library(ggplot2)
library(colorspace)

getwd()

df_genes <- read.csv("~/Desktop/QBS_103/Final_Project/QBS103_GSE157103_genes.csv")
df_matrix <- read.csv("~/Desktop/QBS_103/Final_Project/QBS103_GSE157103_series_matrix-1.csv")

head(rownames(df_genes))
head(df_genes[,1])


base_green <- "#ADCCA3"
base_brown <- "#D2B48C"
pastel_olive <- darken(base_green, amount = 0.2) 
head(df_genes)
colnames(df_genes)[1] <- "Genes" # Names first column of df_genes "Genes"
gene_ABCD3 <- filter(df_genes, df_genes["Genes"]== "ABCD3") # Pulls the row for gene ABCD3 into a data frame by filtering.
hist_values <- as.numeric(gene_ABCD3[, -1]) # I asked chatgpt for help here. I could not figure out how to pull rows, not columns. 
  
hist(hist_values, breaks = 15, 
     main = "ABCD3 Expression", 
     xlab = "Expression Values",
     ylim = c(0, 20), xlim = c(0, 20),
     col= pastel_olive)


samples <- colnames(df_genes)[-1] # Pull column names into a data frame.
df_ABCD3 <- data.frame(participant_id = samples, Expression = hist_values) # Initialize new data frame with expression values and charlson_score, for scatterplot.
# I used chat.dartmout.edu for help on merging the data frames. It suggested the merge() function.
merged_df <- merge(df_ABCD3, df_matrix, by = "participant_id") # Merge Data frames by "participant_id".
head(merged_df)

## Scatter Plot
ggplot(merged_df, aes(merged_df$Expression, merged_df$charlson_score)) +
  geom_point() +
  theme_classic() +
  labs(title="ABCD3 Gene Expression vs. Charlson Score", x = "Expression", y = "Charlson Score")
  








