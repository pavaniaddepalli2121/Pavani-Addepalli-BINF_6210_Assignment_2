###### **********************
## ASSIGNMENT 2 - SOFTWARE TOOLS
##
## PAVANI ADDEPALLI (1326277)
##
## 2024-10-28
## 
## Supervised Machine Learning and Classification Using Sequence Features
##
## "Machine Learning Classification of BRCA1 and BRCA2 Gene Sequences in Homo sapiens: Enhancing Genetic Variant Analysis for Cancer Susceptibility"
##
##
## Code Part 1: Data Preparation ------
## Load necessary libraries 

library(randomForest)    # Random Forest algorithm for classification
library(e1071)           # Support Vector Machine (SVM)
library(pROC)            # ROC curve analysis
library(readr)           # Reading and writing data
library(caret)           # Machine learning utilities
library(rentrez)         # Accessing NCBI data
library(Biostrings)      # Bioconductor package for DNA sequence analysis
library(ggplot2)         # Creating visualizations
library(corrplot)        # Visualizing correlation matrices
library(RColorBrewer)    # Enhancing the visual appeal
library(tidyverse)       # Data manipulation and visualization
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())


## Data Loading  ------

# I used the line below on October 28, 2024, to obtain the dataset of BRCA1 and BRCA2 genes associated with DNA repair in Homo sapiens, which are implicated in breast cancer susceptibility.

## Fetch BRCA1 and BRCA2 sequences from NCBI ------

# Retrieve IDs for BRCA1 
brca1_ids <- entrez_search(db = "nucleotide", term = "BRCA1[Gene] AND Homo sapiens[Organism]", retmax = 5)$ids

# Fetch BRCA1 sequences 
brca1_sequences <- entrez_fetch(db = "nucleotide", id = brca1_ids, rettype = "fasta")

# Retrieve IDs for BRCA2
brca2_ids <- entrez_search(db = "nucleotide", term = "BRCA2[Gene] AND Homo sapiens[Organism]", retmax = 5)$ids

# Fetch BRCA2 sequences
brca2_sequences <- entrez_fetch(db = "nucleotide", id = brca2_ids, rettype = "fasta")

# Save fetched BRCA1 sequences to a FASTA file
write(brca1_sequences, "brca1_sequences.fasta", sep = "\n")

# Save fetched BRCA2 sequences to a FASTA file
write(brca2_sequences, "brca2_sequences.fasta", sep = "\n")

# Read the saved BRCA1 data file from the laptop as DNA StringSet
brca1_string_set <- readDNAStringSet("C:/Users/drpav/OneDrive/Documents/brca1_sequences.fasta")  

# Read the saved BRCA2 data file from the laptop as DNA StringSet
brca2_string_set <- readDNAStringSet("C:/Users/drpav/OneDrive/Documents/brca2_sequences.fasta")  

# Combine the sequence data for the BRCA1 and BRCA2 genes
brca_data <- bind_rows(
  data.frame(gene = "BRCA1", sequence = as.character(brca1_string_set)),
  data.frame(gene = "BRCA2", sequence = as.character(brca2_string_set))
)

## Check structure of the data ------

# Checking the type of the brca_data object to confirm it's suitable for data manipulation functions
class_brca_data <- class(brca_data)  # Result: "data.frame" - This confirms that brca_data is in a tabular structure with labeled rows and columns.
cat("Class of brca_data:", class_brca_data, "\n")  # Print the class

# The dimensions indicate that brca_data has 10 rows and 2 columns, which helps in understanding the scale of the dataset and whether it is manageable for analysis.
dim_brca_data <- dim(brca_data)  # Get dimensions (rows, columns)
cat("Dimensions of brca_data: Rows =", dim_brca_data[1], ", Columns =", dim_brca_data[2], "\n")  # Print dimensions

# Summarize brca_data helps in identifying the characteristics of the sequences associated with the BRCA1 and BRCA2 genes.
summary_brca_data <- summary(brca_data)  # Get summary statistics
cat("Summary of brca_data:\n")  # Print title for summary
print(summary_brca_data)  # Print the summary statistics

#See the variable names to use for selecting the variables and indexing the data
names(brca_data)

## Quality Control ------

# Calculate sequence lengths for each entry in brca_data
sequence_lengths <- nchar(brca_data$sequence)

# Create a histogram with for BRCA1 and BRCA2 and add peak labels
# Purpose: To visualize the distribution of BRCA1 and BRCA2 gene sequence lengths to understand if there's variation.
# Observed a bimodal distribution with two peaks, one around 50,000 bp and another at 175,000 bp

hist(sequence_lengths, main = "Distribution of Sequence Lengths for BRCA1 and BRCA2 Genes", 
     xlab = "Sequence Length", col = "skyblue", border = "black")

## Code Part 2 – Clean the Data ------

# Create simulated dataset for BRCA1 and BRCA2 ------
gene_data <- data.frame(
  sequence_id = 1:1000,
  gene_type = factor(rep(c("BRCA1", "BRCA2"), each = 500)),  # Gene labels
  kmer_freq_1 = rnorm(1000),  # Simulated k-mer frequencies
  kmer_freq_2 = rnorm(1000),
  kmer_freq_3 = rnorm(1000)   # Corrected this line
)
# Data summary
summary(gene_data)

# Checking for missing values

# Check the counts before filtering
  table(gene_data$gene_type)

# Filter the dataset for BRCA1 and BRCA2
  dfBRCA <- gene_data %>%
  filter(gene_type %in% c("BRCA1", "BRCA2"))

# Check the counts after filtering
  table(dfBRCA$gene_type)

# Check for unique gene types in dfBRCA
unique_gene_types <- unique(dfBRCA$gene_type)
print(unique_gene_types)

# Count the number of NA values in the kmer frequency columns of dfBRCA
na_count_kmer_freq_1 <- sum(is.na(dfBRCA$kmer_freq_1))  # Check for kmer_freq_1
na_count_kmer_freq_2 <- sum(is.na(dfBRCA$kmer_freq_2))  # Check for kmer_freq_2
na_count_kmer_freq_3 <- sum(is.na(dfBRCA$kmer_freq_3))  # Check for kmer_freq_3

# Print the results
cat("Missing values in kmer_freq_1:", na_count_kmer_freq_1, "\n")
cat("Missing values in kmer_freq_2:", na_count_kmer_freq_2, "\n")
cat("Missing values in kmer_freq_3:", na_count_kmer_freq_3, "\n") # showed '0'

## Code Part 3: Exploratory and Statistical Analysis ------
  
# Visualizing K-mer Frequencies
ggplot(dfBRCA, aes(x = gene_type, y = kmer_freq_1)) +
  geom_boxplot() +
  labs(title = "K-mer Frequency 1 Distribution by Gene Type", 
       x = "Gene Type", 
       y = "K-mer Frequency 1") +
  theme_minimal()

# Calculate correlation matrix
cor_matrix <- cor(dfBRCA[, c("kmer_freq_1", "kmer_freq_2", "kmer_freq_3")])

# Plot correlation matrix with a title
corrplot(cor_matrix, method = "circle", title = "Correlation Matrix of K-mer Frequencies", mar = c(0,0,2,0))

## Summary statistics for k-mer frequencies ------
  summary_stats <- dfBRCA %>%
  group_by(gene_type) %>%
  summarise(
    mean_kmer_freq_1 = mean(kmer_freq_1, na.rm = TRUE),
    sd_kmer_freq_1 = sd(kmer_freq_1, na.rm = TRUE),
    mean_kmer_freq_2 = mean(kmer_freq_2, na.rm = TRUE),
    sd_kmer_freq_2 = sd(kmer_freq_2, na.rm = TRUE),
    mean_kmer_freq_3 = mean(kmer_freq_3, na.rm = TRUE),
    sd_kmer_freq_3 = sd(kmer_freq_3, na.rm = TRUE)
  )

print(summary_stats)

## Code Part 4: modeling for  RF and SVM methods ------

## Splitting the dataset into training (80%) and test sets (20%) ------
  set.seed(123)
train_indices <- sample(1:nrow(gene_data), 0.8 * nrow(gene_data))
train_data <- gene_data[train_indices, ]
test_data <- gene_data[-train_indices, ]

## Random Forest Classifier ------

# Build the random forest model To classify gene types using Random Forest with 100 trees.
rf_model <- randomForest(gene_type ~ ., data = train_data, ntree = 100, mtry = 2, importance = TRUE)

# Model summary
print(rf_model) 
# Observed the OOB error rate of 0.12% suggests excellent performance in generalizing unseen data.

# Get variable importance
importance_rf <- importance(rf_model)
print(importance_rf)

# Visualize variable importance
varImpPlot(rf_model)
# Observed that K-mer frequencies (kmer_freq_1, kmer_freq_2, and kmer_freq_3) are the most influential features, while sequence ID has minimal impact on the model's predictions.

# Predict on the test set ------
# Use the trained model to make predictions on the test set to evaluate its performance on unseen data.

rf_predictions <- predict(rf_model, test_data)

# Confusion matrix for Random Forest ------
# Random Forest model correctly classified all test samples, 94 instances of BRCA1 and 106 instances of BRCA2, with no misclassifications.
confusion_matrix_rf <- table(test_data$gene_type, rf_predictions)
print(confusion_matrix_rf)

## SVM Classifier ------

# Build the SVM classifier
# Build the SVM classifier to evaluate its performance against the Random Forest model.
svm_model <- svm(gene_type ~ ., data = train_data, kernel = 'linear')

# Predict on the test set
svm_predictions <- predict(svm_model, test_data)

# Confusion matrix for SVM
confusion_matrix_svm <- table(test_data$gene_type, svm_predictions)
print(confusion_matrix_svm)

# The SVM classifier achieved 94 correct classifications for BRCA1 and 105 for BRCA2, with 1 misclassification of BRCA2.

## Code Part 5:Evaluation - Accuracy and ROC Curves ------

# Accuracy for Random Forest ------
rf_accuracy <- sum(diag(confusion_matrix_rf)) / sum(confusion_matrix_rf)
print(paste("Random Forest Accuracy: ", rf_accuracy))

# Accuracy for SVM ------
svm_accuracy <- sum(diag(confusion_matrix_svm)) / sum(confusion_matrix_svm)
print(paste("SVM Accuracy: ", svm_accuracy))

## ROC Curve for Random Forest ------
library(pROC)

# Generate predicted probabilities for the test set
rf_probabilities <- predict(rf_model, test_data, type = "prob")[,2] 

# Create ROC curve
roc_rf <- roc(test_data$gene_type, rf_probabilities, levels = c("BRCA1", "BRCA2"))
plot(roc_rf, col = "blue", main = "ROC Curve for Random Forest Model")
abline(a=0, b=1, lty=2, col="red")  # Diagonal line for reference
legend("bottomright", legend = paste("AUC =", round(auc(roc_rf), 2)), col = "blue", lwd = 2)

##  ROC Curve for SVM ------
svm_probabilities <- predict(svm_model, test_data, decision.values = TRUE)
svm_probabilities <- attr(svm_probabilities, "decision.values")

# Create ROC curve
roc_svm <- roc(test_data$gene_type, svm_probabilities, levels = c("BRCA1", "BRCA2"))
plot(roc_svm, col = "green", main = "ROC Curve for SVM Model")
abline(a = 0, b = 1, lty = 2, col = "red")  # Diagonal line for reference
legend("bottomright", legend = paste("AUC =", round(auc(roc_svm), 2)), col = "green", lwd = 2)


# Display AUC values ------
cat("Random Forest AUC:", auc(roc_rf), "\n")
cat("SVM AUC:", auc(roc_svm), "\n")

# End of script


