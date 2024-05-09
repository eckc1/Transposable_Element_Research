library(ggplot2)
library(reshape2)


#########################################################################
######################## Revised Heatmap of Single_End_Data ########################
dataset <- read.csv("/Users/connoreck/Desktop/DIGIT/everything/Single_combined_genes.csv")
head(dataset)

# Define the breaks for the bins as per your requirements
num_hits_bins <- c(1, 3, 120, 100000)
percent_diff_bins <- c(-Inf, -1, 0, 1, 3, 7, 15, Inf)

# Create bins for NumHitsForGenome and PercentDifferenceToNextHit using cut()
dataset$NumHitsForGenome_Bin <- cut(dataset$NumHitsForGenome, breaks = num_hits_bins, include.lowest = TRUE)
dataset$PercentDifferenceToNextHit_Bin <- cut(dataset$PercentDifferenceToNextHit, breaks = percent_diff_bins, include.lowest = TRUE)

# Create a table of counts for each combination of bins
heatmap_data <- dataset %>%
  group_by(NumHitsForGenome_Bin, PercentDifferenceToNextHit_Bin) %>%
  summarise(Count = n(), .groups = 'drop')  # Drop the grouping for later steps

# Filter out the [-Inf,-1] bin from PercentDifferenceToNextHit_Bin before creating the heatmap
heatmap_data <- heatmap_data %>%
  filter(PercentDifferenceToNextHit_Bin != "[-Inf,-1]")
#filter(NumHitsForGenome_Bin != "[0.000000001,1]")

# Ensure the acast function is available by loading the reshape2 package
# If not already loaded, uncomment the line below to install and load reshape2
# install.packages("reshape2"); library(reshape2)

# Spread the data for the heatmapA
heatmap_matrix <- acast(heatmap_data, PercentDifferenceToNextHit_Bin ~ NumHitsForGenome_Bin, value.var = "Count")

# Replace NAs with 0s for heatmap
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Ensure the data is in matrix form for plotting
if (!is.matrix(heatmap_matrix)) {
  heatmap_matrix <- as.matrix(heatmap_matrix)
}

# Melt the matrix for use in ggplot2
heatmap_melted <- melt(heatmap_matrix, varnames = c("PercentDifferenceToNextHit_Bin", "NumHitsForGenome_Bin"))

# Custom Labels for NumHitsForGenome Bins
num_hits_labels <- c("2", "3-120", ">120")

# Custom Labels for PercentDifferenceToNextHit Bins
percent_diff_labels <- c("-1 to 0%", "0 to 1%", "1 to 3%", "4 to 7%", "7 to 15%", ">15%")

# Plot the heatmap with custom labels
ggplot(heatmap_melted, aes(x = NumHitsForGenome_Bin, y = PercentDifferenceToNextHit_Bin, fill = value)) +
  geom_tile(color = "white") +  # Add the tiles
  geom_text(aes(label = value), color = "black", size = 3) +  # Add text labels
  scale_fill_gradient(low = "blue", high = "red") +  # Color gradient
  scale_x_discrete(labels = num_hits_labels, name = "NumHitsForGenome Bins") +  # Custom X-axis labels
  scale_y_discrete(labels = percent_diff_labels, name = "PercentDifferenceToNextHit Bins") +  # Custom Y-axis labels
  labs(title = "Gene Distribution Single Insertion", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, vjust = 1))


#########################################################################
library(ggplot2)
library(reshape2)

# Manually input the matrix with the provided proportions

heatmap_matrix <- matrix(
  c(
    0.039, 0.063, 0.095, 
    0.018, 0.030, 0.046, 
    0.068, 0.062, 0.063, 
    0.103, 0.068, 0.023,
    0.114, 0.053, 0.010,
    0.114, 0.028, 0.003
  ), 
  nrow = 6, 
  byrow = TRUE
)

# Define row and column names
rownames(heatmap_matrix) <- c("-1 to 0%", "0 to 1%", "1 to 3%", "4 to 7%", "7 to 15%", ">15%")
colnames(heatmap_matrix) <- c("2","3-120", ">120")

# Convert the matrix to a long format
heatmap_melted <- melt(heatmap_matrix, varnames = c("PercentDifferenceToNextHit_Bin", "NumHitsForGenome_Bin"))

# Plot the heatmap
ggplot(heatmap_melted, aes(x = NumHitsForGenome_Bin, y = PercentDifferenceToNextHit_Bin, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = 0.5, color = "black", size = 3) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap of Proportions", x = "Num Hits For Genome Bin", y = "Percent Difference To Next Hit Bin", fill = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##################################
#Adjusted_Bins3_less
heatmap_matrix <- matrix(
  c(
    0.039, 0.063, 0.095, 
    0.018, 0.030, 0.046, 
    0.068, 0.062, 0.063, 
    0.103, 0.068, 0.023,
    0.114, 0.053, 0.010,
    0.114, 0.028, 0.003
  ), 
  nrow = 6, 
  byrow = TRUE
)


# Adjusted_Bins2_more

heatmap_matrix <- matrix(
  c(
    0.039, 0.021, 0.028,0.013,0.095, 
    0.018, 0.013, 0.011,0.006,0.046, 
    0.068, 0.023, 0.026,0.013,0.063,
    0.103, 0.024, 0.027,0.017,0.023,
    0.114, 0.018, 0.023,0.012,0.010,
    0.114, 0.009, 0.011,0.008,0.003
  ), 
  nrow = 6, 
  byrow = TRUE
)



# Adjusted_Bins3_med

heatmap_matrix <- matrix(
  c(
    0.039, 0.047,0.013,0.095, 
    0.018, 0.024,0.006,0.046, 
    0.068, 0.049,0.013,0.063,
    0.103, 0.050,0.017,0.023,
    0.114, 0.042,0.012,0.010,
    0.114, 0.020,0.008,0.003
  ), 
  nrow = 6, 
  byrow = TRUE
)





############# Heatmap of eValues #############
library(ggplot2)
library(reshape2)

# Load datasets
eValue_above_50 <- read.csv("eValue_above_50_single.csv")
eValue_below_50 <- read.csv("eValue_below_50_single.csv")
#eValue_above_50 <- read.csv("eValue_above_50_double.csv")
#eValue_below_50 <- read.csv("eValue_below_50_double.csv")
# Define bins
bins <- seq(0, -2000, by = -10)

# Bin eValue for each dataset
eValue_above_50$Bin <- cut(eValue_above_50$eValue, breaks = bins, include.lowest = TRUE)
eValue_below_50$Bin <- cut(eValue_below_50$eValue, breaks = bins, include.lowest = TRUE)

# Aggregate Query frequency per bin
freq_above_50 <- aggregate(eValue_above_50$Query, by = list(Bin = eValue_above_50$Bin), FUN = length)
freq_below_50 <- aggregate(eValue_below_50$Query, by = list(Bin = eValue_below_50$Bin), FUN = length)

# Rename columns
colnames(freq_above_50) <- c("Bin", "Frequency_above_50")
colnames(freq_below_50) <- c("Bin", "Frequency_below_50")

# Merge datasets
merged_freq <- merge(freq_above_50, freq_below_50, by = "Bin", all = TRUE)
merged_freq[is.na(merged_freq)] <- 0  # Replace NA with 0

# Melt for plotting
melted_data <- melt(merged_freq, id.vars = "Bin")

# Convert Bin to factor for ordering
melted_data$Bin <- factor(melted_data$Bin, levels = rev(levels(melted_data$Bin)))


# Plot with numbers inside the boxes
ggplot(melted_data, aes(x = variable, y = Bin, fill = value)) + 
  geom_tile() + # Draws the heatmap tiles
  geom_text(aes(label = value), color = "white", size = 3) + # Adds text labels inside the tiles
  scale_fill_gradient(low = "blue", high = "red") + # Color gradient
  labs(title = "Heatmap of eValue Frequency Single Insertion", x = "", y = "eValue Bin", fill = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improves x-axis labels readability

####################################################################################################

######################## Revised Heatmap of Single_End eValues vs NumHitsForGenome ########################
# Load required libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# Read the dataset
dataset <- read.csv("/Users/connoreck/Desktop/DIGIT/everything/eValue_logged_numhits.csv")

# Assuming `NumHitsForGenome` and `eValue` are correct column names in your dataset
# Define the breaks for the bins
num_hits_bins <- c(1, 3, 120, 100000)
evalue_bins <- c(0, -10, -20, -30, -40, -50, -60, -70, -80, -90, -100, -110, -120, -130)

# Create bins 
dataset$NumHitsForGenome_Bin <- cut(dataset$NumHitsForGenome, breaks = num_hits_bins, include.lowest = TRUE)
dataset$eValue_Bin <- cut(dataset$eValue, breaks = evalue_bins, include.lowest = TRUE)

# Create a table of counts for each combination of bins
heatmap_data <- dataset %>%
  group_by(NumHitsForGenome_Bin, eValue_Bin) %>%
  summarise(Count = n(), .groups = 'drop')

# Filter out unnecessary bins before creating the heatmap if needed
# Example: filter(PercentDifferenceToNextHit_Bin != "[-Inf,-1]")

# Spread the data for the heatmap
heatmap_matrix <- acast(heatmap_data, eValue_Bin ~ NumHitsForGenome_Bin, value.var = "Count")

# Replace NAs with 0s for heatmap
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Ensure the data is in matrix form for plotting
if (!is.matrix(heatmap_matrix)) {
  heatmap_matrix <- as.matrix(heatmap_matrix)
}

# Melt the matrix for use in ggplot2
heatmap_melted <- melt(heatmap_matrix)

# Custom Labels for bins (modify as per your actual bin labels)
num_hits_labels <- c("1", "3-120", ">120")
evalue_labels <- c("0", "-10", "-20", "-30", "-40", "-50", "-60", "-70", "-80", "-90", "-100", "-110", "-120", "-130")

# Plot the heatmap with custom labels
ggplot(heatmap_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels = num_hits_labels, name = "NumHitsForGenome Bins") +
  scale_y_discrete(labels = evalue_labels, name = "eValue Bins") +
  labs(title = "eValue Distribution Single Insertion", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, vjust = 1))






