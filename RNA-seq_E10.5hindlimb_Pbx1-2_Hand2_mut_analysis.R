# Load required libraries
library(readxl)   # For reading Excel files
library(dplyr)    # For data manipulation
library(tidyr)    # For reshaping data
library(ggplot2)  # For plotting
library(writexl)  # To write output to excel

# === 1. Load the data ===
data <- read_excel("~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/RNA-seq_Pbx1+Hand2mutant_ctrl.xlsx", 
                   sheet = "pbx1+hand2_mutants_all")

# Ensure expression columns are numeric
data <- data %>%
  mutate(across(-1, as.numeric))

# === 2. Reshape data to long format ===
long_data <- data %>%
  pivot_longer(cols = -1,
               names_to = c("Condition", "Replicate"),
               names_sep = "_",
               values_to = "Expression") %>%
  rename(Gene = ...1)

# === 3. Summary stats ===
summary_data <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(
    mean_expression = mean(Expression, na.rm = TRUE),
    se = sd(Expression, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# === 4. Function to compute t-test p-values ===
perform_t_test <- function(data) {
  ctrl_data <- filter(data, Replicate == "ctrl")
  mutant_data <- filter(data, Replicate == "mutant")
  
  if (nrow(ctrl_data) > 1 && nrow(mutant_data) > 1) {
    if (sd(ctrl_data$Expression) > 0 && sd(mutant_data$Expression) > 0) {
      ttest_result <- t.test(ctrl_data$Expression, mutant_data$Expression)
      return(ttest_result$p.value)
    }
  }
  return(NA_real_)
}

# === 5. Compute p-values for each gene and condition ===
p_values <- long_data %>%
  group_by(Gene, Condition) %>%
  summarise(p_value = perform_t_test(cur_data()), .groups = 'drop')

# === 6. Compute log2 fold change per gene per condition ===
logfc_data <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(mean_expression = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Replicate, values_from = mean_expression) %>%
  mutate(log2FC = log2(mutant / ctrl)) %>%
  select(Gene, Condition, log2FC)

# === 7. Join p-values and logFC to summary_data ===
summary_data <- summary_data %>%
  left_join(p_values, by = c("Gene", "Condition")) %>%
  left_join(logfc_data, by = c("Gene", "Condition"))

# === 8. Prepare full expression + log2FC + p-value tables ===

# Get average expression values for ctrl and mutant
expression_summary <- long_data %>%
  group_by(Gene, Condition, Replicate) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Replicate, values_from = mean_expr) %>%
  rename(ctrl = ctrl, mutant = mutant)

# Combine with log2FC and p-values
full_results <- expression_summary %>%
  left_join(logfc_data, by = c("Gene", "Condition")) %>%
  left_join(p_values, by = c("Gene", "Condition")) %>%
  select(Gene, Condition, log2FC, p_value, ctrl, mutant)

# Subset into up/downregulated gene tables with p-value < 0.05
pbx_up <- full_results %>%
  filter(Condition == "pbx", log2FC > 0, p_value < 0.05) %>%
  arrange(desc(log2FC))

pbx_down <- full_results %>%
  filter(Condition == "pbx", log2FC < 0, p_value < 0.05) %>%
  arrange(log2FC)

hand2_up <- full_results %>%
  filter(Condition == "hand2", log2FC > 0, p_value < 0.05) %>%
  arrange(desc(log2FC))

hand2_down <- full_results %>%
  filter(Condition == "hand2", log2FC < 0, p_value < 0.05) %>%
  arrange(log2FC)


# === 9. Write everything to an Excel workbook ===
write_xlsx(
  list(
    Summary = summary_data,
    pbx_Upregulated = pbx_up,
    pbx_Downregulated = pbx_down,
    hand2_Upregulated = hand2_up,
    hand2_Downregulated = hand2_down
  ),
  path = "~/Documents/postdoc/bioinformatics/data/RNA-seq/RNA-seq_E10.5hindlimb_Pbx1-2_Hand2_mut/analysis/RNA-seq_Pbx1+Hand2_mutVSctrl_analysis.xlsx"
)



# Filter for the gene of interest
gene_of_interest <- "Hand1"  # Replace with your actual gene name
gene_data <- summary_data %>%
  filter(Gene == gene_of_interest)

# Check if gene_data is empty
if (nrow(gene_data) == 0) {
  stop("No data available for the specified gene of interest.")
}

# Create a new variable to combine Replicate and Condition
gene_data <- gene_data %>%
  mutate(Group = interaction(Replicate, Condition))

# Create the bar plot for the specific gene
ggplot(gene_data, aes(x = Group, y = mean_expression, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +  # Black border around bars
  geom_errorbar(aes(ymin = mean_expression - se, ymax = mean_expression + se), 
                width = 0.2, position = position_dodge(0.9)) +
  labs(title = paste("Expression of", gene_of_interest, "in E10.5 hindlimb"),
       x = "Biological context",
       y = "Bulk RNA-seq mean expression [AU]") +
  theme_minimal() +
  scale_fill_manual(values = c("hand2" = "lightblue", "pbx" = "purple"), 
                    labels = c("hand2" = "Hand2", "pbx" = "Pbx")) +  # Custom colors and labels
  geom_text(data = gene_data %>% filter(Replicate == "mutant"),
            aes(label = ifelse(!is.na(p_value), paste("p =", round(p_value, 5)), "")), 
            position = position_dodge(0.5), vjust = -2, size = 3) +  # Display p-values only for mutant condition
  theme(legend.title = element_text(size = 10),  # Customize legend title size
        legend.text = element_text(size = 10),   # Customize legend text size
        legend.position = "right") +              # Position the legend on the right
  guides(fill = guide_legend(title = "Condition"))  # Add legend title

