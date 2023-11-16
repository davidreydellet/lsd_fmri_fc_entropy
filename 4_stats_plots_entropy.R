# PLOT ENTROPY METRICS

# Load necessary library
library(dplyr)
library(ggplot2)
library(svglite)
library(tidyverse)
library(R.matlab)

# Read csv file
data <- read.csv("/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long.csv")
# Select specific columns and remove duplicate rows
data <- data %>%
  select(time_series_complexity, geodesic_entropy, meta_state_complexity, normalized_modularity, participant, session, drugsbefore_binary, sex, age, bmi, dataset) %>%
  distinct()

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Set up colors for each session
session_colors <- c("ses_plc" = "#808080", "ses_lsd" = "#FF69B4") 



# Create the LZ plot
LZplot <- ggplot(data, aes(x = session, y = time_series_complexity, color = session, fill = session)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
  geom_point(position = position_dodge2(width = 0.5), size = 2) + 
  scale_color_manual(values = session_colors) +
  scale_fill_manual(values = session_colors) +
  scale_x_discrete(limits = c("ses_plc", "ses_lsd"), labels = c("Placebo", "LSD")) +
  labs(
    x = "Session",
    y = "Temporal BOLD Complexity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(LZplot)


# Create the path-length distribution plot
PLplot <- ggplot(data, aes(x = session, y = geodesic_entropy, color = session, fill = session)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
  geom_point(position = position_dodge2(width = 0.5), size = 2) + 
  scale_color_manual(values = session_colors) +
  scale_fill_manual(values = session_colors) +
  scale_x_discrete(limits = c("ses_plc", "ses_lsd"), labels = c("Placebo", "LSD")) +
  labs(
    x = "Session",
    y = "Path-length distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(PLplot)


# Create the meta-state complexity plot
MSCplot <- ggplot(data, aes(x = session, y = meta_state_complexity, color = session, fill = session)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
  geom_point(position = position_dodge2(width = 0.5), size = 2) + 
  scale_color_manual(values = session_colors) +
  scale_fill_manual(values = session_colors) +
  scale_x_discrete(limits = c("ses_plc", "ses_lsd"), labels = c("Placebo", "LSD")) +
  labs(
    x = "Session",
    y = "Meta-state complexity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(MSCplot)



# Create the normalized modularity plot
MODplot <- ggplot(data, aes(x = session, y = normalized_modularity, color = session, fill = session)) +
  geom_violin(alpha = 0.5, scale = "width", trim = FALSE) +
  geom_point(position = position_dodge2(width = 0.5), size = 2) + 
  scale_color_manual(values = session_colors) +
  scale_fill_manual(values = session_colors) +
  scale_x_discrete(limits = c("ses_plc", "ses_lsd"), labels = c("Placebo", "LSD")) +
  labs(
    x = "Session",
    y = "Normalized modularity"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Display the plot
print(MODplot)



# Define working space
directory <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/"

# Save the plots as SVG
ggsave(paste0(directory, "/PL_plot.svg"), plot = PLplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/LZ_plot.svg"), plot = LZplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/MSC_plot.svg"), plot = MSCplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/MOD_plot.svg"), plot = MODplot, width = 4, height = 4, units = "in")

# Save the plots as PNG
ggsave(paste0(directory, "/PL_plot.png"), plot = PLplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/LZ_plot.png"), plot = LZplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/MSC_plot.png"), plot = MSCplot, width = 4, height = 4, units = "in")
ggsave(paste0(directory, "/MOD_plot.png"), plot = MODplot, width = 4, height = 4, units = "in")


# Function to compute Cohen's d
calculate_cohens_d <- function(df, column_name, group1, group2) {
  
  # Filter data by group and remove NAs
  group1_data <- na.omit(df[df$session == group1, ][[column_name]])
  group2_data <- na.omit(df[df$session == group2, ][[column_name]])
  
  # Calculate means and standard deviations for each group
  mean_group1 <- mean(group1_data)
  mean_group2 <- mean(group2_data)
  
  sd_group1 <- sd(group1_data)
  sd_group2 <- sd(group2_data)
  
  # Calculate sample sizes for each group
  n_group1 <- length(group1_data)
  n_group2 <- length(group2_data)
  
  # Calculate pooled standard deviation
  s_pooled <- sqrt(((n_group1-1) * sd_group1^2 + (n_group2-1) * sd_group2^2) / (n_group1 + n_group2 - 2))
  
  # Calculate Cohen's d
  cohen_d <- (mean_group1 - mean_group2) / s_pooled
  
  return(cohen_d)
}


# Fit LMM + compute Cohen's d
# Load the required packages
library(lme4)
library(lmerTest)
library(emmeans)
library(LMMstar)# Required to fit linear mixed models with Brice's model 


# Fit the model for LZ
modelLZ <- lmer(time_series_complexity ~ session + drugsbefore_binary + age + sex + bmi + dataset + (1|participant), data = data)

# Store statistics in a table_entropy
table_entropy <- data.frame(
  entropy = "time_series_complexity",
  effect_size = calculate_cohens_d(data, "time_series_complexity", "ses_lsd", "ses_plc"),
  estimate = coef(summary(modelLZ))['sessionses_lsd','Estimate'],
  lower_confint = confint(modelLZ)['sessionses_lsd','2.5 %'],
  upper_confint = confint(modelLZ)['sessionses_lsd','97.5 %'],
  pvalue = coef(summary(modelLZ))['sessionses_lsd','Pr(>|t|)']
)


# Fit the model for PL
modelPL <- lmer(geodesic_entropy ~ session + drugsbefore_binary + age + sex + bmi + dataset + (1|participant), data = data)
# Store statistics in a table_entropy
table_entropy <- rbind(table_entropy, 
                       data.frame(
                         entropy = "path_length_distribution",
                         effect_size = calculate_cohens_d(data, "geodesic_entropy", "ses_lsd", "ses_plc"),
                         estimate = coef(summary(modelPL))['sessionses_lsd','Estimate'],
                         lower_confint = confint(modelPL)['sessionses_lsd','2.5 %'],
                         upper_confint = confint(modelPL)['sessionses_lsd','97.5 %'],
                         pvalue = coef(summary(modelPL))['sessionses_lsd','Pr(>|t|)']
                       ))


# Fit the model for MSC
modelMSC <- lmer(meta_state_complexity ~ session + drugsbefore_binary + age + sex + bmi + dataset + (1|participant), data = data)
# Store statistics in a table_entropy
table_entropy <- rbind(table_entropy, 
                       data.frame(
                         entropy = "meta_state_complexity",
                         effect_size = calculate_cohens_d(data, "meta_state_complexity", "ses_lsd", "ses_plc"),
                         estimate = coef(summary(modelMSC))['sessionses_lsd','Estimate'],
                         lower_confint = confint(modelMSC)['sessionses_lsd','2.5 %'],
                         upper_confint = confint(modelMSC)['sessionses_lsd','97.5 %'],
                         pvalue = coef(summary(modelMSC))['sessionses_lsd','Pr(>|t|)']
                       ))



# Fit the model for MOD
modelMOD <- lmer(normalized_modularity ~ session + drugsbefore_binary + age + sex + bmi + dataset + (1|participant), data = data)
# Store statistics in a table_entropy
table_entropy <- rbind(table_entropy, 
                       data.frame(
                         entropy = "normalized_modularity",
                         effect_size = calculate_cohens_d(data, "normalized_modularity", "ses_lsd", "ses_plc"),
                         estimate = coef(summary(modelMOD))['sessionses_lsd','Estimate'],
                         lower_confint = confint(modelMOD)['sessionses_lsd','2.5 %'],
                         upper_confint = confint(modelMOD)['sessionses_lsd','97.5 %'],
                         pvalue = coef(summary(modelMOD))['sessionses_lsd','Pr(>|t|)']
                       ))






# MSSE sample_entropy

# Read csv file
data <- read.csv("/student/davidreydellet/lsd-basel/data/derivative/analysis/data_long_withMSSE.csv")

# Create a new column to loop correctly with mlmm
data <- data %>%
  mutate(scale_network_shen_MSSE = paste0("scale", scale_MSSE, "_", network_shen_MSSE))

# Select specific columns and remove duplicate rows
data <- data %>%
  select(avg_MSSE, scale_network_shen_MSSE, participant, session, drugsbefore_binary, sex, age, bmi, dataset) %>%
  distinct()

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))


# Fit LMM using mlmm
ls.model <- mlmm(avg_MSSE ~ session + drugsbefore_binary + age + sex + bmi + dataset, by = "scale_network_shen_MSSE",
                 repetition = ~ session|participant, structure = "UN", data = data,
                 effects = "sessionses_lsd=0", robust = FALSE)
summary(ls.model, method = "single-step")
plot(ls.model, method = "none")
table <- model.tables(ls.model, columns = add("parameter"), method = "single-step")
table_uncorrected <- model.tables(ls.model, columns = add("parameter"), method = "none")

# Function to compute Cohen's d based on groups
calculate_grouped_cohens_d <- function(df, column_name, group_column, group1, group2) {
  unique_values <- unique(df[[group_column]])
  results <- list()
  
  for (value in unique_values) {
    # Filter dataframe by the current unique value of the grouping column
    sub_df <- df[df[[group_column]] == value, ]
    
    # Compute Cohen's d for this subgroup
    cohen_d_value <- calculate_cohens_d(sub_df, column_name, group1, group2)
    results[[value]] <- cohen_d_value
  }
  
  return(results)
}

# Compute Cohen's d
grouped_effect_sizes <- 
print(grouped_effect_sizes)

# Create table_MSSE
table_MSSE <- table %>%
  left_join(table_uncorrected, by = "by") %>%
  transmute(
    entropy = paste0("MSSE_", by),
    effect_size = calculate_grouped_cohens_d(data, "avg_MSSE", "scale_network_shen_MSSE", "ses_lsd", "ses_plc"),
    estimate = estimate.x, # assuming "estimate" from "table"
    lower_confint = lower.y, # "lower" from "table_uncorrected"
    upper_confint = upper.y, # "upper" from "table_uncorrected"
    lower_confint_singlestep = lower.x, # "lower" from "table"
    upper_confint_singlestep = upper.x, # "upper" from "table"
    pvalue = p.value.y, # "p.value" from "table_uncorrected"
    pvalue_singlestep = p.value.x  # "p.value" from "table"
  )

# Re-organizing columns for table_MSSE
table_MSSE <- table_MSSE %>%
  select(entropy, effect_size, estimate, lower_confint, upper_confint, pvalue, lower_confint_singlestep, upper_confint_singlestep, pvalue_singlestep, everything())

# Convert the effect_size column from list to vector
table_MSSE$effect_size <- unlist(table_MSSE$effect_size)



# Forest plot for MSSE

table_MSSE2 <- table_MSSE %>% 
  separate(entropy, into = c("MSSE", "Scale", "ROI"), sep = "_", remove = TRUE) %>% 
  select(-MSSE) %>% 
  mutate_at(vars(Scale, ROI), as.factor) %>% 
  mutate_at(vars(estimate, lower_confint, upper_confint, pvalue), as.numeric)

table_MSSE2$ROI_num <- as.numeric(gsub("ROI", "", table_MSSE2$ROI))
table_MSSE2$ROI <- factor(table_MSSE2$ROI, levels = na.omit(unique(table_MSSE2$ROI)[order(table_MSSE2$ROI_num)]))

filter(table_MSSE2, pvalue_singlestep<0.05)

dose.labs <- c("Scale 1", "Scale 2", "Scale 3", "Scale 4", "Scale 5")
names(dose.labs) <- c("scale1", "scale2", "scale3", "scale4", "scale5")
ROI_names = gsub("ROI", "", table_MSSE2$ROI)
  
# Define a lookup table for the full network names
lookup_table <- data.frame(ROI_abbrev = c("Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                           ROI_full = c("Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))

table_MSSE2 <- merge(table_MSSE2, lookup_table, by.x = "ROI", by.y = "ROI_abbrev")
ROI_names = gsub("ROI", "", table_MSSE2$ROI_full)

# Add a star column
table_MSSE2 <- table_MSSE2 %>%
  mutate(star = case_when(
    pvalue_singlestep < 0.001 ~ "***",
    pvalue_singlestep < 0.01  ~ "**",
    pvalue_singlestep < 0.05  ~ "*",
    TRUE ~ NA_character_
  ))

# Create the forest plot
MSSEplot <- ggplot(table_MSSE2, aes(x = ROI_full, y = estimate, ymin = lower_confint_singlestep, ymax = upper_confint_singlestep, colour = ROI_full)) +
  geom_errorbar(width = 0.05) +
  geom_point() +
  theme_bw() +
  geom_text(aes(y = upper_confint_singlestep + 0.005, label = star), color = "black", size = 2, inherit.aes = TRUE, na.rm = TRUE) +  
  geom_abline(slope = 0, intercept = 0, color = "black") +
  #scale_y_continuous(limits = c(-0.004, 0.0021), expand = c(0,0), breaks = c(-0.004, -0.002, 0, 0.002)) +
  labs(x = "Network", y = "Estimate \u00B1 95% CI") +
  facet_wrap(~Scale, ncol = 3, labeller = labeller(Scale = dose.labs)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size = 8), legend.title = element_blank()) 
  #+scale_color_discrete(labels = ROI_full) +  guides(color = guide_legend(ncol = 1))

# Save MSSEplot
directory <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/"
ggsave(paste0(directory, "/MSSE_plot.png"), plot = MSSEplot, width = 8, height = 6, dpi = 300)
ggsave(paste0(directory, "/MSSE_plot.svg"), plot = MSSEplot, width = 8, height = 6, units = "in")






# DCC
# Load the necessary packages
library(lme4)   # Required to fit linear mixed models
library(LMMstar)# Required to fit linear mixed models with Brice's model 
library(ggplot2)
library(reshape2)
library(dplyr)

# Load csv file
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data <- read.csv("data_long_withDCC.csv")

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$network <- factor(data$network)
data$time <- factor(data$time)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$drugsbefore <- factor(data$drugsbefore) # nb of drug sessions before the current session, can be either 0, 1, 2 or 3
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Select specific columns and remove duplicate rows
data <- data %>%
  select(DCC, network_shen, participant, session, drugsbefore_binary, sex, age, bmi, dataset) %>%
  distinct()



# Fit LMM
ls.model <- mlmm(DCC ~ session + drugsbefore_binary + age + sex + bmi + dataset, by = "network_shen",
                   repetition = ~ session|participant, structure = "UN", data = data,
                   effects = "sessionses_lsd=0", robust = FALSE)
summary(ls.model, method = "single-step")
plot(ls.model, method = "none")
table <- model.tables(ls.model, columns = add("parameter"), method = "single-step")
table_uncorrected <- model.tables(ls.model, columns = add("parameter"), method = "none")

# Extract unique network names
networks <- unique(unlist(strsplit(as.character(table$by), "_")))

# Create empty matrices for parameter estimate and pvalue
parameter_matrix <- matrix(NA, nrow = length(networks), ncol = length(networks),
                           dimnames = list(networks, networks))
pvalue_matrix <- matrix(NA, nrow = length(networks), ncol = length(networks),
                        dimnames = list(networks, networks))

# Loop through each pair of networks
for (i in 1:length(networks)) {
  for (j in i:length(networks)) {
    # Extract the pair of networks
    network_pair <- paste(networks[i], networks[j], sep = "_")
    
    # Filter the rows of the table corresponding to the network pair
    network_pair_rows <- table[table$by == network_pair, ]
    
    # Compute the parameters and p-values
    parameter_value <- mean(network_pair_rows$estimate, na.rm = TRUE)
    pvalue_value <- mean(network_pair_rows$p.value, na.rm = TRUE)
    
    # Populate the corresponding cells in the matrices
    parameter_matrix[networks[i], networks[j]] <- parameter_value
    parameter_matrix[networks[j], networks[i]] <- parameter_value
    pvalue_matrix[networks[i], networks[j]] <- pvalue_value
    pvalue_matrix[networks[j], networks[i]] <- pvalue_value
  }
}

# Create a dataframe with the parameters and pvalues
df <- melt(parameter_matrix)
df$pvalue <- melt(pvalue_matrix)$value

# Add stars based on pvalues
df$stars <- ifelse(df$pvalue < 0.001, "***",
                   ifelse(df$pvalue < 0.01, "**",
                          ifelse(df$pvalue < 0.05, "*", "")))

# Define a lookup table for the full network names
lookup_table <- data.frame(abbrev = c("Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                           full = c("Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))

# Use left_join() to join the lookup_table to the original data frame
# by matching the "abbrev" column in lookup_table with "Var1" and "Var2" columns in df
df <- df %>%
  left_join(lookup_table, by = c("Var1" = "abbrev")) %>%
  left_join(lookup_table, by = c("Var2" = "abbrev"))

# Rename the new columns to Var1_full and Var2_full
colnames(df)[6:7] <- c("Network1_full", "Network2_full")

# Define the order of the x and y axes
order <- c("Visual 1", "Visual 2", "Visual Association", "Subcortical Cerebellum", 
           "Motor", "Medial Frontal", 
           "Frontoparietal", "Default Mode")

# Subset the data to only keep values below or on the diagonal
df_sub <- subset(df, as.numeric(factor(Network1_full, levels = order)) >= as.numeric(factor(Network2_full, levels = order)))

# Plot beta DCC
DCCplot <- 
  ggplot(df_sub, aes(x = factor(df_sub$Network1_full, levels = order), 
                     y = factor(df_sub$Network2_full, levels = order), 
                     fill = value)) +
  geom_tile(colour = "black", size = 0.1) + # Change the borders around individual tiles
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "β", 
                       guide = guide_colorbar(barwidth = 0.9, barheight = 5)) +
  geom_text(aes(label = stars), color = "black", size = 4) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
        axis.text.y = element_text(hjust = 1, size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.1, 0.9), # Set the position of the legend
        legend.justification = c(0, 0.8), # Set the justification of the legend
        legend.box.just = "left",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) + 
  coord_fixed()  # Tiles are squared

# Save
directory <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/"
ggsave(paste0(directory, "/DCC_plot.png"), plot = DCCplot, width = 6, height = 6, dpi = 300)
ggsave(paste0(directory, "/DCC_plot.svg"), plot = DCCplot, width = 6, height = 6, units = "in")


# DCC Cohen's d
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
# Load csv file
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data <- read.csv("data_long_withDCC.csv")

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$network_shen <- factor(data$network_shen)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Select specific columns and remove duplicate rows
data <- data %>%
  select(DCC, network_shen, participant, session, drugsbefore_binary, sex, age, bmi, dataset) %>%
  distinct()

# Split network column into two separate columns
data_subset <- data %>%
  separate(network_shen, into = c("network_1", "network_2"), sep = "_")

# Calculate sample size, mean, and standard deviation for each pair of networks and session
network_summary <- dplyr::group_by(data_subset, network_1, network_2, session) %>%
  dplyr::summarize(sample_size = sum(!is.na(DCC)), 
                   mean_DCC = mean(DCC, na.rm = TRUE),
                   sd_DCC = sd(DCC, na.rm = TRUE)) %>%
  dplyr::ungroup()


# Calculate effect size for each pair of networks
effect_size <- network_summary %>%
  dplyr::group_by(network_1, network_2) %>%
  dplyr::summarize(cohen_d = (dplyr::first(mean_DCC[session == "ses_lsd"]) - dplyr::first(mean_DCC[session == "ses_plc"])) / 
                     sqrt(((dplyr::first(sample_size[session == "ses_lsd"]) - 1) * (dplyr::first(sd_DCC[session == "ses_lsd"]))^2 +
                             (dplyr::first(sample_size[session == "ses_plc"]) - 1) * (dplyr::first(sd_DCC[session == "ses_plc"]))^2) /
                            (dplyr::first(sample_size[session == "ses_lsd"]) + dplyr::first(sample_size[session == "ses_plc"]) - 2))) # pooled variance

effect_size_temp = effect_size

# Define a lookup table for the full network names
lookup_table <- data.frame(abbrev = c("Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                           full = c("Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))

# Rename networks in effect_size using lookup_table
effect_size$Network1_full <- lookup_table$full[match(effect_size$network_1, lookup_table$abbrev)]
effect_size$Network2_full <- lookup_table$full[match(effect_size$network_2, lookup_table$abbrev)]

# Remove the original abbreviated network name columns
effect_size$network_1 <- NULL
effect_size$network_2 <- NULL

# Define the order of the x and y axes
order <- c("Visual 1", "Visual 2", "Visual Association", "Subcortical Cerebellum", 
           "Motor", "Medial Frontal", 
           "Frontoparietal", "Default Mode")

# Create a new column in effect_size to indicate whether the network names should be reversed
effect_size$reverse_order <- ifelse(match(effect_size$Network1_full, order) < match(effect_size$Network2_full, order), TRUE,
                                    ifelse(match(effect_size$Network1_full, order) > match(effect_size$Network2_full, order), FALSE, NA))
effect_size$reverse_order[which(effect_size$Network1_full == effect_size$Network2_full)] <- FALSE

# Intermediate step
effect_size$Network1_full_reversed <- ifelse(effect_size$reverse_order, effect_size$Network2_full, effect_size$Network1_full)
effect_size$Network2_full_reversed <- ifelse(effect_size$reverse_order, effect_size$Network1_full, effect_size$Network2_full)

# Replace the original columns
effect_size$Network1_full <- effect_size$Network1_full_reversed
effect_size$Network2_full <- effect_size$Network2_full_reversed

# Remove the 'reverse_order' and 'Network1_full_reversed' columns
effect_size$reverse_order <- NULL
effect_size$Network1_full_reversed <- NULL
effect_size$Network2_full_reversed <- NULL

# Add significance stars from df_sub
effect_size <- merge(effect_size, df_sub[c("Network1_full", "Network2_full", "stars")], 
                     by = c("Network1_full", "Network2_full"), all.x = TRUE)


# Plot DCC_cohenplot
DCC_cohenplot <- 
  ggplot(effect_size, aes(x = factor(Network1_full, levels = order), 
                          y = factor(Network2_full, levels = order), 
                          fill = cohen_d)) +
  geom_tile(colour = "black", size = 0.1) + # Change the borders around individual tiles
  geom_text(aes(label = stars), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Cohen's d",
                       breaks = c(0, 0.6, 1.2),
                       limits = c(0, 1.2),
                       guide = guide_colorbar(barwidth = 0.9, barheight = 5)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 15),
        axis.text.y = element_text(hjust = 1, size = 15),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 13),
        legend.position = c(0.1, 0.9), # Set the position of the legend
        legend.justification = c(0, 0.8), # Set the justification of the legend
        legend.box.just = "left") + # Set the justification of the legend box)
  coord_fixed()  # Tiles are squared


# Save
directory <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/"
ggsave(paste0(directory, "/DCC_cohenplot.png"), plot = DCC_cohenplot, width = 6, height = 6, dpi = 300)
ggsave(paste0(directory, "/DCC_cohenplot.svg"), plot = DCC_cohenplot, width = 6, height = 6, units = "in")




# Add Cohen's d to table
# Create a 'by' column in effect_size_temp using network_1 and network_2
effect_size_temp <- effect_size_temp %>%
  mutate(by = paste(network_1, network_2, sep = "_"))
# Left join table with effect_size_temp based on 'by' column
table <- left_join(table, effect_size_temp[, c("by", "cohen_d")], by = "by")


# Create table_DCC
table_DCC <- table %>%
  select(by, cohen_d, estimate, lower, upper, p.value) %>%
  rename(entropy = by,
         effect_size = cohen_d,
         lower_confint_singlestep = lower,
         upper_confint_singlestep = upper,
         pvalue_singlestep = p.value)
# Add uncorrected pvalues
table_DCC <- table_DCC %>%
  left_join(select(table_uncorrected, by, p.value, lower, upper), by = c("entropy" = "by")) %>%
  rename(
    pvalue = p.value,
    lower_confint = lower,
    upper_confint = upper
  ) %>%
  select(-pvalue_singlestep, pvalue, everything())

# Re-organizing columns for table_DCC
table_DCC <- table_DCC %>%
  select(entropy, effect_size, estimate, lower_confint, upper_confint, pvalue, lower_confint_singlestep, upper_confint_singlestep, pvalue_singlestep, everything())

# Rename entropy column add a DCC mention
table_DCC <- table_DCC %>%
  mutate(entropy = paste0("DCC_", entropy))







# Save table_entropy & table_MSSE & table_DCC
directory <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/"
# table_entropy
filename_entropy <- "table_entropy.csv"
full_path_entropy <- paste0(directory, filename_entropy)
write.csv(table_entropy, file = full_path_entropy, row.names = TRUE)
# table_MSSE
filename_MSSE <- "table_MSSE.csv"
full_path_MSSE <- paste0(directory, filename_MSSE)
write.csv(table_MSSE, file = full_path_MSSE, row.names = TRUE)
# table_DCC
filename_DCC <- "table_DCC.csv"
full_path_DCC <- paste0(directory, filename_DCC)
write.csv(table_DCC, file = full_path_DCC, row.names = TRUE)



# Done!








# # Load the necessary libraries
# library(dplyr)
# library(stringr)
# 
# # Read the CSV file into a dataframe
# path <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/table_DCC.csv"
# df <- read.csv(path, stringsAsFactors = FALSE)
# 
# # Manipulate the `entropy` column
# df$entropy <- df$entropy %>%
#   str_replace("DCC_", "") %>%
#   str_replace_all("_", " - ") %>%
#   paste0("\t\t", .)  # prepend two tab characters
# 
# # Save the modified dataframe back to a CSV file
# write.csv(df, "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/table_DCC_temp.csv", row.names = FALSE)
# 
# 
# 
# 
# # Read the CSV file into a dataframe
# path <- "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/table_MSSE.csv"
# df <- read.csv(path, stringsAsFactors = FALSE)
# 
# # Manipulate the dataframe to extract scale and network name
# df <- df %>%
#   # Extract the scale using string manipulation
#   mutate(scale = as.numeric(str_extract(entropy, "(?<=MSSE_scale)\\d")),
#          # Extract everything after the second underscore for the network name
#          `Network name` = paste0("\t\t", str_extract(entropy, "(?<=MSSE_scale\\d_)\\w+$"))) %>%
#   # Reorder the columns so that 'scale' and 'Network name' are right after 'entropy'
#   select(entropy, scale, `Network name`, everything())
# 
# # Save the modified dataframe back to a CSV file
# write.csv(df, "/student/davidreydellet/lsd-basel/data/derivative/analysis/entropy/table_MSSE_temp.csv", row.names = FALSE)