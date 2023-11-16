# LINEAR MIXED MODELS - CONNECTIVITY AND ENTROPY LSD - DREYDELLET 03/11/2023

# Load packages
library(LMMstar)# Required to fit linear mixed models (Brice's model) 
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(tidyr)

# Load csv file
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data <- read.csv("data_long_withShenThal.csv")


# Prepare columns
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$network <- factor(data$network)
data$time <- factor(data$time)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Remove 'connectivity' and 'network' columns (Schaefer atlas) and keep distinct rows
data <- data %>%
  select(-c(connectivity, network)) %>%
  distinct()



## PLOT FC MATRICES FOR LSD AND PLACEBO

# Extract unique network names
networks <- unique(unlist(strsplit(as.character(data$network_shen_thal), "_")))

# Split dataset according to the session
sub_data <- split(data, data$session)

# Loop through the sub-datasets for each session
for (k in unique(data$session)) {
  
  data <- sub_data[[k]]

  # Create empty matrices for parameter estimate and pvalue
  conn_matrix <- matrix(NA, nrow = length(networks), ncol = length(networks),
                             dimnames = list(networks, networks))

  # Loop through each pair of networks
  for (i in 1:length(networks)) {
    for (j in i:length(networks)) {
      # Extract the pair of networks
      network_pair <- paste(networks[i], networks[j], sep = "_")
      
      # Filter the rows of the table corresponding to the network pair
      network_pair_rows <- data[data$network == network_pair, ]
      
      # Compute the parameters and p-values
      conn_value <- mean(network_pair_rows$connectivity_shen_thal, na.rm = TRUE)

      # Populate the corresponding cells in the matrices
      conn_matrix[networks[i], networks[j]] <- conn_value
      conn_matrix[networks[j], networks[i]] <- conn_value

    }
  }

  # Define a lookup table for the full network names
  lookup_table <- data.frame(abbrev = c("Thalamus", "Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                             full = c("Thalamus", "Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))
  
  # Use match() function to map abbreviated network names to full names
  rownames(conn_matrix) <- lookup_table$full[match(rownames(conn_matrix), lookup_table$abbrev)]
  colnames(conn_matrix) <- lookup_table$full[match(colnames(conn_matrix), lookup_table$abbrev)]
  
  # Convert matrix to data frame
  df <- melt(conn_matrix)
  
  # View the min and maximum FC values 
  max_value <- max(df$value)
  print(paste("Maximum FC value is:", max_value))
  min_value <- min(df$value)
  print(paste("Minimum FC value is:", min_value))

  # Define the order of the x and y axes
  order <- c("Visual 1", "Visual 2", "Visual Association", "Subcortical Cerebellum", "Thalamus", 
             "Motor", "Medial Frontal", 
             "Frontoparietal", "Default Mode")
  
  # Subset the data to only keep values below or on the diagonal for plotting purposes
  df_sub <- subset(df, as.numeric(factor(Var1, levels = order)) >= as.numeric(factor(Var2, levels = order)))
  
  # Plot global_rsfcmatrix for plc and lsd
  plot_name <- paste0("global_rsfcmatrix_shen_",k)
  assign(plot_name, 
         ggplot(df_sub, aes(x = factor(Var1, levels = order), 
                            y = factor(Var2, levels = order), 
                            fill = value)) +  # Removed the '+' at the end of this line
           geom_tile(colour = "black", size = 0.1) + 
           scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                midpoint = 0, limits=c(-0.5, 0.5), 
                                name = "Fisher's r-to-z", 
                                breaks = c(-0.5, 0, 0.5),
                                guide = guide_colorbar(barwidth = 0.9, barheight = 5),
                                oob = squish) +
           theme_bw() +
           xlab("") +
           ylab("") +
           theme(plot.title = element_text(hjust = 0.5, size = 16),
                 axis.text.x = element_text(angle = 25, size = 16, hjust = 1),
                 axis.text.y = element_text(size = 16, hjust = 1),
                 panel.border = element_blank(),
                 panel.grid = element_blank(),
                 legend.position = c(0.1, 0.9), 
                 legend.justification = c(0, 0.8),
                 legend.box.just = "left", 
                 legend.title = element_text(size = 16), 
                 legend.text = element_text(size = 14)) + 
           coord_fixed())  # Tiles are squared
}
# Save
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/rsfc/figures/')
svg("global_rsfcmatrix_shen_ses_lsd.svg")
print(global_rsfcmatrix_shen_ses_lsd)
dev.off()
svg("global_rsfcmatrix_shen_ses_plc.svg")
print(global_rsfcmatrix_shen_ses_plc)
dev.off()
 





## LMM GLOBAL

# Load the necessary packages
library(lme4)   # Required to fit linear mixed models
library(LMMstar)# Required to fit linear mixed models with Brice's model 
library(ggplot2)
library(reshape2)
library(dplyr)

# Load csv file
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data <- read.csv("data_long_withShenThal.csv")

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$network <- factor(data$network)
data$network_shen_thal <- factor(data$network_shen_thal)
data$time <- factor(data$time)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$drugsbefore <- factor(data$drugsbefore) # nb of drug sessions before the current session, can be either 0, 1, 2 or 3
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Remove 'connectivity' and 'network' columns and keep distinct rows
data <- data %>%
  select(-c(connectivity, network)) %>%
  distinct()


# Fit LMM
ls.modelCS <- mlmm(connectivity_shen_thal ~ session + drugsbefore_binary + age + sex + bmi + dataset, by = "network_shen_thal",
                repetition = ~ session|participant, structure = "UN", data = data,
                effects = "sessionses_lsd=0", robust = FALSE)
summary(ls.modelCS, method = "single-step")
plot(ls.modelCS, method = "none")
table <- model.tables(ls.modelCS, columns = add("parameter"), method = "single-step")
table_uncorrected <- model.tables(ls.modelCS, columns = add("parameter"), method = "none")




# summary(ls.modelCS$model$Cont_Cont)

# ls.modelCS.DB <- mlmm(connectivity ~ session + drugsbefore + age + sex + bmi + dataset, by = "network",
#                    repetition = ~ session|participant, structure = "UN", data = data,
#                    effects = "drugsbefore1=0", robust = FALSE)
# plot(ls.modelCS.DB, method = "none")


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
lookup_table <- data.frame(abbrev = c("Thalamus", "Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                           full = c("Thalamus", "Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))

# Use left_join() to join the lookup_table to the original data frame
# by matching the "abbrev" column in lookup_table with "Var1" and "Var2" columns in df
df <- df %>%
  left_join(lookup_table, by = c("Var1" = "abbrev")) %>%
  left_join(lookup_table, by = c("Var2" = "abbrev"))

# Rename the new columns to Var1_full and Var2_full
colnames(df)[6:7] <- c("Network1_full", "Network2_full")

# Define the order of the x and y axes
order <- c("Visual 1", "Visual 2", "Visual Association", "Subcortical Cerebellum", "Thalamus", 
           "Motor", "Medial Frontal", 
           "Frontoparietal", "Default Mode")

# Subset the data to only keep values below or on the diagonal
df_sub <- subset(df, as.numeric(factor(Network1_full, levels = order)) >= as.numeric(factor(Network2_full, levels = order)))

# Plot global_rsfcmatrix
global_rsfcmatrix <- 
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
# Save global_rsfcmatrix
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/rsfc/')
svg("global_rsfcmatrix_shen_thal.svg")
print(global_rsfcmatrix)
dev.off()
ggsave("global_rsfcmatrix_shen_thal.png", plot = global_rsfcmatrix, device = "png")
ggsave("global_rsfcmatrix_shen_thal.eps", plot = global_rsfcmatrix, device = "eps")





## EFFECT SIZES

# Load csv file
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/')
data <- read.csv("data_long_withShenThal.csv")

# Convert gender, network, and session to factors
data$sex <- factor(data$sex)
data$participant <- factor(data$participant)
data$network <- factor(data$network)
data$network_shen_thal <- factor(data$network_shen_thal)
data$time <- factor(data$time)
data$drugsbefore_binary <- factor(data$drugsbefore_binary) # was there a drug session before the current session? can be either 0 or 1
data$drugsbefore <- factor(data$drugsbefore) # nb of drug sessions before the current session, can be either 0, 1, 2 or 3
data$session <- factor(data$session, levels = c("ses-plc","ses-lsd"), labels = c("ses_plc","ses_lsd"))

# Remove 'connectivity' and 'network' columns and keep distinct rows
data <- data %>%
  select(-c(connectivity, network)) %>%
  distinct()


# Subset data to only include relevant columns
data_subset <- data %>% select(participant, session, network_shen_thal, connectivity_shen_thal)

# Split network column into two separate columns
data_subset <- data_subset %>%
  separate(network_shen_thal, into = c("network_1", "network_2"), sep = "_")

# Calculate sample size, mean, and standard deviation for each pair of networks and session
network_summary <- data_subset %>%
  dplyr::group_by(network_1, network_2, session) %>%
  dplyr::summarize(sample_size = sum(!is.na(connectivity_shen_thal)),  # Calculate sample size by counting non-missing values
                   mean_connectivity = mean(connectivity_shen_thal, na.rm = TRUE),
                   sd_connectivity = sd(connectivity_shen_thal, na.rm = TRUE)) %>%
  dplyr::ungroup()


# Calculate effect size for each pair of networks
effect_size <- network_summary %>%
  dplyr::group_by(network_1, network_2) %>%
  dplyr::summarize(cohen_d = (dplyr::first(mean_connectivity[session == "ses_lsd"]) - 
                                dplyr::first(mean_connectivity[session == "ses_plc"])) / 
                     sqrt(((dplyr::first(sample_size[session == "ses_lsd"]) - 1) * 
                             (dplyr::first(sd_connectivity[session == "ses_lsd"]))^2 +
                             (dplyr::first(sample_size[session == "ses_plc"]) - 1) * 
                             (dplyr::first(sd_connectivity[session == "ses_plc"]))^2) /
                            (dplyr::first(sample_size[session == "ses_lsd"]) + 
                               dplyr::first(sample_size[session == "ses_plc"]) - 2))) # pooled variance

effect_size_temp = effect_size

# Define a lookup table for the full network names
lookup_table <- data.frame(abbrev = c("Thalamus", "Frontoparietal", "SubcorticalCerebellum", "DefaultMode", "MedialFrontal", "Motor", "VisualAssociation", "Visual1", "Visual2"),
                           full = c("Thalamus", "Frontoparietal", "Subcortical Cerebellum", "Default Mode", "Medial Frontal", "Motor", "Visual Association", "Visual 1", "Visual 2"))

# Rename networks in effect_size using lookup_table
effect_size$Network1_full <- lookup_table$full[match(effect_size$network_1, lookup_table$abbrev)]
effect_size$Network2_full <- lookup_table$full[match(effect_size$network_2, lookup_table$abbrev)]

# Remove the original abbreviated network name columns
effect_size$network_1 <- NULL
effect_size$network_2 <- NULL

# Define the order of the x and y axes
order <- c("Visual 1", "Visual 2", "Visual Association", "Subcortical Cerebellum", "Thalamus",  
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

# Plot global_cohenmatrix
global_cohenmatrix <- 
  ggplot(effect_size, aes(x = factor(effect_size$Network1_full, levels = order), 
                          y = factor(effect_size$Network2_full, levels = order), 
                          fill = effect_size$cohen_d)) +
  geom_tile(colour = "black", size = 0.1) + # Change the borders around individual tiles
  geom_text(aes(label = stars), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Cohen's d",
                       breaks = c(-0.5, 0, 0.5, 1, 1.5),
                       # limits = c(-0.8, 1.6),  # This line ensures the scale goes from -0.8 to 1.6
                       guide = guide_colorbar(barwidth = 0.9, barheight = 5)) +
  theme_bw() +
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
        axis.text.y = element_text(hjust = 1, size = 16),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = c(0.1, 0.9), # Set the position of the legend
        legend.justification = c(0, 0.8), # Set the justification of the legend
        legend.box.just = "left") + # Set the justification of the legend box)
  coord_fixed()  # Tiles are squared
# Save global_cohenmatrix
setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/rsfc/')
svg("global_cohenmatrix_shen_thal.svg")
print(global_cohenmatrix)
dev.off()
ggsave("global_cohenmatrix_shen_thal.png", plot = global_cohenmatrix, device = "png")
ggsave("global_cohenmatrix_shen_thal.eps", plot = global_cohenmatrix, device = "eps")



# Create table_globalFC

# Add Cohen's d to table
# Create a 'by' column in effect_size_temp using network_1 and network_2
effect_size_temp <- effect_size_temp %>%
  mutate(by = paste(network_1, network_2, sep = "_"))
# Left join table with effect_size_temp based on 'by' column
table <- left_join(table, effect_size_temp[, c("by", "cohen_d")], by = "by")


# Create table_globalFC
table_globalFC <- table %>%
  select(by, cohen_d, estimate, lower, upper, p.value) %>%
  rename(network = by,
         effect_size = cohen_d,
         lower_confint_singlestep = lower,
         upper_confint_singlestep = upper,
         pvalue_singlestep = p.value)
# Add uncorrected pvalues
table_globalFC <- table_globalFC %>%
  left_join(select(table_uncorrected, by, p.value, lower, upper), by = c("network" = "by")) %>%
  rename(
    pvalue = p.value,
    lower_confint = lower,
    upper_confint = upper
  ) %>%
  select(-pvalue_singlestep, pvalue, everything())

# Re-organizing columns for table_globalFC
table_globalFC <- table_globalFC %>%
  select(network, effect_size, estimate, lower_confint, upper_confint, pvalue, lower_confint_singlestep, upper_confint_singlestep, pvalue_singlestep, everything())

setwd('/student/davidreydellet/lsd-basel/data/derivative/analysis/rsfc/')
write.csv(table_globalFC, "table_globalFC_shenthal.csv")











# Save svg (for some reason doesn't work the way I did it in loops)
svg("global_rsfcmatrix_ses_plc_shen_thal.svg")
print(global_rsfcmatrix_shen_thal_ses_plc)
dev.off()
svg("global_rsfcmatrix_ses_lsd_shen_thal.svg")
print(global_rsfcmatrix_shen_thal_ses_lsd)
dev.off()
svg("dataset_rsfcmatrix_d1_shen_thal.svg")
print(dataset_rsfcmatrix_d1_shen_thal)
dev.off()
svg("dataset_rsfcmatrix_d2_shen_thal.svg")
print(dataset_rsfcmatrix_d2_shen_thal)
dev.off()
svg("dataset_rsfcmatrix_d3_shen_thal.svg")
print(dataset_rsfcmatrix_d3_shen_thal)
dev.off()
svg("drugnaive_rsfcmatrix_0_shen_thal.svg")
print(drugnaive_rsfcmatrix_0_shen_thal)
dev.off()
svg("drugnaive_rsfcmatrix_1_shen_thal.svg")
print(drugnaive_rsfcmatrix_1_shen_thal)
dev.off()

svg("dataset_cohendmatrix_study1_shen_thal.svg")
print(dataset_cohendmatrix_study1_shen_thal)
dev.off()
svg("dataset_cohendmatrix_study2_shen_thal.svg")
print(dataset_cohendmatrix_study2_shen_thal)
dev.off()
svg("dataset_cohendmatrix_study3_shen_thal.svg")
print(dataset_cohendmatrix_study3_shen_thal)
dev.off()
svg("drugnaive_cohendmatrix_drugexposed_shen_thal.svg")
print(drugnaive_cohendmatrix_drugexposed_shen_thal)
dev.off()
svg("drugnaive_cohendmatrix_drugnaive_shen_thal.svg")
print(drugnaive_cohendmatrix_drugnaive_shen_thal)
dev.off()
