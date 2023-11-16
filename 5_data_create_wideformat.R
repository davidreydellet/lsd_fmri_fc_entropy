# This code transforms the csv files from long formats to a single wide format data_wide.csv
library(tidyr)
library(dplyr)

# Set the working directory to where the files are located
setwd("/student/davidreydellet/lsd-basel/data/derivative/analysis/")

# Read the CSV file
data_long <- read.csv("data_long_withShenThal.csv")

# Remove unwanted columns and transform 'network_shen_thal'
data_long <- data_long %>% 
  select(-network, -connectivity) %>%
  mutate(new_col_name = paste0("FC_", gsub("_", "-", network_shen_thal))) %>% distinct()

# Pivot to wide format using 'participant', 'session', and all other columns as the identifiers
data_wide <- data_long %>%
  pivot_wider(
    id_cols = setdiff(names(data_long), c("network_shen_thal", "connectivity_shen_thal", "new_col_name")),
    names_from = new_col_name,
    values_from = connectivity_shen_thal
  )


# Load a long format with DCC
data_dcc <- read.csv("data_long_withDCC.csv")

# Remove unwanted columns (previous FC analysis with Schaefer)
data_dcc <- data_dcc %>% 
  select(-network, -connectivity) %>% distinct()

# Create new columns based on 'network_shen' strings and fill them with values from 'DCC'
data_dcc_wide <- data_dcc %>%
  mutate(new_col_name = paste0("DCC_", gsub("_", "-", network_shen))) %>%
  select(-network_shen) %>%
  pivot_wider(
    names_from = new_col_name,
    values_from = DCC
  )


# Filter out only the columns that start with "DCC" along with the identifiers
dcc_columns <- c("participant", "session", grep("^DCC_", names(data_dcc_wide), value = TRUE))

# Select only the DCC columns and identifiers from data_dcc_wide
data_dcc_selected <- select(data_dcc_wide, all_of(dcc_columns))

# Join the two data frames based on participant and session
data_wide <- left_join(data_wide, data_dcc_selected, by = c("participant", "session"))


# Load a long format with MSSE
data_msse <- read.csv("data_long_withMSSE.csv")
data_msse <- data_msse %>% 
  select(-network, -connectivity) %>% distinct()

# Create new column names based on 'scale_MSSE' and 'network_shen_MSSE'
data_msse_wide <- data_msse %>%
  mutate(new_col_name = paste0("MSSE_Scale", scale_MSSE, gsub("_", "", network_shen_MSSE))) %>%
  select(-scale_MSSE, -network_shen_MSSE) %>%
  pivot_wider(
    names_from = new_col_name,
    values_from = avg_MSSE
  )

# Filter out only the columns that start with "MSSE" along with the identifiers
msse_columns <- c("participant", "session", grep("^MSSE_", names(data_msse_wide), value = TRUE))

# Select only the MSSE columns and identifiers from data_msse_wide
data_msse_selected <- select(data_msse_wide, all_of(msse_columns))

# Join the two data frames based on participant and session
data_wide <- left_join(data_wide, data_msse_selected, by = c("participant", "session"))

# Rename for clarity
data_wide <- data_wide %>%
  rename(path_length_distribution_entropy = geodesic_entropy)

# Save
library(readr)
write_csv(data_wide, "data_wide.csv")