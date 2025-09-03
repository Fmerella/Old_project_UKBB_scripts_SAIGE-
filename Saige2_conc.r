# Load required libraries
library(data.table)
library(dplyr)

# Define the directory containing the files and the pattern to match
file_dir <- "WGS_Pain_PostJR_FMerella:/Saige2_chrom7/"
file_pattern <- "Saige2_results.*\\.txt$"

# Get a list of all relevant files with full paths
all_files <- list.files(file_dir, pattern = file_pattern, full.names = TRUE)

# Check the number of files in the list
num_files <- length(all_files)

# Print the number of files
print(paste("Number of files:", num_files))

# Initialize an empty data frame to hold concatenated results
ConcSaige2txt <- data.frame()  # Adjust columns as necessary

# Loop over each file, read the data, filter, and combine
for (file_path in all_files){
  # Read the file into a temporary data frame
  temp_dt <- fread(file_path)

  # Filter for low p-values and convert to a data frame if needed
  filtered_dt <- temp_dt %>%
    dplyr::filter(p.value < 5e-8) %>%
    as.data.frame()

  # Concatenate the filtered data into the main results data frame
  ConcSaige2txt <- bind_rows(ConcSaige2txt, filtered_dt)

  # Print the number of rows read for debugging
  print(paste("Read", nrow(filtered_dt), "rows from", file_path))
}

# Print the number of rows in the concatenated data frame
print(paste("Total rows in concatenated data:", nrow(ConcSaige2txt)))

# Optionally, save the concatenated data frame to a file
# write.csv(ConcSaige2txt, "Concatenated_Saige2_Results.csv", row.names = FALSE)
