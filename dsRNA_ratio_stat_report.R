#!/usr/bin/env Rscript

#' dsRNA Ratio Statistical Report Generator
#' 
#' This script calculates and reports statistics for double-stranded RNA (dsRNA) reads
#' from nucleased and total RNA samples. It processes stat files to compute the ratio
#' of dsRNA reads to total reads for each sample.
#'
#' @author Hess Hakimjavadi
#' @version 1.0.0
#' @date 2023-04-06
#'
#' Publication: Trafficking of mitochondrial double-stranded RNA from mitochondria to the cytosol
#' DOI: 10.26508/lsa.202302396

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

#' Configuration settings for the analysis
#' @return List of configuration parameters
get_config <- function() {
  list(
    base_dir = "data",
    nucleased_dir = "nucleased/filtered_stats",
    total_RNA_dir = "total_RNA/filtered_stats",
    output_dir = "results/dsRNA_stats",
    file_pattern = "*_chrM_stats.txt"
  )
}

#' Create necessary output directories
#' @param config Configuration list containing directory paths
create_directories <- function(config) {
  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }
}

#' Get list of stat files for processing
#' @param config Configuration list
#' @param type Type of sample ("nucleased" or "total_RNA")
#' @return Character vector of file paths
get_stat_files <- function(config, type) {
  base_path <- if(type == "nucleased") {
    file.path(config$base_dir, config$nucleased_dir)
  } else {
    file.path(config$base_dir, config$total_RNA_dir)
  }
  
  files <- list.files(
    path = base_path,
    pattern = config$file_pattern,
    full.names = TRUE
  )
  
  if(length(files) == 0) {
    warning(sprintf("No stat files found in %s", base_path))
  }
  
  return(files)
}

#' Extract numeric value from stat file string
#' @param text Text string containing the numeric value
#' @param pattern Pattern to match for extraction
#' @return Numeric value
extract_numeric <- function(text, pattern) {
  tryCatch({
    as.numeric(str_extract(text, "\\d+"))
  }, error = function(e) {
    warning(sprintf("Failed to extract number from: %s", text))
    return(NA)
  })
}

#' Process a single stat file
#' @param file_path Path to the stat file
#' @return Named vector with processed statistics
process_stat_file <- function(file_path) {
  tryCatch({
    # Read the stat file
    data <- read_csv(file_path, col_names = FALSE, show_col_types = FALSE)
    
    # Extract values using regex
    total_reads <- extract_numeric(data[[2]][1], "number of reads:")
    kept_reads <- extract_numeric(data[[3]][1], "number of kept reads:")
    
    # Calculate dsRNA reads and ratio
    ds_reads <- total_reads - kept_reads
    ds_reads_ratio <- ds_reads / total_reads
    
    # Return results as named vector
    c(
      sample_name = basename(file_path),
      total_reads = total_reads,
      ds_reads = ds_reads,
      ds_reads_ratio = ds_reads_ratio,
      sample_type = if(str_detect(file_path, "nucleased")) "nucleased" else "total_RNA"
    )
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", file_path, e$message))
    return(NULL)
  })
}

#' Process multiple stat files and combine results
#' @param stat_files Vector of file paths
#' @return Data frame with combined statistics
process_stat_files <- function(stat_files) {
  # Process all files
  results <- lapply(stat_files, process_stat_file)
  
  # Remove NULL results (failed processing)
  results <- results[!sapply(results, is.null)]
  
  # Combine results into a data frame
  stats_df <- as.data.frame(do.call(rbind, results))
  
  # Convert numeric columns
  stats_df <- stats_df %>%
    mutate(across(c(total_reads, ds_reads, ds_reads_ratio), as.numeric))
  
  return(stats_df)
}

#' Generate summary statistics
#' @param stats_df Data frame containing processed statistics
#' @return Data frame with summary statistics
generate_summary <- function(stats_df) {
  summary_stats <- stats_df %>%
    group_by(sample_type) %>%
    summarise(
      mean_ratio = mean(ds_reads_ratio, na.rm = TRUE),
      sd_ratio = sd(ds_reads_ratio, na.rm = TRUE),
      median_ratio = median(ds_reads_ratio, na.rm = TRUE),
      n_samples = n()
    )
  
  return(summary_stats)
}

#' Save results to files
#' @param stats_df Data frame with processed statistics
#' @param summary_df Data frame with summary statistics
#' @param config Configuration list
save_results <- function(stats_df, summary_df, config) {
  # Create timestamp for unique file names
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Save detailed results
  write_csv(
    stats_df,
    file.path(config$output_dir, sprintf("dsRNA_stats_detailed_%s.csv", timestamp))
  )
  
  # Save summary results
  write_csv(
    summary_df,
    file.path(config$output_dir, sprintf("dsRNA_stats_summary_%s.csv", timestamp))
  )
}

#' Main execution function
#' @return List containing processed results
main <- function() {
  # Get configuration
  config <- get_config()
  
  # Create output directories
  create_directories(config)
  
  # Get stat files
  nucleased_files <- get_stat_files(config, "nucleased")
  total_RNA_files <- get_stat_files(config, "total_RNA")
  all_files <- c(nucleased_files, total_RNA_files)
  
  # Process all files
  message("Processing stat files...")
  stats_df <- process_stat_files(all_files)
  
  # Generate summary statistics
  message("Generating summary statistics...")
  summary_stats <- generate_summary(stats_df)
  
  # Save results
  message("Saving results...")
  save_results(stats_df, summary_stats, config)
  
  # Return results
  list(
    detailed_stats = stats_df,
    summary_stats = summary_stats
  )
}

# Execute main function if script is run directly
if (!interactive()) {
  results <- main()
}