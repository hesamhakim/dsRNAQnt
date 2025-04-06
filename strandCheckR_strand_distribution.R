#!/usr/bin/env Rscript

#' Mitochondrial dsRNA Analysis Pipeline - Strand Distribution Analysis
#' 
#' This script analyzes strand-specific RNA sequencing data to quantify double-stranded RNA
#' from mitochondrial genome using the strandCheckR package. It processes BAM files to determine
#' the proportion of reads mapping to both strands and generates visualizations of the coverage
#' distribution.
#'
#' @author Hess Hakimjavadi
#' @version 1.0.0
#' @date 2023-04-06
#'
#' Publication: Trafficking of mitochondrial double-stranded RNA from mitochondria to the cytosol
#' DOI: 10.26508/lsa.202302396

# Load required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(strandCheckR)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(Homo.sapiens)
})

#' Check if required packages are installed
#' @param pkg_list List of required packages
check_packages <- function(pkg_list) {
  missing_pkgs <- pkg_list[!pkg_list %in% installed.packages()[,"Package"]]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }
}

#' Validate input files existence
#' @param file_paths Vector of file paths to check
#' @return Logical indicating if all files exist
validate_input_files <- function(file_paths) {
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop("Missing input files: ", paste(missing_files, collapse = ", "))
  }
  return(TRUE)
}

#' Configure analysis parameters
#' @return List of configuration parameters
get_config <- function() {
  list(
    window_size = 100,
    step_size = 20,
    strand_threshold = 0.7,
    target_chromosome = "chrM",
    output_dir = "results/strand_distribution"
  )
}

#' Create output directories if they don't exist
#' @param base_dir Base directory path
create_output_dirs <- function(base_dir) {
  dirs <- c("results", "results/strand_distribution", "results/plots")
  for (dir in file.path(base_dir, dirs)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
}

#' Process BAM file using strandCheckR
#' @param bam_file Path to BAM file
#' @param config Configuration parameters
#' @return GRanges object with strand information
process_bam_file <- function(bam_file, config) {
  tryCatch({
    strand_data <- getStrandFromBamFile(
      bamFile = bam_file,
      windowSize = config$window_size,
      stepSize = config$step_size,
      threshold = config$strand_threshold,
      chrs = config$target_chromosome
    )
    return(strand_data)
  }, error = function(e) {
    stop("Error processing BAM file: ", bam_file, "\n", e$message)
  })
}

#' Annotate strand data with gene information
#' @param strand_data GRanges object with strand information
#' @return Annotated data frame
annotate_strand_data <- function(strand_data) {
  # Get transcript annotations
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  transcripts <- transcripts(txdb)
  
  # Intersect with features
  annotated_data <- intersectWithFeature(
    windows = strand_data,
    feature = transcripts
  )
  
  return(annotated_data)
}

#' Generate coverage plots
#' @param strand_data Processed strand data
#' @param sample_name Name of the sample
#' @param output_dir Output directory
plot_coverage <- function(strand_data, sample_name, output_dir) {
  pdf(file.path(output_dir, paste0(sample_name, "_coverage.pdf")))
  plotWin(strand_data, title = paste("Coverage Distribution -", sample_name))
  dev.off()
  
  pdf(file.path(output_dir, paste0(sample_name, "_histogram.pdf")))
  plotHist(strand_data, title = paste("Strand Distribution -", sample_name))
  dev.off()
}

#' Main execution function
#' @param bam_files Vector of BAM file paths
#' @param output_dir Base output directory
main <- function(bam_files, output_dir) {
  # Validate inputs and setup
  validate_input_files(bam_files)
  config <- get_config()
  create_output_dirs(output_dir)
  
  # Process each BAM file
  results <- list()
  for (bam_file in bam_files) {
    sample_name <- tools::file_path_sans_ext(basename(bam_file))
    message("Processing sample: ", sample_name)
    
    # Process BAM file
    strand_data <- process_bam_file(bam_file, config)
    
    # Annotate data
    annotated_data <- annotate_strand_data(strand_data)
    
    # Generate plots
    plot_coverage(strand_data, sample_name, file.path(output_dir, "plots"))
    
    # Store results
    results[[sample_name]] <- annotated_data
    
    # Save processed data
    saveRDS(
      annotated_data,
      file = file.path(output_dir, paste0(sample_name, "_processed.rds"))
    )
  }
  
  return(results)
}

# Script execution
if (!interactive()) {
  # Define input parameters
  bam_files <- c(
    "data/nucleased/filtered_bam/sample1.bam",
    "data/nucleased/filtered_bam/sample2.bam"
    # Add more sample files as needed
  )
  
  output_dir <- "results"
  
  # Execute main function
  results <- main(bam_files, output_dir)
}