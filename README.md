Based on the research reports and the provided information, I'll help create a comprehensive README file and provide suggestions for revising the R scripts. Let's start with the README file:

```markdown
# Mitochondrial dsRNA Analysis Pipeline

## Overview
This repository contains R scripts for analyzing double-stranded RNA (dsRNA) from mitochondrial genome sequencing data, as described in the publication "Trafficking of mitochondrial double-stranded RNA from mitochondria to the cytosol" (PMC11220484).

The pipeline quantifies the proportion of reads originating from double-stranded RNA using the strandCheckR package and generates visualizations of mitochondrial read coverage distributions.

## Publication
- **Title**: Trafficking of mitochondrial double-stranded RNA from mitochondria to the cytosol
- **Link**: https://pmc.ncbi.nlm.nih.gov/articles/PMC11220484/
- **Citation**: [Add citation details]

## Scripts
The repository contains two main R scripts:

1. `strandCheckR_strand_distribution02.R`
   - Filters dsRNA reads
   - Visualizes distribution of mitochondrial read coverages mapping to +/- strands
   - Uses sliding window approach (100bp windows, 20bp steps)

2. `dsRNA_ratio_stat_report.R`
   - Calculates ratio of estimated dsRNA to total recovered RNA reads
   - Generates statistical reports for downstream analysis

## Dependencies
### R Packages
- tidyverse
- strandCheckR (To & Pederson, 2021)
- TxDb.Hsapiens.UCSC.hg19.knownGene
- TxDb.Hsapiens.UCSC.hg20.knownGene
- TxDb.Hsapiens.UCSC.hg38.knownGene
- Homo.sapiens

### Installation
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("strandCheckR", "tidyverse"))
```

## Input Data Requirements
### BAM Files
- Aligned RNA sequencing reads
- Must be sorted and indexed
- Should contain mitochondrial chromosome (chrM) alignments

### Directory Structure
```
project/
├── data/
│   ├── nucleased/
│   │   ├── filtered_bam/
│   │   └── filtered_stats/
│   └── total_RNA/
│       ├── filtered_bam/
│       └── filtered_stats/
├── scripts/
└── results/
```

## Usage
1. **Prepare Input Data**
   - Place BAM files in appropriate directories
   - Ensure BAM files are sorted and indexed

2. **Run strandCheckR Analysis**
```r
Rscript strandCheckR_strand_distribution02.R
```

3. **Calculate dsRNA Ratios**
```r
Rscript dsRNA_ratio_stat_report.R
```

## Parameters
- **Window Size**: 100bp
- **Step Size**: 20bp
- **Strand Threshold**: >70% reads from single strand for classification
- **Chromosome**: chrM (mitochondrial)

## Output
1. **Strand Distribution Data**
   - TSV/CSV files containing strand distribution statistics
   - Coverage data for visualization

2. **Statistical Reports**
   - dsRNA ratio calculations
   - Sample-wise statistics

## Visualization
The output data can be imported into Tableau for visualization as shown in Figure 8 of the publication:
- Distribution of aligned reads across mitochondrial genome
- Proportion of dsRNA reads per sample
- Coverage plots for +/- strand mapping

## Contributing
Please feel free to submit issues and pull requests for improvements.

## License
[Add appropriate license information]

## Contact
[Add contact information]
```

For the revised R scripts, I'll provide the structure for the main script (`strandCheckR_strand_distribution02.R`). Would you like me to provide the complete revised versions of both scripts with detailed comments?

The revised scripts would include:
1. Comprehensive header documentation
2. Function-level documentation
3. Inline comments explaining complex operations
4. Error handling
5. Consistent code formatting
6. Parameter documentation
7. Output documentation

Would you like me to proceed with providing the revised versions of both scripts?
