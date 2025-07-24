# BRB-seq NEXTFLOW Pipeline README

## Overview

The NEXTFLOW pipeline is designed for the comprehensive pre-processing of BRB-seq data. This pipeline facilitates the transition from raw sequencing data to processed outputs that are ready for downstream analysis. The primary functions of the pipeline include demultiplexing, alignment, quality control, data summarization and reporting.

## Major Steps in the Pipeline

1. **Read Alignment and Demultiplexing**:
   - Tool: **STAR solo** (version 2.7.9a)
   - Description: Aligns raw sequencing reads and demultiplexes them to individual samples.

2. **BAM Indexing**:
   - Tool: **samtools** (version 1.12)
   - Description: Indexes the BAM files generated during alignment to facilitate quick data retrieval.

3. **Raw Sequencing Read QC**:
   - Tool: **FastQC** (version 0.11.9)
   - Description: Performs quality control checks on raw sequencing reads to ensure data integrity.

4. **Additional QC**:
   - Tool: **RSeQC** (version 4.0.0)
   - Description: Conducts read distribution abalysis on the aligned reads.

5. **BigWigs Generation**:
   - Tool: **deepTools** (version 3.5.1)
   - Description: Converts BAM files to BigWig format for easy visualization of read coverage.

6. **Counts Matrix Generation and Deduplication**:
   - Process: Generates a raw read count matrix and uses unique molecular identifiers (UMIs) to create a deduplicated counts matrix.

7. **Alignment Summary**:
   - Custom Scripts: Summarizes alignment statistics and results.

8. **HTML Report Generation**:
   - Custom Scripts: Builds an official Alithea Genomics' HTML report summarizing all results.

## Tools and Versions

The pipeline utilizes the following tools and software versions:

- **Python**: 3.9.6
- **STAR**: 2.7.9a
- **R**: 4.1.1
- **tidyverse**: 1.3.1
- **writexl**: 1.4.0
- **readxl**: 1.4.0
- **samtools**: 1.12
- **FastQC**: 0.11.9
- **OpenJDK**: 8.0.152
- **ggplot2**: 3.4.0
- **stringr**: 1.4.0
- **data.table**: 1.14.2
- **reshape2**: 1.4.4
- **RSeQC**: 4.0.0
- **deepTools**: 3.5.1
- **markdown**: 1.1
- **plotly**: 4.10.0
- **kableExtra**: 1.3.4
- **tinytex**: 0.35
- **openpyxl**: 3.0.9
- **ggrepel**: 0.9.3
- **ggpubr**: 0.4.0

## Containerized Environment

The NEXTFLOW pipeline runs in a containerized Miniconda3 environment, ensuring full reproducibility of results. The containerization guarantees that all dependencies and software versions are consistent across different runs, eliminating discrepancies caused by differing system configurations.

## Running the Pipeline

To run the NEXTFLOW pipeline, follow these general steps:

1. **Setup the Containerized Environment**:
   Ensure Miniconda3 is installed and the appropriate environment is configured with the specified tools and versions.

2. **Input Data**:
   Provide the raw sequencing data in fastq format.

3. **Execute the Pipeline**:
   Run the pipeline script, which will carry out all the steps from read alignment to report generation.

4. **Output Files**:

Upon completion, the pipeline generates the following output folders and files:

- **bam/**: Aligned reads on BAM format, STARsolo output.

- **bigwigs/**: Contains BigWig files for each sample.
  
- **fastqc/**: Contains FastQC reports and meta data for raw sequencing reads quality control.

- **plots/**: Contains plate view plots with information on sample demultiplexing rate in percentage.

- **reports/**: Contains HTML report(s) about the quality control of sequencing data.

- **count_matrix/**: Contains count matrices in various formats:

    - `[FASTQ_NAME]_matrix.umi.counts.txt`: Deduplicated UMI counts across all wells in a plate, including both pooled and non-pooled samples, with well IDs used for header.
    
    - `[FASTQ_NAME].umi.counts.sampleIds.txt`: Deduplicated UMI counts across pooled samples mentioned in the sample submission form (SSF), with sample IDs used for header.
    
    - `[FASTQ_NAME].umi.counts.wells.txt`: Deduplicated UMI counts across pooled samples mentioned in the sample submission form (SSF), with well IDs used for header.
    
    - `[FASTQ_NAME].read.counts.wells.detailed.txt`: Non-deduplicated read counts across all wells in a plate, including both pooled and non-pooled samples, with well IDs used for header. Includes details about alignment per sample.
    
    - `[FASTQ_NAME].read.counts.txt`: Non-deduplicated read counts across all wells in a plate, including both pooled and non-pooled samples, with well IDs used for header.
    
    - `[FASTQ_NAME].read.counts.wells.txt`: Non-deduplicated read counts across pooled samples mentioned in the sample submission form (SSF), with well IDs used for header.
    
    - `[FASTQ_NAME].read.counts.sampleIds.detailed.txt`: Non-deduplicated read counts across pooled samples mentioned in the sample submission form (SSF), with sample IDs used for header. Includes details about alignment per sample.
    
    - `[FASTQ_NAME].read.counts.sampleIds.txt`: Non-deduplicated read counts across pooled samples mentioned in the sample submission form (SSF), with sample IDs used for header.

- **qc/**: Contains quality control in 

    - `[FASTQ_NAME]_matrix.alignment.by.sample.txt`: Per-sample meta-information and quality control statistics.

    - `distribution.merged.txt` and `distribution_reads.pdf`: Per-plate RSeQC read distribution analysis output.

    - `summary.sequencing.txt`: Per-plate summary file containing alignment statistics generated by custom scripts..

- **Summary.Alignment.txt**: Per-plate summary file containing alignment statistics generated by STARsolo.

- **md5sums.txt**: Contains checksums for all generated files.




## Contact and Support

For any questions or support regarding the BRB-seq NEXTFLOW pipeline, please contact the Alithea Genomics Bioinformatics team.

---

This README provides a concise overview of the NEXTFLOW pipeline for BRB-seq data preprocessing. The pipeline's structured approach and use of robust tools ensure high-quality and reproducible results, facilitating efficient downstream analysis.