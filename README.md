# LactoPRECISE
LactoPRECISE: A comprehensive framework exploring Lactobacillus reuteri's genetic diversity, regulatory networks, and functional annotations for personalized medicine and microbial therapeutics.

This repository presents a computational workflow to compute and characterize all iModulons for a selected organism. This occurs in five steps:

Gather all publicly available RNA-seq data for the organism (Step 1)
Process the RNA-seq data (Step 2)
Inspect data to identify high-quality datasets (Step 3)
Compute iModulons (Step 4)
Characterize iModulons using PyModulon (Step 5)

iModulons are independently-modulated groups of genes that are computed through Independent Component Analysis (ICA) of a gene expression dataset. To learn more about iModulons or explore published iModulons, visit iModulonDB or see our publications for Escherichia coli, Staphylococcus aureus, or Bacillus subtilis.

Setup
Pre-requisite software is listed within each step of the workflow. In addition, we have provided pre-built Docker containers with all the necessary software.

To begin, install Docker and Nextflow.

Cite
A pre-print is being prepared for this tutorial workflow.
