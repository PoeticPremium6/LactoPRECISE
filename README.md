# LactoPRECISE
LactoPRECISE: A comprehensive framework exploring Lactobacillus reuteri's genetic diversity, regulatory networks, and functional annotations for personalized medicine and microbial therapeutics.

![msystems 01257-23 f001](https://github.com/user-attachments/assets/c9f89cc5-7539-4d2a-b807-b7b0e5e7a420)


This repository presents a computational workflow to compute and characterize all iModulons for a selected organism. This occurs in five steps:

Gather all publicly available RNA-seq data for the organism (Step 1)
Process the RNA-seq data (Step 2)
Inspect data to identify high-quality datasets (Step 3)
Compute iModulons (Step 4)
Characterize iModulons using PyModulon (Step 5)

iModulons are independently-modulated groups of genes that are computed through Independent Component Analysis (ICA) of a gene expression dataset. To learn more about iModulons or explore published iModulons, visit iModulonDB or see our publications for Escherichia coli, Staphylococcus aureus, or Bacillus subtilis.

The original research to reconstruct LactoPRECISE has been published in mSystems (2024):
https://journals.asm.org/doi/10.1128/msystems.01257-23
