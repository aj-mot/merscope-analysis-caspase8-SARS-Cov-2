# merscope-analysis-caspase8-SARS-Cov-2

This is code for analysis of MERFISH spatial transcriptomic data related to the 
manuscript: 

Bader et al, Non-apoptotic caspase-8 is critical for orchestrating exaggerated 
inflammation during severe SARS-CoV-2 infection.

The repository contains:

- cellpose_segmentation_vpt - code and data relating to cell segmentation using
  a custom Cellpose 2.0 model.
- r-seurat-analysis - R scripts for downstream analysis of segmented MERSCOPE 
  data.
- covid_infect_distance_calculation - Python code associated with calculating distances between cells and infected areas using output from Qupath

Please see the README files contained in each of the above directories.
