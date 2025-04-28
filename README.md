# Viral evolution prediction identifies broadly neutralizing antibodies against existing and prospective SARS-CoV-2 variants

This repo contains scripts for reproduce all analyses and figures in the manuscript "Viral evolution prediction identifies broadly neutralizing antibodies against existing and prospective SARS-CoV-2 variants" (Jian et al. bioRxiv 2024)

To reproduce the figures in the manuscript (**Note it is always assumed that the current working directory the directory containing the script while running each script.**):

1. Clone this repo and download the source data files from Zenedo (link to be updated).
2. Decompress the downloaded source data into `source_data/`.
3. Run `scripts/prepare_all_data.py`. This script generates the processed source data for each figure and save them to `processed_source_data` (which should be consistent with the Source Data published along with the manuscript).
4. Run the corresponding plotting script in `plot_scripts/`.

This repo may be updated in the future, but you can always download the orignial version from Zenedo.

The following dependencies are required to run the scripts:

- python: pandas, numpy, scipy, statsmodels, logomaker
- R: ggplot2, tidyverse, ggrastr, ComplexHeatmap

Some Extended Data Figures are generated at the same time as the corresponding main figures.
