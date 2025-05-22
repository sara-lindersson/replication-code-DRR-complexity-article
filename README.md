# Replication Code for "Unpacking and navigating complexity in disaster risk reduction"
This repository contains replication code for the article:
__Lindersson et al. (2025).__ _Unpacking and Navigating Complexity in Disaster Risk Reduction._ Manuscript in review.

## Repository Overview
This repository includes R notebooks for generating the statistical analysis and figures presented in the article. To enhance readability, corresponding markdown files (`.md`) are provided as well.

## Data Requirements
To run script `01`, the following files are needed in the raw-data directory: 
+ `initial-sample.csv`: An data table of the initial sample of disaster risk reduction cases (n=76). This dataset will be made available as supplementary information with the article as `dataset S1`.
+ __Geospatial Data__: Shapefiles from [GADM](https://gadm.org/). This study used GADM v3.6, downloaded on _2024-10-03_.

To run scripts numbered as `02` and `03`, ensure you have the `source_data`-files located in the main directory. These files are outputs from script `01`. `source_data.csv` will also be made available as supplementary information with the article, as `dataset S2`.

## Contact
For questions, please contact: __Sara Lindersson__
