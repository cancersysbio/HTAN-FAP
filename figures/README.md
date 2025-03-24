# Scripts to generate manuscript figures

This repository contains scripts used to generate various figures for the study. The figures are categorized into main figures, extended data figures, and supplementary figures. Each category has its own folder containing the relevant scripts and a README file explaining the contents and usage of each script, including any data files each script requires. 

In general, most scripts in this do not generate additional data files that other scripts depend on, but often require processed data from the `analysis` scripts directory. One notable exception is the exonic mutation burden tables (saved in this repo in `data/_exonic_mutation_burden.tsv`) which are generated in `figure1a.Rmd` and `extended_data_figure1a_e.Rmd` and are required in other figure scripts. In general we recommend downloading/generating the full dataset before running the figure scripts.

## Folder Structure

```
figures
├── README.md
├── extended_data
│ └── README.md
├── main
│ └── README.md
└── supplementary
  └── README.md
```