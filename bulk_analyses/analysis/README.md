# Scripts for data analysis

This folder contains scripts used to process and analyze sequencing data. Code for directly generating figures is not here, but stored in the top-level `figures` directory. Each subdirectory has its own folder containing the relevant scripts and a README file explaining the contents and usage of each script, including any data files each script requires. Note that all analysis scripts here are for the **bulk** HTAN sequencing data, unless specifically stated that they process the single-crypt WGS data.

## Folder Structure

```
analysis
├── README.md
├── statistics_in_maintext.ipynb
├── copy_number
│ └── README.md
├── mutation_postprocessing
│ └── README.md
└── ppVAF_calculations
│ └──  README.md
└── telomere
```

## Description of scripts

`statistics_in_maintext.ipynb`

Requires:
* filtered HTAN WGS and WES mafs (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf` and `HTAN_WES_filtered_ppVAFs.maf`)
* monoclonal/polyclonal classifications (HTAN calls in this repo at `data/clonal_count_estimation/clonal_SNVs_WES_WGS.csv`)
* annotated table of fraction genome altered in HTAN samples, saved in `data/copy_number/fraction_genome_altered/HTAN_FGA_combined.csv`
* gene-level copy number calls for the HTAN datasets (in this repo at `data/copy_number/gene_CN_calls/HTAN_WGS_gene_CNs.tsv` and `data/copy_number/gene_CN_calls/HTAN_WES_gene_CNs.tsv`)
* WGD information (FACETS HTAN calls in this repo at `data/genome_doubling/*_doubled.tsv`)
* COADREAD driver gene list (in this repo at `data/resource/PanCanDrivers_COADREAD_Cell2018.csv`)

Generates:
No files, but prints out the numerical data presented in the main text of the manuscript

