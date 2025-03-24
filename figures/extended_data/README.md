## Extended Data Figures

This folder contains scripts for generating the Extended Data figures. The scripts in this directory can be run in any order, assuming the dependencies generated in previous analysis steps.


### Description of scripts

`extended_data_figure3.ipynb`
Plots and saves all panels to make Extended Data Fig. 3
Requires:
* filtered and merged copy number calls (HTAN calls in this repo at `data/copy_number/*_CN_filtered_merged.tsv`)
* fraction genome altered tables (HTAN calls in this repo at `data/fraction_genome_altered/*_fga.tsv`)
* WGD information (FACETS HTAN calls in this repo at `data/genome_doubling/*_doubled.tsv`)
* GRCh38 centromere/telomere coordinates (in this repo at `data/resource/hg38.UCSC.centromere.telomere.encode.bed`)

Generates:
* plots for ED Fig. 3a-f

