## Extended Data Figures

This folder contains scripts for generating the Extended Data figures. The scripts in this directory can be run in any order, assuming the dependencies generated in previous analysis steps.


### Description of scripts

`extended_data_figure1a_e.Rmd`

Plots and saves oncoplots to make Extended Data Fig. 1a and lollipop plots for Extended Data Fig. 1e (for external PUTH and SCORT cohorts).

Requires:
* filtered PUTH and SCORT mafs (not provided, would have to generate these yourself following the preprocessing pipeline described in the Methods)
* PUTH and SCORT metadata (in this repo at `data/metadata/PUTH_metadata.tsv` and `data/metadata/SCORT_metadata.tsv`)

Generates:
* PUTH and SCORT oncoplots (ED Fig. 1a)
* PUTH and SCORT APC lollipop plots (ED Fig. 1e)
* PUTH and SCORT exonic mutation burden tables (`data/PUTH_exonic_mutation_burden.tsv` and `data/SCORT_exonic_mutation_burden.tsv`)

`extended_data_figure1b_d.Rmd`

Plots and saves exonic mutation burden boxplots in Extended Data Fig. 1b-d and computes p-values with the Wilcoxon rank sums test.

Requires:
* HTAN exonic mutation burden table at `data/HTAN_exonic_mutation_burden.tsv`
* PUTH and SCORT exonic mutation burden tables (`data/PUTH_exonic_mutation_burden.tsv` and `data/SCORT_exonic_mutation_burden.tsv`)

Generates:
* HTAN, PUTH, and SCORT somatic mutation boxplots with p-values included (ED Fig. 1b-d)

`extended_data_figure2.Rmd` 

Plots and saves APC and KRAS mutation information bar plots for ED Fig 2.

Requires:
* HTAN APC and KRAS mutation information table (in this repo at `data/HTAN_APC_KRAS_hits_matrix_wgs_wes.tsv`)

Generates:
* HTAN bar plots of APC and KRAS mutation prevalence (ED Fig. 2)

`extended_data_figure3.ipynb`

Plots and saves all panels to make Extended Data Fig. 3

Requires:
* filtered and merged copy number calls (HTAN calls in this repo at `data/copy_number/*_CN_filtered_merged.tsv`)
* fraction genome altered tables (HTAN calls in this repo at `data/fraction_genome_altered/*_fga.tsv`)
* WGD information (FACETS HTAN calls in this repo at `data/genome_doubling/*_doubled.tsv`)
* GRCh38 centromere/telomere coordinates (in this repo at `data/resource/hg38.UCSC.centromere.telomere.encode.bed`)
* Table S2 with HTAN metadata (in this repo at `data/Table_S2.csv`)
* External cohort exonic mutation burden tables (in this repo at `data/PUTH_exonic_mutation_burden.tsv` and `data/SCORT_exonic_mutation_burden.tsv`)

Generates:
* plots for ED Fig. 3a-f
* Annotated table of fraction genome altered in HTAN samples, saved as HTAN_FGA_combined.csv

`extended_data_figure4.ipynb`

Plots and saves all panels to make Extended Data Fig. 4

Requires:
* filtered HTAN WGS and WES mafs (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf` and `HTAN_WES_filtered_ppVAFs.maf`)
* monoclonal/polyclonal classifications (HTAN calls in this repo at `data/clonal_count_estimation/clonal_SNVs_WES_WGS.csv`)
* Sequenza and FACETS estimated sample purity tables (in this repo at `data/HTAN_WGS_facets_sequenza_purities.csv` and `data/HTAN_WES_facets_sequenza_purities.csv`)
* scATAC-seq cell fraction table (in this repo at `data/scATACseq_annotations/scATAC_celltype_fracs.csv`)

Generates:
* all plots for ED Fig. 4

`extended_data_figure5.Rmd` and `extended_data_figure5.html`

Computes COSMIC SBS signature decompositions and plots trinucleotide context mutation frequency plots and signature decompositions in Extended Data Fig. 5.

Requires:
* filtered HTAN WGS maf (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf`)

Generates:
* all panels for ED Fig. 5
* Trinucleotide context data tables for the HTAN WGS data (in this repo at `data/htan_fap_tnn_matrix_wgs_per-sample.mat` and `data/htan_fap_tnn_matrix_wgs_per-stage.mat`)
