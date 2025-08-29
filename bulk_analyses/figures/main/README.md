## Main Figures

This folder contains scripts for generating the main text figures presented in the study.

### Description of scripts

`figure1b.Rmd` and `figure1b.html`

Plots and saves oncoplot to make main text Fig. 1b. NB: also generates and saves the HTAN APC mutation lollipop plot in Extended Data Fig. 1e as well as the exonic burden table for the HTAN cohort.

Requires:
* filtered HTAN WGS and WES mafs (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf` and `HTAN_WES_filtered_ppVAFs.maf`)
* gene-level HTAN WGS and WES copy number calls (in this repo at `data/copy_number/gene_CN_calls/HTAN_WGS_gene_CNs.tsv` and `data/copy_number/gene_CN_calls/HTAN_WES_gene_CNs.tsv`)
* FACETS whole genome doubling calls for the HTAN datasets (in this repo at `data/copy_number/genome_doubling/HTAN_WGS_doubled.tsv/` and `data/copy_number/genome_doubling/HTAN_WES_doubled.tsv/`)
* Table S2 with HTAN metadata (in this repo at `data/Table_S2.csv`)

Generates:
* HTAN oncoplot (Fig. 1b)
* HTAN APC lollipop plot (ED Fig. 1e)
* HTAN exonic mutation burden table at `data/HTAN_exonic_mutation_burden.tsv`

`figure1c.Rmd` and `figure1c.html`

Plots and saves upset plot shown in main text Fig. 1c. 

Requires:
* filtered HTAN WGS and WES mafs (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf` and `HTAN_WES_filtered_ppVAFs.maf`)
* Table S2 with HTAN metadata (in this repo at `data/Table_S2.csv`)

Generates:
* HTAN upset plot (Fig. 1c)
* HTAN APC and KRAS mutation information table (in this repo at `data/HTAN_APC_KRAS_hits_matrix_wgs_wes.tsv`)

`figure1dg.ipynb`

Plots and saves panels to make main text Fig. 1d-g

Requires:
* filtered HTAN WGS and WES mafs (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf` and `HTAN_WES_filtered_ppVAFs.maf`)
* APC/KRAS mutation marginalized posterior ppVAF distributions (found in Zenodo dataset as `*_marginals_WGS.csv` and `*_marginals_WES.csv`)

Generates:
* plots for Fig. 1d-g

`figure2.ipynb`

Plots and saves panels to make main text Fig. 2 that are not schematics (i.e., Fig. 2c-f)

Requires:
* filtered HTAN WGS maf (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf`)
* monoclonal/polyclonal classifications (HTAN calls in this repo at `data/clonal_count_estimation/clonal_SNVs_WES_WGS.csv`)

Generates:
* plots for Fig. 2c-f

`figure3.ipynb`

Plots and saves all panels to make main text Fig. 3. Also generates and prints p-values comparing APC/KRAS ppVAFs between monoclonal and polyclonal HTAN WGS samples featured in the main text.

Requires:
* filtered HTAN WGS maf (found in Zenodo dataset as `HTAN_WGS_filtered_ppVAFs.maf`)
* monoclonal/polyclonal classifications (HTAN calls in this repo at `data/clonal_count_estimation/clonal_SNVs_WES_WGS.csv`)

Generates:
* plots for Fig. 3a-d
* Wilcoxon rank sum p-values comparing APC/KRAS ppVAFs between monoclonal and polyclonal samples

