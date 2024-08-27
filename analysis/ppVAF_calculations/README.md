# Calculating ppVAF values using single-cell derived purity distributions

The scripts in this directory take maf files from WGS or WES data and annotate them with purity and ploidy adjusted variant allele frequencies (ppVAFs) for each mutation. To reproduce the analyses and figures in the manuscript, you will use these scripts to process WGS and WES maf files from our FAP cohort that have been filtered to remove mutations shared between samples.

#### Execution

To reproduce the manuscript results, the script `calculate_ppVAF_posteriors.py` should be executed on the WGS and WES mafs for each patient in the dataset to produce the patient-specific unnormalized ppVAF posterior probability matrices (saved as .npy files) and trimmed mafs. These matrices are large and not provided in the Zenodo data distribution. These matrices and mafs are required inputs to the script `process_unnormalized_ppVAFs.ipynb`, which annotates the mafs with the final ppVAF values and generate the combined and filtered maf files in the Zenodo data distribution (`combined_noshared_FILTERED_muts_WGS.maf` and `combined_noshared_FILTERED_muts_WES.maf`). The expected clonal SNV counts and poly/monoclonal calls for each sample are also calculated from these matrices and stored as `clonal_noshared_WES_WGS_polycalls.csv`. These mafs and poly/mono calls are used to generate most of the downstream figures in the manuscript. This notebook also generates the purity-marginalized posterior probability matrices used as inputs to generate the joyplots in Fig. 1f.

The notebook `scATACseq_sample_purity.ipynb` takes the published scATAC-seq data from our FAP cohort and creates smoothed epithelial cell fraction distributions that are used to marginalize the ppVAF posterior probabilities over the estimated purity distributions. The inputs and outputs to this script are distributed in this repo in `data/scATACseq_annotations`.

#### Requirements

This was developed using `Python 3.9`.

Python modules required:
```
pandas
numpy
seaborn
matplotlib
scipy
```