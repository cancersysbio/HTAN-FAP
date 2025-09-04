### Code repo for: *Polyclonal origins of human premalignant colorectal lesions*

### :file_folder: Repository content ###
This GitHub repository contains the following folders:
```
| 
└─ bulk_analyses/: directory containing code to analyze bulk sequencing data (WES/WGS)
|   └─ analysis/: code used for data analysis
|   └─ data/: summary level data to generate main, extended and SI figures
|   └─ figures/: scripts, notebooks and other code to reproduce figures
|       └─ main/: scripts and notebooks for main figures
|       └─ extended_data/: scripts and notebooks for extended data figures
|       └─ supplementary/: scripts and notebooks for supplementary figures
| 
└─ sgWGS_analyses/: directory containing single-crypt whole genome sequencing analyses
    └─ sgutils/: utility module for filtering and phylogenetic analysis
    └─ data/: sgWGS-specific data files
    └─ binomial_mixture_filter.ipynb: filtering and analysis notebook
    └─ generate_figures.ipynb: main figure generation notebook
    └─ supplemental_figures.ipynb: supplementary figure generation notebook
```

### :file_folder: Access to data
- DNA sequencing data and metadata have been deposited at the HTAN portal: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002371.v3.p1.

- The raw single-crypt WGS data have been deposited in NGDC under the accession number PRJCA023981 (https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA023981).
  
- WES data for the FAP multi-region and sporadic CRC cohorts should be requested through the original authors/organizations. We have not provided raw data or processed data files but have provided the code needed to generate our results from the raw data. Note that throughout this repo in READMEs and filenames, the external FAP patient multiregion WES data from Li et al (https://pubmed.ncbi.nlm.nih.gov/31744909/) is referred to as **PUTH** and the sporadic multiregion WES data from Cross et al (https://pubmed.ncbi.nlm.nih.gov/30177804/) is referred to as **SCORT**.

- Additional larger intermediate data files that may be helpful to run the code here and reproduce our findings (e.g., somatic MAF files describing the mutations in each sample) can be downloaded from Zenodo at https://doi.org/10.5281/zenodo.13228021.

### :white_check_mark: Citation

