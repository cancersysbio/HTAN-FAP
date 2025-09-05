# sgWGS Analyses

This subfolder contains analyses and tools for processing and analyzing
single-gland whole genome sequencing (sgWGS) data from FAP (Familial Adenomatous
Polyposis) samples, with a focus on phylogenetic reconstruction.

## Key Components

### Utility Module (`sgutils/`)

#### `clone_detection.py`

#### `phylo_utils.py`

### Analysis Notebooks

- **`binomial_mixture_filter.ipynb`**: Implements clonal peak detection using truncated binomial mixture models with configurable VAF thresholds
- **`generate_figures.ipynb`**: Main analysis pipeline for generating publication figures, phylogenetic trees, and statistical summaries
- **`supplemental_figures.ipynb`**: Additional analyses and figures for supplementary materials

### Data Organization

**NOTE:** The contents of the `data/` folder are downloadable from Zenodo at https://doi.org/10.5281/zenodo.13228021, by extracting the `sgWGS.zip` archive.

- **`data/`**: Contains input mutation data, metadata, and reference files (*from Zenodo*)
  - **`v3/`**: Latest version of somatic mutation data (MAF format)
  - **`apc_depth/`**: APC gene depth coverage data for quality control
  - **`TableS1_sgWGS_mutMeta.csv`**: Sample metadata including TMB, driver mutations, and sequencing metrics
  - **`PanCanDrivers_COADREAD_Cell2018.csv`**: Reference driver gene list for colorectal cancer
  - **`sgWGS_black_list.txt`**: Samples excluded from analysis

- **`output/`**: Generated analysis results and model outputs

- **`figures/`**: Publication figures and visualizations

## Usage

To reproduce the sgWGS figures, you will first have to run the binomial mixture
filter notebook, so that subsequent notebooks have access to the information of
which samples passed the filter:

### Run clonal detection analysis

``` bash
jupyter nbconvert --to notebook --execute binomial_mixture_filter.ipynb --output binomial_mixture_filter_output.ipynb
```

Once this has completed, you will be able to reproduce all sgWGS figures and analyses presented in the paper:

### Generate main figures
``` bash
jupyter nbconvert --to notebook --execute generate_figures.ipynb --output generate_figures_output.ipynb
```

### Generate supplementary figures
``` bash
jupyter nbconvert --to notebook --execute supplemental_figures.ipynb --output supplemental_figures_output.ipynb
```
