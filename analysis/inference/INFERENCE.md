# Inferences

The inference procedure can be ran using the `exe_abc.py` script within this directory.

# Summary statistics extraction

We first use `extract_summary_stats.py` to generate and save summary stats for each sample in our MAF file from the `purity_CCF` column.

```
python extract_summary_stats.py
```

This generates the outputs from the data used in the inference method.

```bash
# From within the repos parent directory:
(pyabc) [rschenck@smsh11dsu-srcf-d15-37 model]$ python ./analysis/inference/exe_abc.py --help
usage: exe_abc.py [-h] [--npoints NPOINTS] [--epsilon_diff EPSILON_DIFF] [--max_nr_pop MAXITERS] [--sample_name SAMPLE_NAME] [--patient_id PATIENT] [--wxs WXS]
                        [--outputdir OUTPUTDIR] [--data_sum_stats D_SUM_STATS] [--continue_inf CONTINUE_INF]

Run inference.

options:
  -h, --help            show this help message and exit
  --npoints NPOINTS     Number of simulations per iteration (default:250)
  --epsilon_diff EPSILON_DIFF
                        Difference in between runs (default:0.000000001)
  --max_nr_pop MAXITERS
                        Maximum number of inference iterations (default:10)
  --sample_name SAMPLE_NAME
                        Maximum number of inference iterations (default:10)
  --patient_id PATIENT  Output path for results
  --wxs Denote whether this is WGS or WES (for appropriate labeling only)
  --outputdir OUTPUTDIR
                        Output path for results, for a group of inferences this should be the top level directory.
  --data_sum_stats D_SUM_STATS
                        Path to data summary statistics
  --continue_inf CONTINUE_INF
                        Continue inference already started.
```

For inferences, the constant model parameters are placed in the `./analysis/inference/util/conf.yaml` file.
