import pandas as pd
import numpy as np
from scipy.stats import binom
import sys

'''
Estimates ppVAF posterior probabilities given sample purity values. Takes in a tab-separated maf file with a "Patient" column and computes the unnormalized ppVAF posterior values for all mutations in that patient, specified by the first input argument when running the script. Requires "tcn.sequenza" column, which is the total copy number at each mutation locus. Estimates the expected number of mutant copies as in FACETS-SUITE. Saves a .npy matrix with dimensions 100 x 1000 x [[number of mutations]] dimensions, where the first dimension is the purity value (0.01 to 1, evenly spaced), and the second dimension is the ppVAF value (0.001 to 1 evenly spaced). Also saves a patient specific mutation file as a csv with the same columns as in the original maf, with mutations in the same order as in the 3rd dimension of the .npy matrix.

SCRIPT USAGE:
python3 calculate_ppVAF_posteriors.py [[PATIENT ID]] [[INPUT MAF FILEPATH]] [[OUTPUT NUMPY MATRIX FILEPATH]] [[OUTPUT MUTATION TABLE FILEPATH]]
'''

sample = sys.argv[1]
maf_in = sys.argv[2]
numpy_out = sys.argv[3]
maf_out = sys.argv[4]

def expected_mutant_copies(t_var_freq, total_copies, purity):
    '''
    From FACETS-SUITE (https://github.com/mskcc/facets-suite/tree/master), translated from R from ccf-annotate-maf.R.
    '''
    if np.isnan(total_copies) or total_copies == 0:
        return np.nan
    else:
        mu = t_var_freq * (1 / purity) * (purity * total_copies + (1 - purity) * 2)
        if mu < 1:
            return 1
        else:
            return round(abs(mu))

def estimate_ccf_purity(total_copies, t_alt_count, t_depth):
    '''
    Estimates posterior probability of ppVAF values given sequencing data at a single mutation locus across a mesh sweeping across sample purity values from 0.01 to 1.
    params:
    total_copies: int total copy number at mutation locus
    t_alt_count: int number of mutant reads at locus
    t_depth: int total number of reads at locus
    
    returns:
    100 x 1000 element numpy matrix of unnormalized ppVAF posterior probabilties at the mutation locus, where rows are sample purities from 0.01 to 1 and columns are ppVAF values from 0.001 to 1 (both evenly spaced)
    
    '''
    ccfs = np.linspace(0.001, 1, 1000).reshape(1, -1)
    purities = np.linspace(0.01, 1, 100)
    
    mutant_copies = np.array([expected_mutant_copies(t_alt_count/t_depth, total_copies, purity) for purity in purities]).reshape(-1, 1)
    
    purities = purities.reshape(-1, 1)
    
    expected_vafs = purities * mutant_copies * ccfs / (2 * (1 - purities) + purities * total_copies)
    return(binom.pmf(t_alt_count, t_depth, expected_vafs))

FAP_maf = pd.read_csv(maf_in, sep="\t")

to_ccf = FAP_maf[FAP_maf["Patient"] == sample]
to_return = np.zeros((100,1000,len(to_ccf)))

for i in range(len(to_ccf)):
    row = to_ccf.iloc[i]
    total_cn = row["tcn.sequenza"]
    if np.isnan(total_cn) or total_cn == 0:
        total_cn = 2
    to_return[:, :, i] = estimate_ccf_purity(total_cn, row["t_alt_count"], row["t_depth"])
    
np.save(numpy_out, to_return)
to_ccf.to_csv(maf_out, index=False)
