import pandas as pd
import numpy as np
from scipy.stats import binom
import sys

'''
Estimates ppVAF posterior probabilities given sample purity values. Takes in a tab-separated maf file with a "Patient" column and computes the unnormalized ppVAF posterior values for all mutations in that patient, specified by the first input argument when running the script. Requires "tcn.sequenza" and "ln.sequenza" columns, which are the total copy number and minor allele copy number at each mutation locus, respectively. Assumes mutation is on the major allele. Saves a .npy matrix with dimensions 100 x 1000 x [[number of mutations]] dimensions, where the first dimension is the purity value (0.01 to 1, evenly spaced), and the second dimension is the ppVAF value (0.001 to 1 evenly spaced). Also saves a patient specific mutation file as a csv with the same columns as in the original maf, with mutations in the same order as in the 3rd dimension of the .npy matrix.

SCRIPT USAGE:
python3 calculate_ppVAF_posteriors.py [[PATIENT ID]] [[INPUT MAF FILEPATH]] [[OUTPUT NUMPY MATRIX FILEPATH]] [[OUTPUT MUTATION TABLE FILEPATH]]
'''

sample = sys.argv[1]
maf_in = sys.argv[2]
numpy_out = sys.argv[3]
maf_out = sys.argv[4]

def estimate_ccf_purity(total_copies, mutant_copies, t_alt_count, t_depth):
    '''
    Estimates posterior probability of ppVAF values given sequencing data at a single mutation locus across a mesh sweeping across sample purity values from 0.01 to 1.
    params:
    total_copies: int total copy number at mutation locus
    mutant_copies: int minor allele copy number at mutation locus
    t_alt_count: int number of mutant reads at locus
    t_depth: int total number of reads at locus
    
    returns:
    100 x 1000 element numpy matrix of unnormalized ppVAF posterior probabilties at the mutation locus, where rows are sample purities from 0.01 to 1 and columns are ppVAF values from 0.001 to 1 (both evenly spaced)
    '''
    ccfs = np.linspace(0.001, 1, 1000).reshape(1, -1)
    purities = np.linspace(0.01, 1, 100).reshape(-1, 1)
    expected_vafs = purities * ccfs * mutant_copies / (2 * (1 - purities) + purities * total_copies)
    return(binom.pmf(t_alt_count, t_depth, expected_vafs))

FAP_maf = pd.read_csv(maf_in, sep="\t")

to_ccf = FAP_maf[FAP_maf["Patient"] == sample]
to_return = np.zeros((100,1000,len(to_ccf)))

for i in range(len(to_ccf)):
    row = to_ccf.iloc[i]
    total_cn = row["tcn.sequenza"]
    if np.isnan(total_cn):
        total_cn = 2
    mut_cn = row["tcn.sequenza"] - row["ln.sequenza"]
    if np.isnan(mut_cn):
        mut_cn = 1
    to_return[:, :, i] = estimate_ccf_purity(total_cn, mut_cn, row["t_alt_count"], row["t_depth"])
    
np.save(numpy_out, to_return)
to_ccf.to_csv(maf_out, index=False)
