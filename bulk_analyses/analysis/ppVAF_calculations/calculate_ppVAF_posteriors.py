import pandas as pd
import numpy as np
from scipy.stats import binom
import sys
from ppVAF_utils import *

'''
Estimates ppVAF posterior probabilities given sample purity values. Takes in a tab-separated maf file with a "Patient" column and computes the unnormalized ppVAF posterior values for all mutations in that patient, specified by the first input argument when running the script. Requires "tcn" and "lcn" columns, which are the total copy number and minor copy number at each mutation locus. Assumes mutations are on the major allele. Saves a .npy matrix with dimensions 100 x 1000 x [[number of mutations]] dimensions, where the first dimension is the purity value (0.01 to 1, evenly spaced), and the second dimension is the ppVAF value (0.001 to 1 evenly spaced). Also saves a patient specific mutation file as a csv with the same columns as in the original maf, with mutations in the same order as in the 3rd dimension of the .npy matrix.

SCRIPT USAGE:
python3 calculate_ppVAF_posteriors.py [[PATIENT ID]] [[INPUT MAF FILEPATH]] [[OUTPUT NUMPY MATRIX FILEPATH]] [[OUTPUT MUTATION TABLE FILEPATH]]
'''

sample = sys.argv[1]
maf_in = sys.argv[2]
numpy_out = sys.argv[3]
maf_out = sys.argv[4]

FAP_maf = pd.read_csv(maf_in, sep="\t")

if sample != "none":
    to_ccf = FAP_maf[FAP_maf["Patient"] == sample]
else:
    to_ccf = FAP_maf

to_return = np.zeros((100,1000,len(to_ccf)))

for i in range(len(to_ccf)):
    row = to_ccf.iloc[i]
    total_cn = row["tcn"]
    minor_cn = row["lcn"]
    if np.isnan(minor_cn):
        minor_cn = max(total_cn - 1, 0)
    mut_cn = total_cn - minor_cn
    to_return[:, :, i] = estimate_ccf_purity(total_cn, mut_cn, row["t_alt_count"], row["t_depth"])
    
np.save(numpy_out, to_return)
to_ccf.to_csv(maf_out, index=False)
