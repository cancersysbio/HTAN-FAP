# Must use h5py=3.2.0
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import os, sys
import warnings
# import pyabc
from scipy.stats import gaussian_kde
from statsmodels.distributions.empirical_distribution import ECDF

warnings.filterwarnings("ignore")
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

sys.path.append('./util')
from util.IO import read_yaml_to_dict, get_parameter_output, run_model_job, h5_tree, read_h5, process_data

def calculate_densities(arr, x_range = np.linspace(0,1,21)):
    # Define the range
    x_range = np.linspace(0, 1, 50)

    # Calculate the KDE and density
    kde = gaussian_kde(arr)
    density = kde(x_range)

    return density

def calculate_ecdf(arr, points = np.linspace(0,1,21)):
    ecdf = ECDF(arr)
    
    # Eval ECDF at points
    ecdf_values = ecdf(points)
    
    return ecdf_values

def get_sample_summary_data(muts):
    # Take only autosomes
    # Return num total mutations, num subclonal mutations, num clonal mutations, array of mutation frequencies as vaf and ccfs
    muts = muts[muts['Chromosome'].isin(['chrY', 'chrX'])==False].reset_index(drop=True)

    # Filter out all Zeros for purity_ccf
    muts = muts[muts['purity_ccf']>0.].reset_index(drop=True)

    muts['clonality'] = muts['purity_clonal'].str.replace('+','').str.replace('-','')
    
    c = muts.groupby('clonality')['purity_clonal'].size().reset_index(drop=False)
    
    num_clonal_subclonal_total = np.zeros(3)
    clonal = c[c['clonality']=='CLONAL']['purity_clonal'].to_list()
    subclonal = c[c['clonality']=='SUBCLONAL']['purity_clonal'].to_list()
    if len(clonal)>0:
        num_clonal_subclonal_total[0] = clonal[0]
    if len(subclonal)>0:
        num_clonal_subclonal_total[1] = subclonal[0]

    num_clonal_subclonal_total[2] = muts.shape[0]

    ccfs = muts['purity_ccf'].to_numpy()
    vafs = muts['vaf'].to_numpy()

    densities = calculate_densities(ccfs)

    probs = calculate_ecdf(ccfs)

    densities_vaf = calculate_densities(vafs)

    return(num_clonal_subclonal_total, ccfs, vafs, probs, densities, densities_vaf)

# data = "/home/rschenck/oak_dir/cluster/data/combined_noshared_FILTERED_muts_WES.csv"
data = "./data/combined_noshared_FILTERED_muts_WGS.csv"
ages = "./data/FAP_age_ethnicity.tsv" # Contains the colectomy age
wxs = 'wgs'

data = pd.read_csv(data, sep=',', low_memory=False)
ages = pd.read_csv(ages, sep='\t', low_memory=False)

muts = data[data['Chromosome'].isin(['chrY', 'chrX'])==False].reset_index(drop=True)

for s, sdf in muts.groupby(['Patient','Tumor_Sample_Barcode', 'Stage']):
    p, s, stage = s[0], s[1], s[2]
    # p = p.replace('XX','F').replace('XX','G')
    # s = s.replace('XX','F').replace('XX','G')
    pat_age = ages[ages['patient_id'] == p]["Colectomy_age"].iloc[0]
    median_depth = sdf['t_depth'].median()
    with h5py.File(f"/home/rschenck/oak_dir/cluster/data/summary_stats/{wxs}/{s}.hdf5", "w") as f:
        if sdf.shape[0] > 1:
            summ_stats = get_sample_summary_data(sdf)
            # Output as hdf5
            counts = f.create_dataset(f"{s}/num_clonal_subclonal_total", data=summ_stats[0])
            ccfs = f.create_dataset(f"{s}/ccfs", data=summ_stats[1])
            vafs = f.create_dataset(f"{s}/vafs", data=summ_stats[2])
            probs = f.create_dataset(f"{s}/probs", data=summ_stats[3])
            densities = f.create_dataset(f"{s}/density", data=summ_stats[4])
            densities_vaf = f.create_dataset(f"{s}/density_vaf", data=summ_stats[5])
            age = f.create_dataset(f"{s}/colectomy_age", data=np.asarray([pat_age]))
            stage = f.create_dataset(f"{s}/stage", data=stage)
            purity = f.create_dataset(f"{s}/purity", data=np.asarray([0.8]))
            depth = f.create_dataset(f"{s}/depth", data=np.asarray([median_depth]))
        else:
            print(f"Not enough mutations: {s}")
