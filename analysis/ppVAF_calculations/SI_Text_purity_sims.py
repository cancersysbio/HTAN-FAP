import pandas as pd
import numpy as np
from scipy.stats import binom
import sys, pickle, os
from ppVAF_utils import *

# root directory where processed mutation data are stored
# if using the zenodo data, this should point to the location of the base zenodo directory (including wes, wgs, etc subfolders)
data_dir = ""
filtered_maf_WGS = pd.read_csv(data_dir+"HTAN_WGS_filtered_ppVAFs.maf", sep="\t")

# filepath of the github repo data directory
repo_data_dir = ""

only_SNVs_WGS = filtered_maf_WGS[filtered_maf_WGS["Variant_Type"]=="SNP"]
true_subclonal_WGS = only_SNVs_WGS[only_SNVs_WGS["Stage"]=="Mucosa"]

clonal_count_dir = repo_data_dir+"clonal_count_estimation/"
expected_clonal = pd.read_csv(clonal_count_dir+"clonal_SNVs_WES_WGS.csv")
simulated_clonal_WGS = pd.read_csv(clonal_count_dir+"simulated_clonal_WGS.csv")

output_dir = repo_data_dir+"purity_sims/"
os.makedirs(output_dir, exist_ok=True)

def lower_purity(maf, new_purity):
    maf["t_ref_count"] = np.round(maf["t_ref_count"] * new_purity + maf["t_depth"] * (1 - new_purity))
    maf["t_alt_count"] = maf["t_depth"] - maf["t_ref_count"]

def simulate_mixture(sample1, sample2, comb_maf, frac1, vaf_threshold=0.01, alt_threshold=2):
    maf1 = comb_maf[comb_maf["Tumor_Sample_Barcode"]==sample1]
    maf2 = comb_maf[comb_maf["Tumor_Sample_Barcode"]==sample2]
    
    duplicates = list(set(maf1["Mut_ID"]).intersection(maf2["Mut_ID"]))
    
    concat_maf1 = maf1[~np.isin(maf1["Mut_ID"], duplicates)]
    concat_maf2 = maf2[~np.isin(maf2["Mut_ID"], duplicates)]
    
    lower_purity(concat_maf1, frac1)
    lower_purity(concat_maf2, 1-frac1)
    dups = maf1[np.isin(maf1["Mut_ID"], duplicates)]
    
    new_refs = []
    new_alts = []
    for mut_id in dups["Mut_ID"]:
        new_refs.append(round(frac1*maf1[maf1["Mut_ID"]==mut_id].iloc[0]["t_ref_count"] + (1-frac1)*maf2[maf2["Mut_ID"]==mut_id].iloc[0]["t_ref_count"]))
        new_alts.append(round(frac1*maf1[maf1["Mut_ID"]==mut_id].iloc[0]["t_alt_count"] + (1-frac1)*maf2[maf2["Mut_ID"]==mut_id].iloc[0]["t_alt_count"]))
    dups["t_ref_count"] = new_refs
    dups["t_alt_count"] = new_alts
    dups["t_depth"] = dups["t_ref_count"] + dups["t_alt_count"]
    
    to_return = pd.concat([concat_maf1, concat_maf2, dups], ignore_index=True)
    to_return["vaf"] = to_return["t_alt_count"]/to_return["t_depth"]
    to_return = to_return[to_return["vaf"] > vaf_threshold]
    to_return = to_return[to_return["t_alt_count"] > alt_threshold]
    return to_return

def get_clonal_pipeline(to_ccf, thresholds, purity_dict):
    to_return = np.zeros((100,1000,len(to_ccf)))

    for i in range(len(to_ccf)):
        row = to_ccf.iloc[i]
        total_cn = row["tcn"]
        minor_cn = row["lcn"]
        if np.isnan(minor_cn):
            minor_cn = max(total_cn - 1, 0)
        mut_cn = total_cn - minor_cn
        to_return[:, :, i] = estimate_ccf_purity(total_cn, mut_cn, row["t_alt_count"], row["t_depth"])
    maf_save, marg = add_ppVAFs(to_return, to_ccf, purity_dict, thresholds)
    poly_calls_WGS, _ = add_count_clonal(maf_save, "clonal_cont_0.8", simulated_clonal_WGS, true_subclonal_WGS, 63)
    return poly_calls_WGS

annot_dir = repo_data_dir+"scATACseq_annotations/"
purity_dict = pickle.load(open(annot_dir+"scATAC_purities.p", "rb"))

mix_frac = float(sys.argv[1])

to_sample = expected_clonal[~expected_clonal["is_poly"]]
to_sample = to_sample[to_sample["has_WGS"]]
to_sample = to_sample[np.isin(to_sample["stage"], ["Dysplasia", "Benign"])]["Tumor_Sample_Barcode"].tolist()

normal_sample = expected_clonal[expected_clonal["has_WGS"]]
normal_sample = normal_sample[normal_sample["stage"]=="Mucosa"]["Tumor_Sample_Barcode"].tolist()

thresholds = [0.6, 0.7, 0.8, 0.9, 0.95]

n_iters = 200
fracs = []
clonal_muts = []
APC_muts = []
for i in range(n_iters):
    sampled = np.random.choice(to_sample)
    normal_sampled = np.random.choice(normal_sample)
    to_ccf = simulate_mixture(sampled, normal_sampled, only_SNVs_WGS, mix_frac)
    print(len(to_ccf))
    APC_muts.append(to_ccf[np.logical_and(to_ccf["Hugo_Symbol"]=="APC", to_ccf["Driver"])]["vaf"].tolist())
    to_ccf["Tumor_Sample_Barcode"] = "mixed"
    to_add = get_clonal_pipeline(to_ccf, thresholds, purity_dict)
    to_add["starting_sample"] = sampled
    clonal_muts.append(to_add)
    print("--")
mixing_sims = pd.concat(clonal_muts, ignore_index=True)
mixing_sims["n_APC"] = [len(x) for x in APC_muts]
mixing_sims["mix_frac"] = mix_frac
mixing_sims.to_csv(output_dir+"/frac_"+str(mix_frac)+".csv", index=False)