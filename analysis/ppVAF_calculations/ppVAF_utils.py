import numpy as np
import pandas as pd
from scipy.stats import binom

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

def estimate_ccf_purity(total_copies, mutant_copies, t_alt_count, t_depth):
    '''
    Estimates posterior probability of ppVAF values given sequencing data at a single mutation locus across a mesh sweeping across sample purity values from 0.01 to 1.
    params:
    total_copies: int total copy number at mutation locus
    mutant_copies: int mutant copy number at mutation locus
    t_alt_count: int number of mutant reads at locus
    t_depth: int total number of reads at locus
    
    returns:
    100 x 1000 element numpy matrix of unnormalized ppVAF posterior probabilties at the mutation locus, where rows are sample purities from 0.01 to 1 and columns are ppVAF values from 0.001 to 1 (both evenly spaced)
    
    '''
    ccfs = np.linspace(0.001, 1, 1000).reshape(1, -1)
    purities = np.linspace(0.01, 1, 100)
    
    #mutant_copies = np.array([expected_mutant_copies(t_alt_count/t_depth, total_copies, purity) for purity in purities]).reshape(-1, 1)
    
    purities = purities.reshape(-1, 1)
    
    expected_vafs = purities * mutant_copies * ccfs / (2 * (1 - purities) + purities * total_copies)
    return(binom.pmf(t_alt_count, t_depth, expected_vafs))

def purity_dist_normalize(prob_mat, maf, purity_dict):
    all_stages = list(set(maf["Stage"]))
    for stage in all_stages:
        if stage not in purity_dict:
            raise AssertionError("invalid stage")
        has_stage = np.nonzero((maf["Stage"] == stage).tolist())[0]
        prob_mat[:, :, has_stage] = np.multiply(prob_mat[:, :, has_stage], purity_dict[stage].reshape((-1, 1, 1)))
    prob_mat = np.divide(prob_mat, np.sum(prob_mat, axis=(0,1)).reshape(1, 1, -1))
    return prob_mat

def get_ccfs_clonality_many(probs):
    num_ccf_grid = np.shape(probs)[0]
    
    ccfs = np.argmax(probs, axis=0)
    ccf_half_max = probs > (np.max(probs, axis=0).reshape((1, -1)) / 2)
    ccf_half_max = np.where(ccf_half_max==0, np.nan, np.arange(num_ccf_grid).reshape(-1,1))
    
    ccf_lower = np.maximum(np.nanmin(ccf_half_max, axis=0) - 1, 1) # closest ccf value before half-max range (within 0-1 range)
    ccf_upper = np.minimum(np.nanmax(ccf_half_max, axis=0) + 1, num_ccf_grid) # closest ccf value after half-max range (within 0-1 range)

    ccf_lower = ccf_lower / num_ccf_grid
    ccf_upper = ccf_upper / num_ccf_grid
    
    ccfs = ccfs/num_ccf_grid
    
    return (ccfs, ccf_lower, ccf_upper)

def get_CCF_MAP(prob_mat, maf, ccf_col="ppVAF", bounds_prefix="ppVAF"):
    #adds new columns to maf df with best ccf estimate marginalized over purity distribution
    marginalized = np.sum(prob_mat, axis=0) * prob_mat.shape[1]
    CCFs, lower, upper = get_ccfs_clonality_many(marginalized)
    maf[ccf_col] = CCFs
    maf[bounds_prefix+"_lower"] = lower
    maf[bounds_prefix+"_upper"] = upper
    return marginalized

def add_ppVAFs(prob_mat, maf, purity_dict, clonal_thresholds):
    maf_save = None
    
    num_ccf_grid = np.shape(prob_mat)[1]
    all_marg = np.zeros((np.shape(prob_mat)[2], np.shape(prob_mat)[1]))
    start_idx = 0
    for sample in list(set(maf["Tumor_Sample_Barcode"])):
        is_sample = np.nonzero((maf["Tumor_Sample_Barcode"] == sample).tolist())[0]
        new_maf = maf.iloc[is_sample]
        new_mat = prob_mat[:, :, is_sample]
        n_muts = len(new_maf)
        new_mat = purity_dist_normalize(new_mat, new_maf, purity_dict)
        marg = get_CCF_MAP(new_mat, new_maf)
        all_marg[start_idx:start_idx+n_muts,:] = np.transpose(marg)
        start_idx += n_muts
        
        for clonal_threshold in clonal_thresholds:
            idx_threshold = int(clonal_threshold*num_ccf_grid)
            probs_clonal = np.sum(new_mat[:, idx_threshold:, :], axis=(0, 1))
            new_maf["clonal_cont_"+str(clonal_threshold)] = probs_clonal
        
        if maf_save is None:
            maf_save = new_maf
        else:
            maf_save = pd.concat([maf_save, new_maf], ignore_index=True)
    return maf_save, all_marg

def expected_count_clonal(prob_mat, maf, clonal_thresholds=[0.95], filter_maf=None):
    num_ccf_grid = np.shape(prob_mat)[1]
    
    if filter_maf is not None:
        filtered_mat = prob_mat[:, :, filter_maf]
        CI_clonal = [np.nansum(maf.iloc[filter_maf]["ppVAF_upper"]==1)]
    else:
        filtered_mat = prob_mat
        CI_clonal = [np.nansum(maf["ppVAF_upper"]==1)]
    
    
    df_return = pd.DataFrame({"CI_clonal":CI_clonal})
    
    
    for clonal_threshold in clonal_thresholds:
        idx_threshold = int(clonal_threshold*num_ccf_grid)
        probs_clonal = np.sum(prob_mat[:, idx_threshold:, :], axis=(0, 1))
        maf["clonal_cont_"+str(clonal_threshold)] = probs_clonal
        df_return["exp_clonal_"+str(clonal_threshold)] = np.nansum(filtered_mat[:, idx_threshold:, :])
            
    return df_return

def add_ccfs_count_clonal(prob_mat, maf, purity_dict, clonal_thresholds):
    maf_save = None
    clonal = None
    all_marg = np.zeros((np.shape(prob_mat)[2], np.shape(prob_mat)[1]))
    start_idx = 0
    for sample in list(set(maf["Tumor_Sample_Barcode"])):
        is_sample = np.nonzero((maf["Tumor_Sample_Barcode"] == sample).tolist())[0]
        new_maf = maf.iloc[is_sample]
        new_mat = prob_mat[:, :, is_sample]
        n_muts = len(new_maf)
        new_mat = purity_dist_normalize(new_mat, new_maf, purity_dict)
        marg = get_CCF_MAP(new_mat, new_maf)
        all_marg[start_idx:start_idx+n_muts,:] = np.transpose(marg)
        start_idx += n_muts
        
        filter_maf = new_maf["Variant_Type"]=="SNP"
        clonal_add = expected_count_clonal(new_mat, new_maf, clonal_thresholds=clonal_thresholds, filter_maf=np.nonzero((filter_maf).tolist())[0])
        clonal_add["sample"] = sample

        if clonal is None:
            clonal = clonal_add
        else:
            clonal = pd.concat([clonal, clonal_add])
        if maf_save is None:
            maf_save = new_maf
        else:
            maf_save = pd.concat([maf_save, new_maf], ignore_index=True)
    return maf_save, clonal, all_marg

def sample_to_patient(sample):
    if sample[0] == "A":
        patient = sample[:4]
    else:
        patient = sample[:1] + "001"
    return patient

def make_ID(components):
    return '.'.join([str(x) for x in components])

def extract_full_marginal(maf, marg, gene_symbol):
    only_genes_idx = maf["Hugo_Symbol"] == gene_symbol
    mut_names = maf[only_genes_idx][['Tumor_Sample_Barcode', 'Mut_ID']].agg(make_ID, axis=1).tolist()
    marg_dist = marg[only_genes_idx, :].transpose()
    return pd.DataFrame(marg_dist, columns=mut_names)

def optimize_cutoff(posterior_colname, true_clonal, true_subclonal, fraction_clonal, return_diagnostic=False):
    cutoffs = np.linspace(0, 1, 101)
    
    subclonal_accuracies = []
    clonal_accuracies = []
    for cutoff in cutoffs:
        subclonal_accuracies.append(np.mean(true_subclonal[posterior_colname] < cutoff))
        clonal_accuracies.append(np.mean(true_clonal[posterior_colname] > cutoff))
    
    subclonal_accuracies = np.array(subclonal_accuracies)

    mean_accuracies = fraction_clonal*np.array(clonal_accuracies) + (1-fraction_clonal)*subclonal_accuracies
    expected_error = (1-fraction_clonal)*(1-subclonal_accuracies) - fraction_clonal*(1-np.array(clonal_accuracies))
    
    max_index = np.argmin(np.abs(expected_error))
    #max_index = np.argmax(mean_accuracies)
    
    if return_diagnostic:
        return cutoffs[max_index], expected_error, subclonal_accuracies, clonal_accuracies, mean_accuracies
    else:
        return cutoffs[max_index]

def add_count_clonal(maf, posterior_colname, true_clonal, true_subclonal, poly_cutoff, tolerance=1):
    sample_to_stage = dict(zip(maf["Tumor_Sample_Barcode"], maf["Stage"]))
    
    initial_cutoff = optimize_cutoff(posterior_colname, true_clonal, true_subclonal, 0.5)
    
    maf["init_clonal"] = maf[posterior_colname] > initial_cutoff
    maf["total_muts"] = 1

    to_return = maf[["init_clonal", "Tumor_Sample_Barcode", "total_muts"]].groupby("Tumor_Sample_Barcode").sum()
    to_return["stage"] = [sample_to_stage[x] for x in to_return.index]
      
    clonal_counts = []
    annotated_mafs = []
    for sample in to_return.index:
        old_clonal = to_return.loc[sample]["init_clonal"]
        total_muts = to_return.loc[sample]["total_muts"]
        only_sample = maf[maf["Tumor_Sample_Barcode"]==sample]
        curr_diff = np.inf
        
        while curr_diff > tolerance:
            optimal_cutoff = optimize_cutoff(posterior_colname, true_clonal, true_subclonal, old_clonal/total_muts)
            new_clonal = np.sum(only_sample[posterior_colname] > optimal_cutoff)
            curr_diff = np.abs(old_clonal-new_clonal)
            old_clonal = new_clonal
        
        only_sample["final_clonal"] = only_sample[posterior_colname] > optimal_cutoff
        clonal_counts.append(np.sum(only_sample["final_clonal"]))
        only_sample.drop(columns=["total_muts"])
        annotated_mafs.append(only_sample)
    to_return["final_clonal"] = clonal_counts
    to_return["is_poly"] = to_return["final_clonal"] <= poly_cutoff
    
    return to_return, pd.concat(annotated_mafs, ignore_index=True)

def single_diploid_mutation_posterior(t_alt, t_depth, purity):
    prob_mat = estimate_ccf_purity(2, 1, t_alt, t_depth)

    purity_dist = np.zeros(100)
    purity_dist[int(purity*100)] = 1

    normed_prob_mat = np.multiply(prob_mat, purity_dist.reshape((-1, 1)))
    normed_prob_mat = np.divide(normed_prob_mat, np.sum(normed_prob_mat, axis=None))
    marginalized = np.sum(normed_prob_mat, axis=0) * 1000
    return marginalized