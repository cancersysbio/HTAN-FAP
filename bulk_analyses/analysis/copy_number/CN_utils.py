import pandas as pd
import numpy as np

def deduplicate_gene_calls(sample_gene_CNs):
    sample_gene_CNs["length"] = sample_gene_CNs["gene_end"] - sample_gene_CNs["gene_start"]
    duplicates = sample_gene_CNs.value_counts("gene_name")
    duplicates = duplicates[duplicates > 1]
    
    if len(duplicates) == 0:
        to_return = sample_gene_CNs
    
    else:
        deduplicated = []
        for gene in duplicates.index:
            both_listings = sample_gene_CNs[sample_gene_CNs["gene_name"]==gene]
            to_add = pd.DataFrame(both_listings.iloc[0]).transpose()
            to_add["gene_start"] = np.min(both_listings["gene_start"])
            to_add["gene_end"] = np.max(both_listings["gene_end"])
            longest_segment = both_listings.sort_values("length", ascending=False).iloc[0]
            to_add["tcn_em"] = longest_segment["tcn_em"]
            to_add["lcn_em"] = longest_segment["lcn_em"]
            deduplicated.append(to_add)

        deduplicated = pd.concat(deduplicated)
        to_return = sample_gene_CNs[~np.isin(sample_gene_CNs["gene_name"], duplicates.index)]
        to_return = pd.concat([to_return, deduplicated])
    return to_return[["gene_name", "gene_chrom", "gene_start", "gene_end", "tcn_em", "lcn_em"]]

def get_CN_type(chrom, tcn, lcn, is_male, is_doubled):
    assert chrom != "chrY", "no chrY allowed"
    
    normal_tcn = 2
    if chrom == "chrX":
        if is_male:
            normal_tcn = 1
    
    if is_doubled:
        normal_tcn *= 2
    
    tcn_call = "NORMAL"
    if tcn != normal_tcn:
        if tcn > normal_tcn:
            tcn_call = "GAIN"
        else:
            tcn_call = "LOSS"
            
    het_call = "NO LOH"
    if lcn == 0:
        het_call = "LOH"
    return tcn_call, het_call

def add_CN_types(sample_genes_CNs, is_male, is_doubled):
    # ADDS columns "CN_call" and "LOH_call" IN PLACE- does not return anything
    cn_types = []
    loh_types = []
    for i in range(len(sample_genes_CNs)):
        row = sample_genes_CNs.iloc[i]
        call, het = get_CN_type(row["gene_chrom"], row["tcn_em"], row["lcn_em"], is_male, is_doubled)
        cn_types.append(call)
        loh_types.append(het)
    sample_genes_CNs["CN_call"] = cn_types
    sample_genes_CNs["LOH_call"] = loh_types

def check_if_male(patient):
    if patient in ["A001", "G001", "PUTH_FAP5", "SCORT_A02", "SCORT_A03", "SCORT_A06", "SCORT_A07", "SCORT_A08", "SCORT_C06", "SCORT_C07", "SCORT_C08"]:
        return True
    else:
        return False

def check_same_CN(seg1, seg2):
    return seg1["tcn_em"] == seg2["tcn_em"] and seg1["lcn_em"] == seg2["lcn_em"]

def extend_merge_onesample(segs, telomeres, is_male, is_doubled):
    # telomeres dict {chr_name: [start_of_chr, end_of_chr]}
    all_chrom = ["chr"+str(x) for x in range(1,23)]
    all_chrom += ["chrX"]
    to_return = []
    for chrom in all_chrom:
        chrom_segs = segs[segs["chrom"]==chrom].sort_values("loc_start")
        
        if chrom != "chrX" or not is_male:
            normal_tcn = 2
            normal_lcn = 1
        else:
            normal_tcn = 1
            normal_lcn = 0
            
        if is_doubled:
            normal_tcn *= 2
            normal_lcn *= 2
        
        if len(chrom_segs) == 0:
            to_add = segs.iloc[0]
            to_add["chrom"] = chrom
            to_add["cf_em"] = 1
            to_add["loc_start"] = telomeres[chrom][0]
            to_add["loc_end"] = telomeres[chrom][1]
            to_add["tcn_em"] = normal_tcn
            to_add["lcn_em"] = normal_lcn
            to_return.append(to_add)
            continue
        
        curr_loc = telomeres[chrom][0]
        next_seg = chrom_segs.iloc[0]        
        next_seg["loc_start"] = telomeres[chrom][0]
        for i in range(len(chrom_segs)):
            curr_seg = next_seg
            if i == len(chrom_segs)-1:
                curr_seg["loc_end"] = telomeres[chrom][1]
                to_return.append(curr_seg)
                break
            next_seg = chrom_segs.iloc[i+1]
            if check_same_CN(curr_seg, next_seg):
                next_seg["loc_start"] = curr_seg["loc_start"]
                next_seg["cf_em"] = max(next_seg["cf_em"], curr_seg["cf_em"])
            elif next_seg["loc_start"] == curr_seg["loc_end"]:
                to_return.append(curr_seg)
            else:
                curr_diff = abs(curr_seg["tcn_em"] - normal_tcn)
                next_diff = abs(next_seg["tcn_em"] - normal_tcn)
                curr_diff_lcn = abs(curr_seg["lcn_em"] - normal_lcn)
                next_diff_lcn = abs(next_seg["lcn_em"] - normal_lcn)
                
                if curr_diff > next_diff or (curr_diff == next_diff and curr_diff_lcn > next_diff_lcn):
                    new_breakpoint = curr_seg["loc_end"]
                else:
                    new_breakpoint = next_seg["loc_start"]
                curr_seg["loc_end"] = new_breakpoint
                to_return.append(curr_seg)
                next_seg["loc_start"] = new_breakpoint
    return pd.concat([pd.DataFrame(x).transpose() for x in to_return], ignore_index=True)

def extend_merge_all(segs, telomeres, doubled_dict):
    to_return = []
    for sample in list(set(segs["sample_id"])):
        to_merge = segs[segs["sample_id"] == sample]
        to_return.append(extend_merge_onesample(to_merge, telomeres, check_if_male(to_merge.iloc[0]["patient"]), doubled_dict[sample]))
    to_return = pd.concat(to_return, ignore_index=True)
    return to_return[["chrom", "loc_start", "loc_end", "tcn_em", "lcn_em", "cf_em", "sample_id", "patient"]]

def fraction_genome_altered(segs, genome_len, doubled_dict):
    all_samples = list(set(segs["sample_id"]))
    to_return = []
    
    autosomes = segs[~np.isin(segs["chrom"], ["chrX", "chrY"])]
    for sample in all_samples:
        only_sample = autosomes[autosomes["sample_id"] == sample]
        
        normal_tcn = 2
        normal_lcn = 1
        if doubled_dict[sample]:
            normal_tcn = 4
            normal_lcn = 2
        
        not_normal_cn = only_sample[~np.logical_and(np.logical_or(only_sample["tcn_em"] == normal_tcn, pd.isna(only_sample["tcn_em"])), np.logical_or(only_sample["lcn_em"] == normal_lcn, pd.isna(only_sample["lcn_em"])))]
        not_normal_cn["length"] = not_normal_cn["loc_end"] - not_normal_cn["loc_start"]
        to_return.append(np.sum(not_normal_cn["length"]/genome_len))
    return pd.DataFrame({"FGA":to_return}, index=all_samples)
