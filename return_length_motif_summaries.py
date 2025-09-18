#!/usr/bin/env python
#usage: python return_candidate_denovos.py -i <input_vcf> -p <ped_file> -s <motif_size> -dp <min_depth>

from cyvcf2 import VCF
import pandas as pd
import os
from collections import Counter, defaultdict


def load_pedigree(ped_file):
    """
    load in the ped file and return a DataFrame with pedigree information
    Args:
        ped_file (str): Path to the pedigree file in .ped format
    Returns:
        pd.DataFrame: DataFrame containing pedigree information
    """
    ped_cols = ['famID', 'indID', 'patID', 'matID', 'sex', 'pheno']
    ped_df = pd.read_csv(ped_file, sep=r'\s+', names=ped_cols)
    
    return ped_df
    
def get_trios(ped_df):
    """
    Identify trios within the ped file. return trios.
    Args:
        ped_df (pd.DataFrame): DataFrame containing pedigree information
    Returns:
        list of trios in format: (child_id, mother_id, father_id)
    """
    trios = []
    for idx, row in ped_df.iterrows():
        if row['patID'] != '0' and row['matID'] != '0':
            trios.append((row['indID'], row['matID'], row['patID']))
    return trios
    
def get_sample_index(mom, dad, kid, samples):
    """
    Returns index positions for mom, dad, and kid in the cyvcf2 samples list
    Args:
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2
    Returns:
        tuple: index positions for mom, dad, and kid in the format: (mom_col, dad_col, kid_col)
    """
    #get index positions for mom, dad, kid
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)
    
    return mom_col, dad_col, kid_col
    
def get_variant_type(variant, motif_size):
    """
    Returns True if the variant is desired motif size,
    which HipSTR annotates under PERIOD (1 = homopolymer, 2 = dinucleotide, etc).
    Args:
        variant: cyvcf2 Variant object (single VCF record)
        motif_size (int): Desired motif size (e.g., 1 for homopolymers, 2 for dinucleotide repeats, etc.)
    Returns:
        str: "motif_match" if the variant matches the desired motif size, otherwise "no_motif_match"
    """
    if variant.INFO.get('PERIOD') == motif_size:
        return "motif_match"
    else:
        return "no_motif_match"
        
def is_sex_chromosome(variant):
    """
    Check whether variant is on a sex chromosome
    Args:
        variant: cyvcf2 Variant object (single VCF record)
    Returns:
        bool: True if variant is on chrX or chrY, otherwise False
    """
    chrom = variant.CHROM
    if chrom == "chrX" or chrom == "chrY":
        return True
    else:
        return False

def get_missing_genotype_info(variant, mom, dad, kid, samples):
    """
    Identifies whether any individual in the trio has a missing gt, specifies which individual(s) have missing gts
    Args:
        variant: cyvcf2 Variant object (single VCF record)
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2
    Returns:
        dict: counts for missing maternal, paternal, or child gts and total number of missing genotypes per site
    """

    #get index positions for mom, dad, kid
    mom_col, dad_col, kid_col = get_sample_index(mom, dad, kid, samples)
    gts = variant.gt_types
    
    #initialize a dict to store counts of missing gts
    missing_gt_counts = {
        "mom_missing": 0,
        "dad_missing": 0,
        "mom_and_dad_missing": 0,
        "kid_missing": 0,
        "total_missing": 0,
    }
    
    #check if any gt is missing
    if gts[mom_col] == 2 and gts[dad_col] == 2:
        missing_gt_counts["mom_and_dad_missing"] += 1
        missing_gt_counts["total_missing"] += 2
    else:
        if gts[mom_col] == 2:
            missing_gt_counts["mom_missing"] += 1
            missing_gt_counts["total_missing"] += 1
        if gts[dad_col] == 2:
            missing_gt_counts["dad_missing"] += 1
            missing_gt_counts["total_missing"] += 1
    
    if gts[kid_col] == 2:
        missing_gt_counts["kid_missing"] += 1
        missing_gt_counts["total_missing"] += 1
        
    return missing_gt_counts

def get_coverage_info(variant, mom, dad, kid, samples, min_dp):
    """
    Access DP (depth of coverage) for mom, dad, kid and determine if all individuals meet minimum depth requirement
    Args:
        variant: cyvcf2 Variant object (single VCF record)
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2 indexing
        min_dp (int): Minimum depth of coverage required for each individual in the trio (main argument)
    Returns:
        str: "sufficient_covg" if all individuals meet minimum depth requirement, otherwise "insufficient_covg"
        tuple: (mom_DP, dad_DP, kid_DP) if sufficient coverage, otherwise not returned
    """
    
    #get index positions for mom, dad, kid
    mom_col, dad_col, kid_col = get_sample_index(mom, dad, kid, samples)
    
    #get DP and PDP for mom, dad, kid
    mom_DP = variant.format('DP')[mom_col]
    dad_DP = variant.format('DP')[dad_col]
    kid_DP = variant.format('DP')[kid_col]
    
    if mom_DP < min_dp or dad_DP < min_dp or kid_DP < min_dp:
        return "insufficient_covg"
    else:
        return "sufficient_covg", mom_DP, dad_DP, kid_DP
        
def is_queried_for_dnm(variant):
    """
    Check if variant meets criteria to be queried for DNM status (sufficient coverage, no missing genotypes, correct motif size)
    Args:
        variant: cyvcf2 Variant object (single VCF record)
    Returns:
        bool: True if variant has already been queried for DNM status, otherwise False
    
    """

    return True

def analyze_trio(input_vcf, kid, mom, dad, motif_size, min_dp, bin_size, min_length):
    """
    Analyze a single trio and return child_id, total_queried, bin_counts (dict), max_ref_len, motif_counts (dict)
    Args:
        input_vcf (str): Path to the HipSTR VCF file (main argument)
        kid (str): Child's sample ID
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        motif_size (int): Motif size (e.g., 1 for homopolymers, 2 for dinucleotide repeats, etc.) (main argument)
        min_dp (int): Minimum depth of coverage required for each individual in the trio (main argument)
        bin_size (int): Size of bins for reference allele lengths (e.g., 15 bp) (main argument)
        min_length (int): Minimum reference allele length to consider (e.g., 10 bp) (main argument)
        output_dir (str): Directory to save output files (main argument)
    Returns:
        tuple: (child_id, total_queried, bin_counts (dict), max_ref_len, motif_counts (dict))        
    """
    #load VCF
    vcf = VCF(input_vcf)

    #initialize counters and storage
    queried = 0
    bin_counts = defaultdict(int)
    motif_counts = Counter()
    max_ref_len = 0

    #function to determine bin label
    def get_bin(length):
        if length < min_length:
            return None
        start = min_length + bin_size * ((length - min_length) // bin_size)
        end = start + bin_size - 1
        return f"{start}-{end}"

    #iterate through VCF and analyze each variant
    for variant in vcf:
        if get_variant_type(variant, motif_size) != "motif_match":
            continue
        if is_sex_chromosome(variant):
            continue
        if any(get_missing_genotype_info(variant, mom, dad, kid, vcf.samples).values()):
            continue
        if get_coverage_info(variant, mom, dad, kid, vcf.samples, min_dp) == "insufficient_covg":
            continue
        
        #site is queried
        queried += 1

        #length bin tallying
        ref_len = len(variant.REF)
        max_ref_len = max(max_ref_len, ref_len)
        bin_label = get_bin(ref_len)
        if bin_label:
            bin_counts[bin_label] += 1

        #motif tallying
        period = variant.INFO.get('PERIOD')
        ref = variant.REF
        motif_substrings = [ref[i:i+period] for i in range(len(ref) - period + 1)]
        if motif_substrings:
            most_common_motif, _ = Counter(motif_substrings).most_common(1)[0]
            motif_counts[most_common_motif] += 1

    return kid, queried, bin_counts, max_ref_len, motif_counts

def main(input_vcf, ped_file, motif_size, min_dp, output_dir, bin_size, min_length):
    """
    Main function to analyze all trios in the pedigree and summarize reference allele length distributions and motifs.
    Saves length distribution summary and motif summary to output directory.
    Args:
        input_vcf (str): Path to the HipSTR VCF file (main argument)
        ped_file (str): Path to the pedigree file in .ped format (main argument)
        motif_size (int): Motif size (e.g., 1 for homopolymers, 2 for dinucleotide repeats, etc.) (main argument)
        min_dp (int): Minimum depth of coverage required for each individual in the trio (main argument)
        output_dir (str): Directory to save output files (main argument)
        bin_size (int): Size of bins for reference allele lengths (e.g., 15 bp) (main argument)
        min_length (int): Minimum reference allele length to consider (e.g., 10 bp) (main argument)
    Returns:
        None
    """
    pedigree = load_pedigree(ped_file)
    trios = get_trios(pedigree)

    #for length summary
    all_length_rows = []
    all_bins = set()
    max_seen_len = 0

    #for motif summary
    all_motif_rows = []
    all_motifs = set()

    for kid, mom, dad in trios:
        child, queried, bin_counts, max_ref_len, motif_counts = analyze_trio(
            input_vcf, kid, mom, dad, motif_size, min_dp, bin_size, min_length
        )
        
        #track bins and max length
        max_seen_len = max(max_seen_len, max_ref_len)
        all_bins.update(bin_counts.keys())

        #build length row
        length_row = {"child": child, "total_queried": queried}
        length_row.update(bin_counts)
        all_length_rows.append(length_row)

        #track motif set and build motif row
        all_motifs.update(motif_counts.keys())
        motif_row = {"child": child, "total_queried": queried}
        for motif, count in motif_counts.items():
            motif_row[f"motif_{motif}"] = count
        all_motif_rows.append(motif_row)

    # consistent columns for length bins
    bin_labels = []
    if max_seen_len >= min_length:
        for start in range(min_length, max_seen_len + 1, bin_size):
            end = start + bin_size - 1
            bin_labels.append(f"{start}-{end}")

    length_df = pd.DataFrame(all_length_rows)
    for bin in bin_labels:
        if bin not in length_df.columns:
            length_df[bin] = 0
    length_df = length_df[["child", "total_queried"] + bin_labels]
    length_df.to_csv(
        f"{output_dir}/query_length_bin_summary_{motif_size}.csv",
        sep="\t", index=False
    )

    #consistent columns for motifs
    motif_df = pd.DataFrame(all_motif_rows)
    motif_cols = [f"motif_{motif}" for motif in sorted(all_motifs)]
    for c in motif_cols:
        if c not in motif_df.columns:
            motif_df[c] = 0
    motif_df = motif_df[["child", "total_queried"] + motif_cols]
    motif_df.to_csv(
        f"{output_dir}/query_motif_summary_{motif_size}.csv",
        sep="\t", index=False
    )

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_vcf", required=True, help="HipSTR VCF file")
    ap.add_argument("-p", "--ped_file", required=True, help="Pedigree file in .ped format")
    ap.add_argument("-s", "--motif_size", required=True, type=int, help="Motif size (e.g., 1 for homopolymers, 2 for dinucleotide repeats, etc.)")
    ap.add_argument("-dp", "--minimum_depth", required=True, type=int)
    ap.add_argument("-o", "--output_dir", required=True, help="Output directory")
    ap.add_argument("-b", "--bin_size", type=int, default=15, help="Bin size for reference allele length (bp)")
    ap.add_argument("-m", "--min_length", type=int, default=10, help="Minimum reference allele length to consider (bp)")

    args = ap.parse_args()
    
    #run main function
    main(args.input_vcf, args.ped_file, args.motif_size, args.minimum_depth, args.output_dir, args.bin_size, args.min_length)
    


