#!/usr/bin/env python
"""
Return candidate de novo mutations (DNMs) from a HipSTR VCF for all trios in a pedigree.

This script scans a HipSTR VCF, applies basic queryability filters (motif size,
autosomes only, no missing genotypes, per-sample DP >= threshold), and flags
sites where the child's allele(s) are absent in both parents as candidate DNMs.
It writes a per-child TSV of candidate DNMs and an exclusions/queried summary TSV.

Notes
-----
- Motif size is determined from HipSTR's INFO/PERIOD (1 = homopolymer, 2 = dinucleotide, ...).
- Genotypes are taken from `variant.gt_types` and alleles from `variant.gt_bases`.
- This is a liberal candidate DNM definition intended for downstream filtering/validation.
"""

from typing import Dict, Tuple, List
from cyvcf2 import VCF
import pandas as pd
import os

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

def get_allele_lens(variant, mom, dad, kid, samples):
    """
    Returns allele lengths for mom, dad, and kid
    Args:
        variant: cyvcf2 Variant object (single VCF record)
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2
    Returns:
        tuple[str, str, str]: allele lengths for mom, dad, and kid in the format: (mom_
    """
    
    #get index positions for mom, dad, kid
    mom_col, dad_col, kid_col = get_sample_index(mom, dad, kid, samples)
    
    #access gt bases for mom, get length of 2 alleles
    mom_gt_bases = variant.gt_bases[mom_col].split('|')
    mom_allele1_len = len(mom_gt_bases[0])
    mom_allele2_len = len(mom_gt_bases[1])
    mom_allele_lens = str(mom_allele1_len) + "," + str(mom_allele2_len)
    
    #access gt bases for dad, get length of 2 alleles
    dad_gt_bases = variant.gt_bases[dad_col].split('|')
    dad_allele1_len = len(dad_gt_bases[0])
    dad_allele2_len = len(dad_gt_bases[1])
    dad_allele_lens = str(dad_allele1_len) + "," + str(dad_allele2_len)
    
    #access gt bases for kid, get length of 2 alleles
    kid_gt_bases = variant.gt_bases[kid_col].split('|')
    kid_allele1_len = len(kid_gt_bases[0])
    kid_allele2_len = len(kid_gt_bases[1])
    kid_allele_lens = str(kid_allele1_len) + "," + str(kid_allele2_len)
    
    return mom_allele_lens, dad_allele_lens, kid_allele_lens

def get_alleles(variant, mom, dad, kid, samples):
    """
    Returns raw allele strings for mom, dad, and kid
    Args:
        variant: cyvcf2 Variant object (single VCF record)
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2
    Returns:
        tuple[str, str, str]: raw allele strings for mom, dad, and kid in the format: (mom_alleles, dad_alleles, kid_alleles) (e.g., "AAAA|AAAAAA")
    """
    #get index positions for mom, dad, kid
    mom_col, dad_col, kid_col = get_sample_index(mom, dad, kid, samples)

    mom_alleles = variant.gt_bases[mom_col]
    dad_alleles = variant.gt_bases[dad_col]
    kid_alleles = variant.gt_bases[kid_col]    
    
    return mom_alleles, dad_alleles, kid_alleles

def get_candidate_dnms(variant, mom, dad, kid, samples):
    """
    Flag candidate DNMs by checking whether each child allele is absent in both parents.

    Args:
        variant: cyvcf2 Variant object.
        mom (str): Mother's sample ID
        dad (str): Father's sample ID
        kid (str): Child's sample ID
        samples (list): needed for cyvcf2 indexing

    Returns:
        str: "candidate_dnm" if any child allele is absent in both parents; else "inherited".

    Notes
    -----
    - This is deliberately liberal and may capture Mendelian errors and artifacts.
      Downstream filters (read-level evidence, stutter, strand balance, etc.) are expected.
    """
    
    #get index positions for mom, dad, kid
    mom_col, dad_col, kid_col = get_sample_index(mom, dad, kid, samples)
    
    #access gt bases for mom, dad, kid
    mom_gt_bases = variant.gt_bases[mom_col].split('|')
    dad_gt_bases = variant.gt_bases[dad_col].split('|')
    kid_gt_bases = variant.gt_bases[kid_col].split('|')
    
    #make a list to append if criteria are met
    kid_dnms = []
    for allele in kid_gt_bases:
        #if allele differs from mom and differs from dad, append to dnm list
        #if one allele differs from mom/dad, dnm list will have a length of 1
        #if both alleles differ from mom and dad, dnm list will have a length of 2
        if allele not in mom_gt_bases and allele not in dad_gt_bases:
            kid_dnms.append(allele)
    
    if len(kid_dnms) > 0:
        return "candidate_dnm"
    else:
        return "inherited"
        
def is_perfect_repeat(variant):
    """
    Determine if the reference allele is a perfect (uninterrupted) repeat.
    Args:
        variant: cyvcf2 Variant object.
    Returns:
        bool: True if REF consists of the same base repeated; otherwise False.
    """
    unique_characters = set(variant.REF)
    if len(unique_characters) == 1:
        return True
    else:
        return False

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

def analyze_trio(input_vcf, kid, mom, dad, motif_size, min_dp, output_dir):
    """
    
    """
    vcf = VCF(input_vcf)
    queried = 0

    exclusion_counts = {
        "sex_chromosome": 0,
        "mom_missing": 0,
        "dad_missing": 0,
        "kid_missing": 0,
        "one_gt_unknown": 0,
        "two_gt_unknown": 0,
        "three_gt_unknown": 0,
        "insufficient_covg": 0,
        "other_exclusions": 0,
    }
    
    #initialize an empty DF to store results
    dnm_df_cols = [
        "chrom",
        "pos",
        "family",
        "child",
        "child_gt",
        "child_DP",
        "mat_gt",
        "mat_DP",
        "pat_gt",
        "pat_DP",
        "len_repeat_unit",
        "ref_allele_len",
        "is_perfect",
        "child_alleles",
        "mat_alleles",
        "pat_alleles"
    ]
    
    dnms = []
        
    for variant in vcf:
        variant_type = get_variant_type(variant, motif_size)
        chromosome = is_sex_chromosome(variant)
        missing_gts = get_missing_genotype_info(variant, mom, dad, kid, vcf.samples)

        #impose exclusion rules
        #first, skip motif sizes you aren't interested in (will exclude the majority)
        if variant_type != "motif_match":
            exclusion_counts["other_exclusions"] += 1
            continue
        #skip if sex chromosome
        if chromosome == True:
            exclusion_counts["sex_chromosome"] += 1
            continue
        #skip if any of the genotypes are missing
        if missing_gts["total_missing"] > 0:
            #track individual missing gts
            exclusion_counts["mom_missing"] += missing_gts["mom_missing"]
            exclusion_counts["dad_missing"] += missing_gts["dad_missing"]
            exclusion_counts["kid_missing"] += missing_gts["kid_missing"]

            #track total missing counts per site
            if missing_gts["total_missing"] == 1:
                exclusion_counts["one_gt_unknown"] += 1
            elif missing_gts["total_missing"] == 2:
                exclusion_counts["two_gt_unknown"] += 1
            elif missing_gts["total_missing"] == 3:
                exclusion_counts["three_gt_unknown"] += 1
            continue
        #skip if number of reads used for hipSTR genotyping in mom, dad, or kid < min_dp
        covg = get_coverage_info(variant, mom, dad, kid, vcf.samples, min_dp)
        if covg == "insufficient_covg":
            exclusion_counts["insufficient_covg"] += 1
            continue

        #if you reach this point, the variant has been queried for dnm status
        if is_queried_for_dnm(variant) == True:
            queried += 1
        
        
        #if one of proband's alleles is not present in mom or dad, this is a candidate dnm
        #this is very liberal and candidate dnms will need additional filtering
        if get_candidate_dnms(variant, mom, dad, kid, vcf.samples) == "candidate_dnm":
            mat_gt, pat_gt, child_gt = get_allele_lens(variant, mom, dad, kid, vcf.samples)
            chrom = variant.CHROM
            pos = variant.POS
            family = kid.split("-")[0]
            child = kid
            len_repeat_unit = variant.INFO.get('PERIOD')
            ref_allele_len = len(variant.REF)
            is_perfect = is_perfect_repeat(variant)
            mat_alleles, pat_alleles, child_alleles = get_alleles(variant, mom, dad, kid, vcf.samples)
            _, mat_DP, pat_DP, child_DP = get_coverage_info(variant, mom, dad, kid, vcf.samples, min_dp)
            
            #create a dictionary for row data
            dnm_df_row = {
                "chrom": chrom,
                "pos": pos,
                "family": family,
                "child": child,
                "child_gt": child_gt,
                "child_DP": child_DP,
                "mat_gt": mat_gt,
                "mat_DP": mat_DP,
                "pat_gt": pat_gt,
                "pat_DP": pat_DP,
                "len_repeat_unit": len_repeat_unit,
                "ref_allele_len": ref_allele_len,
                "is_perfect": is_perfect,
                "child_alleles": child_alleles,
                "mat_alleles": mat_alleles,
                "pat_alleles": pat_alleles,
            }
            
            dnms.append(dnm_df_row)
            
    candidate_dnms = pd.DataFrame(dnms, columns=dnm_df_cols)
    candidate_dnms.to_csv(f"{output_dir}/{kid}_candidate_dnms_{motif_size}bp.csv", sep='\t', index=False)
        
    #return queried count for this trio
    return {
        "child": kid,
        "queried_count": queried,
        "sex_chromosome": exclusion_counts["sex_chromosome"],
        "mom_missing": exclusion_counts["mom_missing"],
        "dad_missing": exclusion_counts["dad_missing"],
        "kid_missing": exclusion_counts["kid_missing"],
        "one_gt_unknown": exclusion_counts["one_gt_unknown"],
        "two_gt_unknown": exclusion_counts["two_gt_unknown"],
        "three_gt_unknown": exclusion_counts["three_gt_unknown"],
        "insufficient_covg": exclusion_counts["insufficient_covg"],
        "other_exclusions": exclusion_counts["other_exclusions"],    
    }
        
def main(input_vcf, ped_file, motif_size, min_dp, output_dir):
    """
    Analyze all trios and save query counts to a summary CSV file
    """
    #load pedigree file
    pedigree = load_pedigree(ped_file)
    
    #extract trios
    trios = get_trios(pedigree)
    
    #collect queried counts for all trios
    queried_summary = []
    
    for kid, mom, dad in trios:
        trio_summary = analyze_trio(input_vcf, kid, mom, dad, motif_size, min_dp, output_dir)
        queried_summary.append(trio_summary)
        
    queried_summary_df = pd.DataFrame(queried_summary)
    queried_summary_df.to_csv(f"{output_dir}/query_summary_{motif_size}.csv", sep='\t', index=False)

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_vcf", required=True, help="HipSTR VCF file")
    ap.add_argument("-p", "--ped_file", required=True, help="Pedigree file in .ped format")
    ap.add_argument("-s", "--motif_size", required=True, type=int, help="Motif size (e.g., 1 for homopolymers, 2 for dinucleotide repeats, etc.)")
    ap.add_argument("-dp", "--minimum_depth", required=True, type=int)
    ap.add_argument("-o", "--output_dir", required=True, help="Output directory")
    
    args = ap.parse_args()
    
    #run main function
    main(args.input_vcf, args.ped_file, args.motif_size, args.minimum_depth, args.output_dir)