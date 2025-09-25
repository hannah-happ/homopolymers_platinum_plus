from cyvcf2 import VCF
import pysam
import argparse
from tqdm import tqdm
from collections import Counter
import pandas as pd
import numpy as np
import gzip
from typing import List, Tuple, Union

MATCH, INS, DEL = range(3)
OP2DIFF = {MATCH: 0, INS: 1, DEL: -1}
# from SAM spec (https://samtools.github.io/hts-specs/SAMv1.pdf)
# page 8 -- 1 if the operation consumes reference sequence
CONSUMES_REF = dict(zip(range(9), [1, 0, 1, 1, 0, 0, 0, 1, 1]))

def get_bam(sample, manifest):
    """
    Retrieves path to sample's bam file.
    
    Args:
        sample: Sample to find the bam file for.
    
    Returns:
        AlignmentFile object for the sample's bam.
    """
    row = manifest[manifest['sample_id'] == sample]
    bam_path = row.bam.values[0]
    bam = pysam.AlignmentFile(bam_path, "rc")
    
    return bam
    
def get_family_bams(manifest, sample_list, family_id):
    """
    Reads in all the alignment files for a family.
    
    Args:
        manifest: Pandas DataFrame with family_id, sample_id, and path to bam or cram file for each sample.
        sample_list: List of IDs for each family member.
        family_ID: integer of the unique ID for each family.
    
    Returns:
        List of (name of member, AlignmentFile for bam) for mom, dad, proband
    """   
    pro_id = [x for x in sample_list if 'p' in x][0]
    
    mom = ("mom", get_bam("%s.mo" % family_id, manifest))
    dad = ("dad", get_bam("%s.fa" % family_id, manifest))
    pro = (pro_id, get_bam(pro_id, manifest))
    
    return [mom, dad, pro]
    
    
def bp_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
    """
    Simple utility to determine amount of overlap between 2 interval coordinates
    
    Args:
        s1 (int): start of interval 1
        e1 (int): end of interval 1
        s2 (int): start of interval 2
        e2 (int): end of interval 2

    Returns:
        int: size of overlap
    """
    return max(
        max((e2 - s1), 0) - max((e2 - e1), 0) - max((s2 - s1), 0),
        0,
    )

def count_indel_in_read(
    ct: List[Tuple[int, int]],
    rs: int,
    vs: int,
    ve: int,
    slop: int = 1,
) -> int:
    """
    Count up inserted and deleted sequence in a read using the pysam cigartuples object.
    a cigartuples object is a list of tuples -- each tuple stores a CIGAR operation as 
    its first element and the number of bases attributed to that operation as its second element.
    exact match is (0, N), where N is the number of bases that match the reference, insertions are
    (1, N), and deletions are (2, N). We only count CIGAR operation that completely overlap the
    expected STR locus interval.
    
    Args:
        ct (Tuple[int, int]): pysam cigartuples object
        rs (int): start of read w/r/t ref
        vs (int): start of TR locus in reference
        ve (int): end of TR locus in reference
        
    Returns:
        Int: net ins/del sequence in read w/r/t reference
    """
    
    # keep a running total of net ins/del operations in the read.
    cigar_op_total = 0
    # count from the first position in the read w/r/t the reference
    cur_pos = rs
    
    # loop over the cigartuples object, operation by operation
    for op, bp in ct:
        # check if this operation type consumes reference sequence
        # if so, we'll keep track of the bp associated with the op
        # if not,
        if bool(CONSUMES_REF[op]):
            op_length = bp
        else:
            op_length = 0
        # if the operation isn't an INS, DEL, or MATCH, we can just increment the running
        # start and move on. e.g., if the operation is a mismatch or a soft clip, we're 
        # not interested in incrementing our net CIGAR totals by the number of bp affected
        # by that operation. however, we *do* need to increment our current position.
        if op not in (INS, DEL, MATCH):
            cur_pos += op_length
            continue
        else:
            # keep track of the *relative* start and end of the operation,
            # given the number of consumable base pairs that have been encountered in the
            # iteration so far.
            op_s, op_e = cur_pos, cur_pos + max([1, op_length])
            # figure out the amount of overlap between this operation and our TR locus.
            # increment our counter of net CIGAR operations by this overlap.
            overlapping_bp = bp_overlap(op_s, op_e, vs - slop, ve + slop)
            if overlapping_bp > 0:
                  cigar_op_total += (bp * OP2DIFF[op])
            # increment our current position counter regardless of whether the operation
            # overlaps our STR locus of interest
            cur_pos += op_length
            
    return cigar_op_total
    
def get_read_diff(
    read: pysam.AlignedSegment,
    start: int,
    end: int,
    min_mapq: int = 60,
    slop: int = 1,
) -> Union[None, int]:
    """
    Compare a single sequencing read to the reference. then, count up the net inserted/
    deleted sequence in the read.
    
    Args:
        read (pysam.AlignedSegment): pysam read (aligned segment) object.
        start (int): start position of the TR locus in the reference.
        end (int): end position of the TR locus.
        min_mapq (int, optional): minimum mapping quality for a read to be considered. Defaults to 60.
        slop (int, optional): amount of slop around the start and end of the variant reported site. Defaults to 1.
        
    Returns:
        Union[None, int]: either None (if the read fails basic checks) or the net ins/del in read w/r/t reference
        """
    # initial filter on mapping quality
    if read.mapping_quality < min_mapq:
        return None
        
    # get the start and end positions of the read w/r/t the reference
    qs, qe = read.reference_start, read.reference_end
    
    # ensure that this read completely overlaps the TR locus along
    # with slop. if the read starts after the adjusted/slopped start
    # or if it ends before the adjusted/slopped end, skip the read
    adj_start, adj_end = start - slop, end + slop
    if qs > adj_start or qe < adj_end:
        return None
        
    # query the CIGAR string in the read.
    diff = count_indel_in_read(
        read.cigartuples,
        qs,
        start,
        end,
        slop=slop,
    )
    
    return diff

def extract_diffs_from_bam(
    bam_info,
    chrom: str,
    start: int,
    end: int,
    min_mapq: int = 60,
) -> List[Tuple[int, int]]:
    """
    Gather information from all reads aligned to the specified region. Extract
    the net ins/del in each read w/r/t the reference. Count up the number of reads
    with each net ins/del value.
    
    Args:
        bam_info: Tuple of sample name and pysam.AlignmentFile object
        chrom (str): chromosome we want to query
        start (int): start of TR locus
        end (int): end of TR locus
        min_mapq (int, optional): min MAPQ required for reads. Defaults to 60.
        
    Returns:
        List[Tuple[int, int]]: List of (X, Y) tuples where X is the read "diff" and Y is the
        number of reads with that "diff"
    """
    sample_name, bam = bam_info
    diffs = []
    
    if bam is None:
        diffs.append(0)
    else:
        for read in bam.fetch(chrom, start, end):
            diff = get_read_diff(
                read,
                start,
                end,
                slop = int(0.1 * (end - start)),
                min_mapq=min_mapq,
            )
            
            if diff is None:
                continue
            else:
                diffs.append(diff)
                
    # count up all recorded diffs between reads and reference allele
    diff_counts = Counter(diffs).most_common()
    return diff_counts   
       
def validate(mom_diffs, dad_diffs, child_diffs, threshold_diffs, dnm_min_reads):
    """
    Evaluate diff_counts from mom, dad, and child to determine if the child has a denovo diff_count.
    Ensure the dnm is supported by a minimum number of reads.

    Args:
        mom_diffs: List of diff_counts tuples for mom [(diff_size, read_count), ...]
        dad_diffs: List of diff_counts tuples for dad [(diff_size, read_count), ...]
        child_diffs: List of diff_counts tuples for child [(diff_size, read_count), ...]
        threshold_diffs: Threshold number of reads with the de novo allele to say that a parent definitely has it.
    
    Returns:
        Tuple (validation_status, denovo_diff_size) where validation_status is either 'true_de_novo' or 'false_positive', and denovo_diff_size is the size of the de novo diff_count if applicable, otherwise None.
    """
    
    # Convert the parents' diff_counts into dictionaries for easier lookup
    mom_diff_dict = dict(mom_diffs)
    dad_diff_dict = dict(dad_diffs)
    
    # Iterate over the child's diff_counts
    for diff_size, child_count in child_diffs:
        # Check if this diff_size is present in mom or dad with read count > threshold_diffs
        mom_count = mom_diff_dict.get(diff_size, 0)
        dad_count = dad_diff_dict.get(diff_size, 0)
        
        if mom_count <= threshold_diffs and dad_count <= threshold_diffs and child_count >= dnm_min_reads:
            # If neither parent has this diff_size above the threshold, it's a potential de novo mutation
            return 'true_de_novo', diff_size, child_count
    
    # If no de novo mutation is found, return false_positive
    return 'false_positive', None, None

def infer_repeat_motif(denovo_allele, motif_length):
    """
    Identify the most common repeat motif in an STR allele sequence.
    Args:
        allele: STR allele sequence (string)
        motif_length: repeat unit size to consider
    Returns:
        Most common repeat motif
    """
    # Convert the denovo allele to a string
    denovo_allele_str = str(denovo_allele)
    
    motif_counts = Counter()
    
    # Count occurrences of motifs of the specified length
    for i in range(len(denovo_allele_str) - motif_length + 1):
        motif = denovo_allele_str[i:i + motif_length]
        motif_counts[motif] += 1

    # Get the most common motif
    if motif_counts:
        most_common_motif, _ = motif_counts.most_common(1)[0]
        return most_common_motif
    else:
        return None  # Return None if no motifs are found
   
def process_family_candidate_variants(family, family_variants, bams, ref, threshold, dnm_min_reads, motif_period):
    """
    Process all candidate variants and print results
    
    Args:
        family: family ID
        family_variants: DataFrame of candidate variants for the family
        bams: DataFrame of bam file manifest
        ref: pysam.FastaFile object for reference genome
        threshold: max allowable read evidence in a parent
        dnm_min_reads: min number of reads supporting de novo allele
        motif_period: repeat unit size to consider
    
    Returns:
        None (writes results to output file)
    """    
    
    validated_dnm_df_cols = ["chrom", "start", "end", "child", "family", "child_DP", "mom_DP", "dad_DP", "child_diffcounts", "mom_diffcounts", "dad_diffcounts", "denovo_allele", "motif", "validation_status", "denovo_diff_size", "dnm_read_support"]
    
    dnms = []
    samples = bams.loc[bams['sample_id'].str.contains(str(family))]['sample_id'].values
    family_bams = get_family_bams(bams, samples, family)
    
    # loop through families and process each variant
    for idx, row in family_variants.iterrows():
        denovo_child = bams.loc[bams['ceph_id'] == str(row['child'])]['sample_id'].values[0]
        
        # extract chrom, start, end, child, family
        chrom, ref_start = row['chrom'], row['pos']
        ref_length = row['ref_allele_len']
        start = int(ref_start) - 1
        ref_end = start + int(ref_length)
        end = int(ref_end)
        child = row['child']
        
        #extract denovo allele
        child_alleles = set(row['child_alleles'].split('|'))
        parental_alleles = set(row['mat_alleles'].split('|') + row['pat_alleles'].split('|'))
        denovo_allele = child_alleles - parental_alleles

        #extract motif
        motif = infer_repeat_motif(denovo_allele, motif_period)

        #extract depths
        child_dp = row['child_DP'].replace("[", "").replace("]", "").strip()
        mom_dp = row['mat_DP'].replace("[", "").replace("]", "").strip()
        dad_dp = row['pat_DP'].replace("[", "").replace("]", "").strip()
            
        # Initialize variables for each family member's diffcounts
        child_diffcounts = []
        mom_diffcounts = []
        dad_diffcounts = []
        
        for bam_info in family_bams:
            sample_name, bam = bam_info

            diffcounts = extract_diffs_from_bam(bam_info, chrom, start, end, min_mapq=60)
            
            if sample_name == "mom":
                mom_diffcounts = diffcounts
            elif sample_name == "dad":
                dad_diffcounts = diffcounts
            else:
                child_diffcounts = diffcounts
                    
        # Validate whether the child has a true denovo mutation
        validation_status, denovo_diff_size, dnm_read_support = validate(mom_diffcounts, dad_diffcounts, child_diffcounts, args.threshold, args.dnm_min_reads)
                  
        validated_dnm_rows = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "child": child,
            "family": family,
            "child_DP": child_dp,
            "mom_DP": mom_dp,
            "dad_DP": dad_dp,
            "child_diffcounts": str(child_diffcounts),
            "mom_diffcounts": str(mom_diffcounts),
            "dad_diffcounts": str(dad_diffcounts),
            "denovo_allele": list(denovo_allele),
            "motif": motif,
            "validation_status": validation_status,
            "denovo_diff_size": denovo_diff_size,
            "dnm_read_support": dnm_read_support
        }
        
        dnms.append(validated_dnm_rows)
        
    validated_dnms = pd.DataFrame(dnms, columns=validated_dnm_df_cols)
    # validated_dnms.to_csv(f"{child}_validated_dnms_t{threshold}.csv", sep='\t', index=False)
    validated_dnms.to_csv(args.output, sep='\t', index=False)

def main(bams, candidate_dnm_file, ref, threshold, dnm_min_reads, output, motif_period):
    bams = pd.read_csv(args.bams, sep='\t')
    variants = pd.read_csv(args.input, sep='\t')
    ref = pysam.FastaFile(args.reference)

    families = list(set(variants['family']))

    # loop through families and process each variant
    for family in tqdm(families, total=len(families)):
        family_variants = variants.loc[variants['family'] == family]
        process_family_candidate_variants(family, family_variants, bams, ref, args.threshold, args.dnm_min_reads, args.motif_period)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
    ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
    ap.add_argument("-t", "--threshold", required=False, help="max allowable read evidence in a parent",
                    default=1, type=int)
    ap.add_argument("-r", "--reference", required=False, help="Reference genome fasta")
    ap.add_argument("-f", "--flanking_length", required=False, help="Length of flanking sequence",
                    default=5, type=int)
    ap.add_argument("-d", "--dnm_min_reads", required=False, help="min number of reads supporting de novo allele",
                    default=5, type=int)
    ap.add_argument("-o", "--output", required=False, help="Output file name")
    ap.add_argument("-p", "--motif_period", required=True, help="Motif period length", type=int)
    args = ap.parse_args()
    
    #run main function
    main(args.bams, args.input, args.reference, args.threshold, args.dnm_min_reads, args.output, args.motif_period)