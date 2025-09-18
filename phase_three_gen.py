from cyvcf2 import VCF
import pandas as pd
from tqdm import tqdm
from typing import List, Dict
from collections import Counter
import argparse
import ast
import numpy as np
    
def get_str_inheritance(
    str_vcf: VCF,
    chrom: str,
    start: int,
    focal: str,
    kids: List[str],
    denovo_allele,
    smp2idx: Dict[str, int]
):
    """
    Determine which kids inherited the STR allele matching the focal individual's STR genotype?
    Args:
        str_vcf (VCF): HipSTR VCF file
        chrom (str): chromosome of the str site
        start (int): start position of the STR site
        focal (str): Sample ID of the focal individual
        kids (List[str]): list of children sample IDs
        denovo_allele: de novo STR allele genotype (From validated_dnms file)
        smp2idx (Dict[str, int]): Dictionary mapping sample IDs to VCF sample indexes.

    Returns:
       List[str]: list of children who inherited the STR allele.
    """
    
    inherited_kids = []
    
    try:
        denovo_allele_l = list(ast.literal_eval(denovo_allele))  # Safely parse the string to a list
    except (ValueError, SyntaxError) as e:
        raise ValueError(f"Failed to parse denovo_allele: {denovo_allele}") from e
    
    for variant in str_vcf(f"{chrom}:{start}-{start}"):
       # get focal individual's genotype
       focal_idx = smp2idx[focal]
       focal_alleles = variant.gt_bases[focal_idx].split('|')
       
       #debugging
#       print(f"Focal alleles: {focal_alleles}")
#       print(f"De novo allele: {denovo_allele_l}")
#        print(focal_idx)
#        focal_gt = variant.gt_types[focal_idx]
#        print(focal_gt)
#        print(focal_alleles)
#        print(denovo_allele_l)
       
       # confirm the de novo allele is one of the focal individual's alleles
       if not any(allele in focal_alleles for allele in denovo_allele_l):
           raise ValueError(f"De novo allele {denovo_allele_l} not found in focal individual's genotype: {focal_alleles}")
       
#        if denovo_allele not in focal_alleles:
#            raise ValueError(f"de novo allele {denovo_allele} not found in focal individual's genotype.")
       
       #check each kid for inheritance of the de novo allele	
       for kid in kids:
           kid_idx = smp2idx[kid] # get VCF index for the child
           kid_alleles = variant.gt_bases[kid_idx].split('|')
           #compare kid's alleles to de novo allele
           # if denovo_allele_l in kid_alleles:
#                inherited_kids.append(kid)
           if any(allele in kid_alleles for allele in denovo_allele_l):  # Check for any overlap
                inherited_kids.append(kid)
    
    return inherited_kids
  
  
def catalog_informative_sites(
    vcf: VCF,
    region: str,
    mom: str,
    dad: str,
    focal: str,
    focal_spouse: str,
    kids: List[str],
    kids_with_str: List[str],
    smp2idx: Dict[str, int],
    min_gq: int = 20,
    min_dp: int = 10,
):
    """
    Identify informative SNPs for phasing in the specified region
    """
    all_smps = [mom, dad, focal, focal_spouse] + kids
    all_idxs = np.array([smp2idx[s] for s in all_smps])
    print(all_idxs)

    mom_idx, dad_idx = smp2idx[mom], smp2idx[dad]
    focal_idx, spouse_idx = smp2idx[focal], smp2idx[focal_spouse]

    informative_sites = []
    for v in vcf(region):
        # Filter SNPs
        if v.var_type != "snp" or v.call_rate < 1.:
            continue
        if any(v.gt_quals[idx] < min_gq for idx in all_idxs):
            continue
        ref_depths, alt_depths = v.gt_ref_depths, v.gt_alt_depths
        total_depths = ref_depths + alt_depths
        if any(total_depths[idx] < min_dp for idx in all_idxs):
            continue
    
        # If parental genotypes don't match AND focal is het AND focal_spouse is HOM_REF,
        # we have an informative site
        
        # make sure parental genotypes don't match
        dad_gt, mom_gt = v.gt_types[dad_idx], v.gt_types[mom_idx]
        focal_gt, spouse_gt = v.gt_types[focal_idx], v.gt_types[spouse_idx]
        if mom_gt == dad_gt or focal_gt != 1 or spouse_gt != 0:
            continue
            
        #identify informative genotype (from kids with STR DNM)
        kid_gts = [v.gt_types[smp2idx[kid]] for kid in kids]
        informative_gt = None
        for kid, gt in zip(kids, kid_gts):
            #check if kid inherited STR DNM
            if kid in kids_with_str:
                if informative_gt is None:
                    informative_gt = gt # set informative genotype
#                elif informative_gt != gt:
 #                   informative_gt = None # perfect segregation is violated if a kid with DNM has a different genotype
  #                  break
       
        if informative_gt is None:
            continue #no consistent informative genotype found
        
        informative_parent = None
        if informative_gt == mom_gt:
            informative_parent = "mom"
        elif informative_gt == dad_gt:
            informative_parent = "dad"
        if informative_parent is None:
            continue
           

#         informative_parent = None
#         # figure out informative parent. as long as parental genotypes
#         # don't match AND the kid is HET, we have an informative site.
#         # we know for a fact that if the kid is HET, the parent with more ALTs
#         # donated the ALT allele (e.g., if kid is 0/1, dad is 1/1, and mom is 0/1,
#         # dad donated the 1 and mom donated the 0).
#         if dad_gt > mom_gt:
#             informative_parent = "dad"
#         if mom_gt > dad_gt:
#             informative_parent = "mom"
#         if informative_parent is None:
#             continue


        # loop over kids to catalog inheritance
        inf_string = []
        for kid, gt in zip(kids, kid_gts):
#             kid_idx = smp2idx[kid]
#             kid_gt = v.gt_types[kid_idx]
           
            # if kid has str, has_str = "Y", otherwise has_str = "N"
            has_str = "Y" if kid in kids_with_str else "N"
            if gt == informative_gt:
                inf_string.append(f"{kid}-INF-{has_str}")
            else:
                inf_string.append(f"{kid}-NON-INF-{has_str}")
        
        if len(inf_string) == 0:
            continue
        
        informative_sites.append(f"{v.CHROM}:{v.POS}:{informative_parent}:{'|'.join(inf_string)}")

    return informative_sites   


def check_for_dnm_inheritance(inf: str):
    """
    Check if the informative site is relevant to the inheritance of the DNM allele
    """
#     children = inf.split(":")[-1]  # Extract the children part of the string (e.g., "kid1-REF-Y|kid2-ALT-N")
#     inherited_with_dnm = [kid.split("-")[-1] == "Y" for kid in children.split("|")]  # Check inheritance of the DNM allele
#     return any(inherited_with_dnm)  # Return True if any child inherited the DNM allele

    children = inf.split(":")[-1] # Extract the children part of the string (e.g., "kid1-Y|kid2-N|kid3-Y")
    has_dnm = [kid.split("-")[-1] for kid in children.split("|")] # Extract the 'Y' or 'N' for each child
    # only interested in the inheritance patterns for which at least one kid inherited
    # the DNM allele
    return all([h == "N" for h in has_dnm]) # Return True if all children have 'N' (did not inherit the DNM)


def main(args):
    # Load validated DNMs
#    print("Loading validated DNMS...")
    validated_dnms = pd.read_csv(args.validated_dnms, sep='\t')
    validated_dnms = validated_dnms[validated_dnms['validation_status'] == 'true_de_novo']
    print(validated_dnms)
    
    # load SNP VCF
#    print("loading SNP VCF...")
    SNV_VCF = VCF(args.snv_vcf, gts012=True)
    SMP2IDX_SNV = dict(zip(SNV_VCF.samples, range(len(SNV_VCF.samples))))
    print(SMP2IDX_SNV)
#     SMP2IDX = {sample: idx for idx, sample in enumerate(SNV_VCF.samples)}
    
    #Load STR VCF
    STR_VCF = VCF(args.str_vcf)
    SMP2IDX_STR = dict(zip(STR_VCF.samples, range(len(STR_VCF.samples))))
    
    # Load pedigree
#    print("Loading pedigree...")
    pedigree = pd.read_csv(args.pedigree, sep='\t')
#    print(pedigree)
    
    results = []
    
#    print("Processing DNMS...")
    for _, row in tqdm(validated_dnms.iterrows(), total=validated_dnms.shape[0]):
        chrom = row['chrom']
        start = int(row['start']) + 1
        end = int(row['end'])
        focal = str(row['child'])
        print(focal)
        focal_type = type(focal)
        print(focal_type)
        denovo_allele = row['denovo_allele']
#         denovo_allele = row['denovo_allele'].replace("{", "").replace("}", "").strip()
        

#        test = pedigree['sample_id'] == focal
#        print(test)
        #Extract family relationships
        family = pedigree[pedigree['sample_id'] == focal]
        print(family)
        if family.empty:
            print(f"Skipping {focal}: Not found in pedigree.")
            continue
        	
        mom, dad = family['maternal_id'].values[0], family['paternal_id'].values[0]
        
        #Identify spouse of focal individual
#        spouse = None
#        spouse_row = pedigree[(pedigree['maternal_id'] == focal) | (pedigree['paternal_id'] == focal)]
#        if not spouse_row.empty:
#            spouse = spouse_row['sample_id'].values[0]
#        print(spouse)
#        if spouse is None or spouse not in SMP2IDX_SNV:
#            print(f"Skipping {focal}: Spouse not found or missing in VCF")
#            continue
                
        # identify children of the focal individual and spouse
        kids = pedigree[(pedigree['maternal_id'] == focal) | (pedigree['paternal_id'] == focal)]['sample_id'].tolist()
        
        ex_kid = pedigree[pedigree['sample_id'] == kids[0]]
        spouse = None
        spouse = ex_kid['paternal_id'].values[0] if ex_kid['paternal_id'].values[0] != focal else ex_kid['maternal_id'].values[0]
        print(spouse)
        if spouse is None or spouse not in SMP2IDX_SNV:
            print(f"Skipping {focal}: Spouse not found or missing in VCF")
            continue   

        # check kids for inheritance
        kids_with_str = get_str_inheritance(STR_VCF, chrom, start, focal, kids, denovo_allele, SMP2IDX_STR)
        
        # Define region with slop
        slop = 250_000
        region = f"{chrom}:{max(1, start - slop)}-{end + slop}"

        # Catalog informative SNPs
        informative_sites = catalog_informative_sites(
            vcf=SNV_VCF,
            region=region,
            mom=mom,
            dad=dad,
            focal=focal,
            focal_spouse=spouse,
            kids=kids,
            kids_with_str=kids_with_str, 
            smp2idx=SMP2IDX_SNV,
        )

        # Phase de novo mutation based on informative SNPs
        phase_result = {
            "chrom": chrom,
            "pos": start,
            "end": end,
            "family": row['family'],
            "child": focal,
            "denovo_allele": denovo_allele,
            "kids_with_str": ",".join(kids_with_str),
            "informative_sites": informative_sites,
        }

        if len(informative_sites) > 0:
            relevant_phase_infos = [p for p in informative_sites if not check_for_dnm_inheritance(p)]
#             relevant_phase_infos = [p for p in informative_sites if check_for_dnm_inheritance(p)]
            if len(relevant_phase_infos) == 0:
                phase_result.update({"most_common_hap": "UNK", "most_common_freq": 0.0, "candidate_postzygotic": False})
            else:
#                 inheritance_combos = Counter([":".join(p.split(":")[2:]) for p in relevant_phase_infos]).most_common()
                inheritance_combos = Counter([":".join(p.split(":")[2:]) for p in relevant_phase_infos]).most_common()  # Keep the `REF` or `ALT` and the STR DNM info
                most_common_hap = inheritance_combos[0][0] #returns most common haplotype
                most_common_freq = inheritance_combos[0][1] / len(relevant_phase_infos) #count of most common haplotype divided by total number of relevant informative sites
                # most common hap is formatted like
                # dad:NA12886-Y|NA12887-Y|NA12882-N|NA12879-Y|NA12883-Y
                # keep track of whether all of the children that inherited the informative
                # allele at this site also inherited the DNM
                is_pz = len(set(most_common_hap.split(":")[-1].split("|"))) == 2
                
                
                
#                 most_common_hap = inheritance_combos[0][0]
#                 most_common_freq = inheritance_combos[0][1] / total_combos
#                 most_common_children = [c.split("-")[-1] for c in most_common_hap.split(":")[-1].split("|")]
#                 is_pz = len(set(most_common_children)) == 2
                    
                phase_result.update({
                    "most_common_hap": most_common_hap,
                    "most_common_freq": most_common_freq,
                    "candidate_postzygotic": is_pz,
#                     "ref_alt_inheritance": most_common_hap
                })
        
        else:
            phase_result.update({"most_common_hap": "UNK", "most_common_freq": 0.0, "candidate_postzygotic": False})
             
        
        results.append(phase_result)
        
    
    #save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep='\t', index=False)       
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Phase validated de novo mutations using SNP data")
    parser.add_argument("-dnms", "--validated_dnms", required=True, help="validated DNMs file")
    parser.add_argument("-snvs", "--snv_vcf", required=True, help="SNV VCF file")
    parser.add_argument("-strs", "--str_vcf", required=True, help="HipSTR vcf")
    parser.add_argument("-ped", "--pedigree", required=True, help="pedigree file")
    parser.add_argument("-o", "--output", required=True, help="Output file for phased results")
    args = parser.parse_args()

    main(args)



