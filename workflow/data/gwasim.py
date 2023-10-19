import sys
from pybgen import PyBGEN
import numpy as np
import pandas as pd
import random
import argparse
import time

def set_settings_from_command_line_args(S):
    """
    Set settings in the SettingsSingleton (S) based on command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Set simulation parameters.")
    parser.add_argument("--bed_regions_path", type=str, default="cds_test.bed", help="Path to the BED file.")
    parser.add_argument("--bgen_file_path", type=str, default="Rum_recoded_repos_norel_rnd3000_chr22.bgen", help="Path to the BGEN file.")
    parser.add_argument("--threshold", type=float, default=0.05, help="Threshold for MAF rare vs. common.")
    parser.add_argument("--my_chr", type=str, default="1", help="Chromosome of interest.")
    parser.add_argument("--num_common", type=int, default=1, help="Number of common variants for simulation.")
    parser.add_argument("--num_rare", type=int, default=2, help="Number of rare variants for simulation.")
    parser.add_argument("--rare_eff_mean", type=float, default=2, help="Mean effect size for rare variants.")
    parser.add_argument("--rare_eff_std_dev", type=float, default=0.3, help="Std. deviation of effect size for rare variants.")
    parser.add_argument("--common_eff_mean", type=float, default=0.3, help="Mean effect size for common variants.")
    parser.add_argument("--common_eff_std_dev", type=float, default=0.01, help="Std. deviation of effect size for common variants.")
    parser.add_argument("--err_mean", type=float, default=0.05, help="Mean for the error term in phenotype simulation.")
    parser.add_argument("--err_sd", type=float, default=0.01, help="Std. deviation for the error term in phenotype simulation.")
    parser.add_argument("--num_sim", type=int, default=1, help="Number of simulations.")
    
    args = parser.parse_args()

    # Access the singleton instance and set values
    S.bed_regions_path = args.bed_regions_path
    S.bgen_file_path = args.bgen_file_path
    S.threshold = args.threshold
    S.my_chr = args.my_chr
    S.num_common = args.num_common
    S.num_rare = args.num_rare
    S.rare_eff_mean = args.rare_eff_mean
    S.rare_eff_std_dev = args.rare_eff_std_dev
    S.common_eff_mean = args.common_eff_mean
    S.common_eff_std_dev = args.common_eff_std_dev
    S.err_mean = args.err_mean
    S.err_sd = args.err_sd
    S.num_sim = args.num_sim

class SettingsSingleton:
    """
    Singleton class to store user-defined settings.

    Attributes:
        bed_regions_path (str): Path to the BED file containing region data.
        bgen_file_path (str): Path to the BGEN file.
        threshold (float): Alleles with a frequency below this threshold are considered rare.
        my_chr (str): Chromosome of interest.
        num_common (int): Number of common variant markers used for simulation.
        num_rare (int): Number of rare variant markers used for simulation.
        rare_eff_mean (float): Mean effect size for rare variants.
        rare_eff_std_dev (float): Standard deviation of effect size for rare variants.
        common_eff_mean (float): Mean effect size for common variants.
        common_eff_std_dev (float): Standard deviation of effect size for common variants.
        err_mean (float): Mean for the error term in phenotype simulation.
        err_sd (float): Standard deviation for the error term in phenotype simulation.
        num_sim(int): Number of simulations to perform.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(SettingsSingleton, cls).__new__(cls)
            cls._instance.init_settings()
        return cls._instance

    def init_settings(self):
        # Define the path to your regions bed file and BGEN file
        #self.bed_regions_path = "cds_chr_fixed_pruned.bed"
        self.bed_regions_path = "cds_test.bed"
        self.bgen_file_path = "Rum_recoded_repos_norel_rnd3000_chr22.bgen"
        self.threshold = 0.05
        self.my_chr = "1"
        self.num_common = 1
        self.num_rare = 2
        self.rare_eff_mean = 2
        self.rare_eff_std_dev = 0.3
        self.common_eff_mean = 0.3
        self.common_eff_std_dev = 0.01
        self.err_mean = 0.05
        self.err_sd = 0.01
        self.num_sim = 1

        
def get_maf(variant):
    """
    Calculate Minor Allele Frequency (MAF) from variant data.

    Args:
        variant (object): Variant data.

    Returns:
        float: MAF value.
    """
    gt = list(variant[1])
    num_alleles = 2 * len(gt)
    num_aa = gt.count(0.0)
    num_aA = gt.count(1.0)
    num_AA = gt.count(2.0)
    f_a1 = (2 * num_aa + num_aA) / num_alleles
    f_a2 = (2 * num_AA + num_aA) / num_alleles
    f_minor = min(f_a1, f_a2)
    f_major = max(f_a1, f_a2)
    #if (f_minor != threshold):
    #    print("f_minor: ", f_minor, ", f_major: ", f_major)
    return(f_minor)

def read_bed_file(S):
    """
    Read and parse region data from a BED file.

    Args:
        S (SettingsSingleton): Singleton object containing user-defined settings.

    Returns:
        list: A list of dictionaries containing region information.
    """
    filter_chr = S.my_chr
    file_path = S.bed_regions_path
    data = []  # Initialize a list to store the data
  
    with open(file_path, 'r') as bed_file:
        for line in bed_file:
            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Ensure the line has the correct number of fields (4)
            if len(fields) == 4:
                chromosome, start, stop, region_name = fields
                if (filter_chr == "all" or chromosome == filter_chr): 
	                # Convert start and stop to integers if needed
	                start = int(start)
	                stop = int(stop)
            	    # Append the data as a tuple or dictionary, depending on your preference
                	# Here, it's appended as a dictionary
                	data.append({
                    	'chromosome': chromosome,
                    	'start': start,
                    	'stop': stop,
                    	'region_name': region_name,
                    	'num_rare': 0,
                    	'num_common': 0,
                    	'is_valid': 'no'
	                })
    #print(data)
    return data

def validate_regions(regions, S):
    """
    Validate regions based on simulation criteria.

    Args:
        regions (list): List of region dictionaries.
        S (SettingsSingleton): Singleton object containing user-defined settings.

    Returns:
        list: A filtered list of validated regions.
    """
    bgen_file_path = S.bgen_file_path
    my_chr = S.my_chr
    threshold = S.threshold
    num_common = S.num_common
    num_rare = S.num_rare
    
    with PyBGEN(bgen_file_path) as bgen:
        #print("There are", bgen.nb_variants, "variants and", bgen.nb_samples, "samples in the file.")
        while True:
            try:
                variant = bgen.next()
                chr = variant[0].chrom

                # HACK
                if variant[0].pos == 0:
                    variant[0].pos = int(variant[0].name.split('_')[3])

                position = variant[0].pos

                if chr == my_chr:
                    for region in regions:
                        if region['start'] <= position <= region['stop']:
                            maf = get_maf(variant)
                            if 0 < maf <= threshold:
                                region['num_rare'] += 1
                            elif maf > threshold:
                                region['num_common'] += 1

                        if region['num_common'] >= num_common and region['num_rare'] >= num_rare:
                            region['is_valid'] = 'yes'
            except StopIteration:
                break
            except Exception as e:
                sys.stderr.write(f"ERROR. Chromosome {chr}, position {position}: {e}. Skipping over...\n")
    valid_regions = [item for item in regions if item['is_valid'] != 'no']
    return(valid_regions)

def fix_gt(gt):
    """
    Fix genotype data.

    Args:
        gt (numpy.ndarray): Genotype data.

    Returns:
        numpy.ndarray: Fixed genotype data.
    """
    # Step 1: Convert the array to integers
    gt = np.array(gt)
    gt = gt.astype(int)

    # Step 2: Check and swap 0 and 2 if necessary, i.e. if minor and major allele are mixed up
    count_0 = np.count_nonzero(gt == 0)
    count_1 = np.count_nonzero(gt == 1)
    count_2 = np.count_nonzero(gt == 2)

    if count_2 > count_0:
        gt[gt == 0] = 2
        gt[gt == 2] = 0

    # Step 3: Impute, i.e. replace any other values (missing data) with the most frequent value among 0, 1, and 2
    most_frequent_value = np.argmax([count_0, count_1, count_2])

    for i in range(len(gt)):
        if gt[i] not in (0, 1, 2):
            gt[i] = most_frequent_value
    return(gt)

def select_variants(bgen, region, S):
    """
    Select variants for simulation.

    Args:
        bgen (PyBGEN): PyBGEN object for reading BGEN data.
        region (dict): Region dictionary containing region information.
        S (SettingsSingleton): Singleton object containing user-defined settings.

    Returns:
        dict: A dictionary containing selected variants and metadata for simulation.
    """
    filter_chr = S.my_chr
    threshold = S.threshold
    selected_variants = []
    num_rare = S.num_rare
    num_common = S.num_common
    start = region['start']
    stop = region['stop']
    rare = []
    common = []
    for variant in bgen.iter_variants_in_region(filter_chr, start, stop):
        #print(variant)
        maf = get_maf(variant)
        if (maf <= threshold):
            rare.append(variant[0].name)
        else:
            common.append(variant[0].name)
            
    # Randomly select variants
    sel_rare = random.sample(rare, num_rare)
    sel_common = random.sample(common, num_common)
    
    result = {}
    result['region'] = region
    result['selected_variants'] = sel_rare + sel_common
    return(result)

def simulate_phenotype(bgen, variants, S):
    """
    Simulate phenotypes based on selected variants.

    Args:
        bgen (PyBGEN): PyBGEN object for reading BGEN data.
        variants (dict): Dictionary containing selected variants for simulation.
        S (SettingsSingleton): Singleton object containing user-defined settings.

    Returns:
        dict: A dictionary containing simulated phenotype data and metadata.
    """

    region = variants['region']
    selected = variants['selected_variants']
    variants_table = pd.DataFrame(columns=['variant_name', 'chr', 'pos', 'type', 'eff', 'maf'])
    y = np.zeros(bgen.nb_samples)    # Create an empty array for phenotypes
    for variant_name in selected:
        variant = bgen.get_variant(variant_name)[0]
        name = variant[0].name
        chr = variant[0].chrom
       
        # HACK
        if variant[0].pos == 0:
            variant[0].pos = int(variant[0].name.split('_')[3])
        pos = variant[0].pos

        maf = get_maf(variant)
        if (maf <= S.threshold):
            type = 'rare'
            eff = np.random.normal(S.rare_eff_mean, S.rare_eff_std_dev)
        else:
            type = 'common'
            eff = np.random.normal(S.common_eff_mean, S.common_eff_std_dev)
           
        tmp = pd.DataFrame([{'variant_name': name, 'chr': chr, 'pos': pos, 'type': type, 'eff': eff, 'maf': maf}])
        variants_table = pd.concat([variants_table, tmp], ignore_index=True)
        gt = fix_gt(variant[1])
        y_contrib = eff * gt
        #print(y_contrib[y_contrib > 0]) # sanity check
        y = y + y_contrib
    y = y + np.random.normal(S.err_mean, S.err_sd)
    result = {'region': region, 'variants_data': variants_table, 'y': y}     
    return(result)
     
if __name__ == '__main__':
    # User-defined variables are stored in a singleton
    S = SettingsSingleton()
    set_settings_from_command_line_args(S)
    
    # Read file with regions
    print("Reading regions file...", end = ' ')
    start_time = time.time()
    regions = read_bed_file(S)
    exec_time = round(time.time() - start_time, 2)
    print(f" done in {exec_time}s.")
    n_reg = len(regions)
    print(f"Validating {n_reg} regions on chromosome {S.my_chr} from {S.bed_regions_path}...", end = ' ')
    start_time = time.time()
    valid_regions = validate_regions(regions, S)
    exec_time = round(time.time() - start_time, 2)
    print(f" done in {exec_time}s.")
    print(f"Found {len(valid_regions)} region(s) matching simulation criteria.")
    print(f"Drawing {S.num_sim} region(s) for simulation (w. replacement)", end = ' ')
    start_time = time.time()
    sim_regions = random.choices(validate_regions(regions, S), k=S.num_sim)
    exec_time = round(time.time() - start_time, 2)
    print(f" done in {exec_time}s.")

    cnt = 1
    for region in sim_regions:
        start_time = time.time()
        print(f"Performing simulation {cnt} of {S.num_sim}...")
        print(f"\t - processing region {region['region_name']}")
        print(f"\t\t - scanning variants in the region")
        with PyBGEN(S.bgen_file_path) as bgen:
            print(f"\t\t - selecting variants")
            variants = select_variants(bgen, region, S)
            #print(variants)
            print(f"\t\t - simulating phenotype")
            sim_result = simulate_phenotype(bgen, variants, S)
            #print(sim_result)
            y_mean = round(np.mean(sim_result['y']), 2)
            y_sd = round(np.std(sim_result['y']), 2)
            print(f"\t\t - simulated phenotype mean = {y_mean}, std. dev = {y_sd}")
            exec_time = round(time.time() - start_time, 2)
            print(f"\t\t - done in {exec_time}s.")
        cnt += 1
