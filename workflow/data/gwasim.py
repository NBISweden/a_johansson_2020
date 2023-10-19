import sys
from pybgen import PyBGEN
import numpy as np
import pandas as pd
import random
from multiprocessing import Pool

class SettingsSingleton:
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
        
def get_maf(variant):
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
    
    # Read file with regions
    print("Reading regions...")
    regions = read_bed_file(S)
    n_reg = len(regions)
    print(f"Validating {n_reg} regions on chromosome {S.my_chr} found in {S.bed_regions_path}...")
    valid_regions = validate_regions(regions, S)
    print(f"Found {len(valid_regions)} region(s) matching simulation criteria.")
    for region in valid_regions:
        print(f"\t - processing region {region['region_name']}")
        print(f"\t\t - scanning variants in the region")
        with PyBGEN(S.bgen_file_path) as bgen:
            print(f"\t\t - selecting variants")
            variants = select_variants(bgen, region, S)
            #print(variants)
            print(f"\t\t - simulating phenotype")
            sim_result = simulate_phenotype(bgen, variants, S)
            print(sim_result)

