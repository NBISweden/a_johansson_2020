import sys
from pybgen import PyBGEN
import numpy as np
import pandas as pd
import random
from multiprocessing import Pool

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

def read_bed_file(file_path, filter_chr = "all"):
    data = []  # Initialize a list to store the data

    with open(file_path, 'r') as bed_file:
        for line in bed_file:
            # Split the line into fields using tab as the delimiter
            fields = line.strip().split('\t')

            # Ensure the line has the correct number of fields (4)
            if len(fields) == 4:
                chromosome, start, stop, region_name = fields
                if filter_chr == "all" or chromosome == filter_chr: 
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
    return data

def validate_regions(bgen_file_path, my_chr, regions, threshold, num_common, num_rare):
    with PyBGEN(bgen_file_path) as bgen:
        #print("There are", bgen.nb_variants, "variants and", bgen.nb_samples, "samples in the file.")
        while True:
            try:
                variant = bgen.next()
                chr = variant[0].chrom

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

if __name__ == '__main__':

    # User-defined variables
    # Define the path to your regions bed file and BGEN file
    #bed_regions_path = "cds_chr_fixed_pruned.bed"
    bed_regions_path = "cds_test.bed"
    bgen_file_path = "Rum_recoded_repos_norel_rnd3000_chr22.bgen"
    threshold = 0.05 	# Alleles with frequency below this threshold are considered rare
    my_chr = "1" 	    # chromosome of interest
    num_common = 1      # number of common variant markers used for simulation
    num_rare = 2        # number of rare variant markers used for simulation 
    rare_eff_mean = 2 
    rare_eff_std_dev = 0.3

    # Read file with regions
    print("Reading regions...")
    regions = read_bed_file(bed_regions_path, filter_chr = 'chr' + my_chr)
    n_reg = len(regions)
    print(f"Validating {n_reg} regions...")
    valid_regions = validate_regions(bgen_file_path, my_chr, regions, threshold, num_common, num_rare)
    for region in valid_regions:
    	print(region)
    start = valid_regions[0]['start']
    end = valid_regions[0]['stop']
    common = []
    rare = []
    with PyBGEN(bgen_file_path) as bgen:
        for variant in bgen.iter_variants_in_region(my_chr, start, end):
            maf = get_maf(variant)
            if (maf <= threshold):
                rare.append(variant[0].name)
            else:
                common.append(variant[0].name)
    result = {}
    result['region'] = region
    selected_rare = random.sample(rare, num_rare)
    selected_common = random.sample(common, num_common)
    selected_variants = pd.DataFrame(columns=['variant_name', 'chr', 'pos', 'type', 'eff', 'maf'])
    with PyBGEN(bgen_file_path) as bgen:
        for variant_name in selected_rare:
            variant = bgen.get_variant(variant_name)[0]
            name = variant[0].name
            chr = variant[0].chrom
            if variant[0].pos == 0:
                variant[0].pos = int(variant[0].name.split('_')[3])
            pos = variant[0].pos
            type = 'rare'
            eff = np.random.normal(rare_eff_mean, rare_eff_std_dev)
            maf = get_maf(variant)
            tmp = pd.DataFrame([{'variant_name': name, 'chr': chr, 'pos': pos, 'type': type, 'eff': eff, 'maf': maf}])
            selected_variants = pd.concat([selected_variants, tmp], ignore_index=True)
        print(selected_variants)    
