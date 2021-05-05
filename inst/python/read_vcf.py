import io
import os
import pandas as pd
import vcf
import pyranges as pr


bed_path = "/Users/kiero/Projects/johansson/promoters_tab_delim.bed"
vcf_path = "/Users/kiero/Projects/johansson/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# read bed file
# df = pr.read_bed(bed_path, as_df=True)
# for chr, start, end in zip(df['Chromosome'], df['Start'], df['End']):
#     region_name = chr + "_" + str(start) + "_" + str(end)
#     try:
#         vcf_data = vcf.Reader(filename = vcf_path).fetch(chr.replace('chr', ''), start, end)
#         print(region_name)
#         print(vcf_data)
#     except:
#         continue
#         #  print('Region not in the supplied vcf file!')
# print("Done!")

vcf_path = "/Users/kiero/Projects/johansson/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
region = '22:51220999-51223201'
def read_region(f, region):
    reader = vcf.Reader(filename = vcf_path).fetch(region)
    df = pd.DataFrame([vars(r) for r in reader])
    out = df.merge(pd.DataFrame(df.INFO.tolist()),
                   left_index=True, right_index=True)
    return out

#data = read_region(vcf_path, region)
#print(data)


import vcf
import pandas as pd

def vcf_region_reader(path, region):
    print("Reading ", region, " from: ", path)
    try:
        reader = vcf.Reader(filename = path).fetch(region)
        var_names = reader.samples
        index_names = []
        df_list = []
        for record in reader:
            line = []
            index_names.append(str(record.CHROM) + "_" + str(record.POS))
            for sample in record.samples:
                enc = ''
                gtype = sample['GT']
                if(gtype == '0|0' or gtype == '0/0'):
                    enc = '0'
                elif(gtype == '1|1' or gtype == '1/1'):
                    enc = '2'
                else:
                    enc = '1'
                line.append(enc)
            df_list.append(line)
        df = pd.DataFrame(df_list, columns=var_names, index = index_names).T
        return(df)
    except Exception as e:
        print(e)
        print('skipping region...')

path = "/Users/kiero/Projects/johansson/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
region = '22:51220999-51223201'
data = vcf_region_reader(path, region)
print(data)
