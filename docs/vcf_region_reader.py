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
