#!/bin/bash -l
 
#SBATCH -A snic2021-22-574
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7:00:00
#SBATCH -J download_ensembl_regulatory_and_extract_promoters

ensembl_path="http://ftp.ensembl.org/pub/current_regulation/homo_sapiens/" 
ensembl_file="homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz"
workdir=`pwd`
cd /home/marcin/hadriens_project/data/
if [ ! -f ${ensembl_file} ]; then
	wget ${ensembl_path}${ensembl_file}
else
	echo "The ${ensembl_file} exists! Skipping download..."
fi
zcat ${ensembl_file} | \
grep "promoter" | \
grep -v "Flanking" | \
awk -F"[\t;:]" '{sub(/ /,"_",$12)} {sub(/description=/,"",$12)} {print "chr"$1, $4+1, $5+2, $12}' | \
sort -k1,1 -k2,2n > all_GCRh38_promoters.bed
cd ${workdir}
