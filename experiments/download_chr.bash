#!/bin/bash -l
 
#SBATCH -A snic2021-22-574
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 7:00:00
#SBATCH -J download_chr22
 

workdir=`pwd`
cd /home/marcin/hadriens_project/data/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20200515_EBI_Freebayescalls/ALL.chr22.freebayes.20200518.snps_indels.NYhc.GRCh38.vcf.gz
cd ${workdir}
