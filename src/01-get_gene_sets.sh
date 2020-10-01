#!/bin/bash
 
# Burden test pipeline -- part 1.
# Retrieval and preparation of promoter and exonic regions info.
# marcin.kierczak@scilifelab.se
# last changes: 21.09.2020
 
# User defined variables
# Currently for hg19 (GRCh37)
FTP_GENCODE="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh37_mapping/gencode.v35lift37.annotation.gff3.gz"
FTP_ENSEMBL="ftp://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20191101.gff.gz"

# Project variables
PROJDIR=`pwd`
RESOURCES_DIR=${PROJDIR}"/resources"
 
# Download and prepare the required resources
mkdir ${RESOURCES_DIR} 2>/dev/null
 
#------------------------------------------ Exons ------------------------------------------
 
# Exon coordinates from GENCODE
GENCODE_FILE=$(basename "${FTP_GENCODE}" .gz)
cd ${RESOURCES_DIR}
if [ ! -f ${GENCODE_FILE} ]; then
    echo "Downloading GENCODE annotations"
    wget ${FTP_GENCODE}
    gunzip *.gz
fi
 
# Extract exonic regions only, filter out pseudogenes and antisense RNAs and save them in a bed file
echo "Processing ${GENCODE_FILE}..."
if [ ! -f "exons.bed" ]; then
    awk '/^chr/ {print}' ${GENCODE_FILE}| \
    awk -F"[\t;:]" '{sub(/gene_type=/,"",$15); sub(/gene_name=/,"",$16); sub(/gene_id=/,"",$13)} $3=="exon" {print $1,$4-1,$5,$7,$3,$15,$16,$13}' | \
    grep -v ".*pseudogene.*" | \
    grep -v "antisense" \
    > exons.bed
fi
 
#------------------------------------------ CDS ------------------------------------------
 
 
if [ ! -f "CDS.bed" ]; then
        awk '/^chr/ {print}' ${GENCODE_FILE}| \
        sed 's/chr//' | \
        awk '$3=="CDS" {OFS=" ";print $1"_"$4"_"$5,$1,$4-1,$5,$7,"CDS",$9}' | \
        grep -v ".*pseudogene.*"  | \
        sed "s/ID=.*gene_id=//" | \
        sed "s/;.*gene_type=/    /" | \
        sed "s/;.*gene_name=/    /" | \
        sed 's/;.*//'  | \
        sort -k1,1 -u | \
        awk '{print $2,$3,$4,$5,$6,$8,$7,$9}' | \
        sort -k1,1 -k2n,2 \
        > CDS.bed
fi
#------------------------------------------ Promoters ------------------------------------------
 
# Promoter coordinates from Ensembl
ENSEMBL_FILE=$(basename "${FTP_ENSEMBL}" .gz)
if [ ! -f ${ENSEMBL_FILE} ]; then
    echo "Downloading Ensembl annotations"
    wget ${FTP_ENSEMBL}
    gunzip *.gz
fi
 
# Extract promoter regions only and make a bed file out of them
# Shift start by -1 due to bed being 0-based
if [ ! -f "promoters.bed" ]; then
    grep "promoter" ${ENSEMBL_FILE} | \
    grep -v "Flanking" | \
    awk -F"[\t;:]" '{sub(/ /,"_",$12)} {sub(/description=/,"",$12)} {print "chr"$1, $4-1, $5, $12}' | \
    sort -k1,1 -k2,2n > promoters.bed
fi
 
#------------------------------------------ Match promoters with exons ------------------------------------------
 
# Define upstream and downstream (of the exon) regions. Usually 500 or 1000bp.
#UPSTREAM=2000
#DOWNSTREAM=0
 
# Matching exons with promoters -- create extended promoter regions
# hrM 15956 16023 - exon Mt_tRNA gene_status=KNOWN ENSG00000210196.2
#cat exons.bed | \
#awk '{ if ($4 == "+") { print $1,$2-'$UPSTREAM', $2+'$DOWNSTREAM', $4, $5, $6,$7,$8 } else if ($4 == "-") { print $1, $3-'$UPSTRERAM', $3+'$DOWNSTREAM', $4,$5,$6,$7,$8 }}' \
#> promoter_regions.bed
