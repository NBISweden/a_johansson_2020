library("optparse")

option_list = list(
  make_option(c("--chr-name"), type="character", default="chr22",
              help="name of the chromosome top process, e.g. chr22"),
  make_option(c("--data-path"), type="character", default="/home/marcin/ExomeSeq2ndRel/marcin/",
              help="path to data files"),
  make_option(c("--regions-bed"), type="character", default="GRCh38_chr22_cds.bed",
              help="bed file defining functional regions of interest, eg. CDS"),
  make_option(c("--afreq"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.afreq",
              help="an afreq file with allele frequencies computed by Plink2"),
  make_option(c("--fam"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.fam",
              help="fam file"),
  make_option(c("--bim"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.bim",
              help="bim file"),
  make_option(c("--plink-prefix"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22",
              help="prefix (no extension) of the plink2 data files"),
  make_option(c("--min-maf"), type="numeric", default = 0.005,
              help = "minimal maf to consider (used to remove ultra-rare/fixed markers)"),
  make_option(c("--rare-maf"), type="numeric", default = 0.01,
              help = "markers above this threshold will be considered common"),
  make_option(c("--num-mrk"), type="numeric", default = 5,
              help = "number of markers per region to be used for simulating phenotype"),
  make_option(c("--num-mrk-neg"), type="numeric", default = 0,
              help = "per region number of markers with negative effect"),
  make_option(c("--mean-eff"), type="numeric", default = 1,
              help = "mean effect of a rare allele"),
  make_option(c("--sd-eff"), type="numeric", default = 0.01,
              help = "std. dev. for the distribution of rare allele effects"),
  make_option(c("--mean-err"), type="numeric", default = 0,
              help = "mean error (residuals)"),
  make_option(c("--sd-err"), type="numeric", default = 0.05,
              help = "std. dev. for the distribution errors (residuals)"),
  make_option(c("--num-sim"), type="numeric", default = 10,
              help = "number of simulations to run (regions to use)")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
