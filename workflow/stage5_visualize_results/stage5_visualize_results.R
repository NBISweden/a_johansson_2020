library(optparse)

option_list = list(
  make_option(c("--title"), type="character", default="Phemulator", help="title of the report"),
  make_option(c("--path"), type="character", default="~/Dropbox/WABI/Projects/Johansson_PP/output_bianca_analyses/2022-01-27_snakemake/",
              help="path to data"),
  make_option(c("--meta"), type="character", default="simulated_metadata.txt",
              help="path to metadata file"),
  make_option(c("--sim"), type="character", default="simulated_phenotypes.txt",
              help="path to file with simulated phenotypes")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

title = opt$title
path <- opt$path
sim_meta <- opt$meta
sim_phenos <- opt$sim

render_report = function(title = "Phemulator", path, meta, sim) {
  rmarkdown::render(
    "stage5_report_template.Rmd", params = list(
      title = title,
      path = path,
      meta = meta,
      sim = sim
    ),
    output_file = paste0("Report-", title, ".pdf")
  )
}

render_report(title = title, path = path, meta = sim_meta, sim = sim_phenos)
