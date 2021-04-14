# # global reference to scipy (will be initialized in .onLoad)
# pandas <- NULL
# pyvcf <- NULL
#
# .onLoad <- function(libname, pkgname) {
#   # use superassignment to update global reference to scipy
#   pandas <<- reticulate::import("pandas", delay_load = TRUE)
#   pyvcf <<- reticulate::import("pyvcf", delay_load = TRUE)
# }
