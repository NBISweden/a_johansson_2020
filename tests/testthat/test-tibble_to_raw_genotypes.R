test_that("tibble to raw genotype conversion works", {
     test_data <- tibble(
       chr = c(1),
       pos = c(1001,1002,1005,1017),
       snp = c('snp1', 'snp2', 'snp3', 'snp4'),
       strand = c('u'),
       ind1 = c(0,1,1,2),
       id2 = c(2,NA,1,1),
       john = c(0,0,1,0),
       jill = c(1,1,1,2),
       elon = c(2,0,NA,1)
     )
  local_edition(3)
  helper <- function() {
    path <- tempfile(fileext = ".raw")
    tibble_to_raw_genotypes(x = test_data, output = path, progress = 0)
    return(path)
  }
  expect_snapshot_file(helper(), "genotypes.raw")
})

test_that("tibble to raw genotype conversion stops when some weird genotype coding is provided", {
  test_data <- tibble(
    chr = c(1),
    pos = c(1001,1002,1005,1017),
    snp = c('snp1', 'snp2', 'snp3', 'snp4'),
    strand = c('u'),
    ind1 = c(0,1,1,2),
    id2 = c(2,NA,7,1),
    john = c(0,0,1,0),
    jill = c(1,1,1,2),
    elon = c(2,0,NA,1)
  )

  expect_error(tibble_to_raw_genotypes(x = test_data, output = path, progress = 1))
})
