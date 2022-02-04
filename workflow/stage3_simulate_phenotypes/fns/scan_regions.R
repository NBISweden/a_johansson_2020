scan_regions <- function(regions, freqs, min_nmrk) {
  # Scan through all the regions and check which ones are valid, i.e. which ones have at least N markers that fulfill
  # specified criteria and thus are suitable for simulation
  num_regions <- dim(cds_uniq)[1]
  scan_data <- data.frame(index=1:num_regions, gene=rep("",times=num_regions), num_valid_mrk=rep(0, times=num_regions))
  valid_markers_lst <- vector("list", length = num_regions)
  for (region_idx in 1:num_regions) {
    cds  <- cds_uniq[region_idx,]
    valid_markers <- freqs %>%
      filter(pos >= cds$start & pos <= cds$stop)
    scan_data[region_idx,] <- c(region_idx, cds$gene, dim(valid_markers)[1])
    valid_markers_lst[[region_idx]] <- valid_markers
  }
  scan_data <- scan_data %>% as_tibble()
  valid_scan_data <- scan_data %>% filter(num_valid_mrk >= min_nmrk)
  return(valid_scan_data)
}
