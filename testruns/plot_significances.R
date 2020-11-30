library('tidyverse')
promoters <- read_table2(file = 'testruns/promoters.bed',
             col_names = c('chr', 'start', 'end', 'descr'),
             col_types = 'ciic')

load('testruns/results/common_runs_23_11_2020/common_run_results_1_marker.Rda')
result1 <- -log10(result)
rm(result)
load('testruns/results/common_runs_23_11_2020/common_run_results_2_markers.Rda')
result2 <- -log10(result)
rm(result)

plot_significances <- function(simulated, lowest) {
  x_lim <- c(0, 1.05 * max(simulated))
  y_lim <- c(0, 1.05 * max(lowest))
  plot(x=simulated, y=lowest, xlim = x_lim, ylim = y_lim, type='n',
      cex.axis = .7, las = 1, xlab='-log10(p_sim)', ylab = '-log10(p_best)',
      main = "Significance of simulated vs. non-simulated.")
  min_sim <- min(simulated)
  min_best <- min(lowest)
  abline(a=0, b=1)
  abline(v=min_sim, h=min_best, col = c('red', 'blue'))
  mtext(text = paste0('noise = ', round(min_best, 2)), side = 2, line = 1, at = min_best, las=1, cex=.5, col='red')
  mtext(text = round(min_sim, 2), side = 1, line = 2, at = min_sim, las=1, cex=.5, col='blue')
  bad <- (lowest > simulated)
  good <- !bad
  points(x = simulated[good], y = lowest[good], pch = 19, cex = .8, col = "olivedrab")
  points(x = simulated[bad], y = lowest[bad], pch = 1, cex = .8, col = "darkgrey")
}

par(mfrow = c(2,1))
plot_significances(simulated = result1[,1], lowest = result1[,2])
plot_significances(simulated = pmax(result2[,1], result2[,2]), lowest = result2[,3])
par(mfrow = c(1,1))

plot(x=result2[,1], y=result2[,2], xlim = c(0, 60), ylim = c(0, 60), col = 'grey',
     type = 'n', cex = .4, las = 1, xlab='-log10(p_sim1)', ylab = '-log10(p_sim2)', pty='m', plot=F)
abline(a=0, b=1, col = 'lightgrey')
points(x=result2[,1], y=result2[,2], pch=1, col='darkgrey', cex=.4)


par(mfrow = c(3,1))
plot(x=result1[,1], y=result1[,2], xlim = c(0, 60), ylim = c(0, 60),
     pch = 19, cex = .4, las = 1, xlab='-log10(p_sim)', ylab = '-log10(p_best)')
plot(x=result2[,1], y=result2[,3], xlim = c(0, 60), ylim = c(0, 60), col = 'slateblue',
     pch = 19, cex = .4, las = 1, xlab='-log10(p_sim1)', ylab = '-log10(p_best)')
plot(x=result2[,2], y=result2[,3], xlim = c(0, 60), ylim = c(0, 60), col = 'tomato',
     pch = 19, cex = .4, las = 1, xlab='-log10(p_sim2)', ylab = '-log10(p_best)')
par(mfrow = c(1,1))

plot_phenos <- function(shape = c(.1, .1),
                          thr = 0.01,
                          factor = c(0.1, 0),
                          maf = seq(0, 0.5, 0.005),
                          mean = 0,
                          N = 1000) {
  ## Plot effects
  old_par <- par()
  par(cex.axis = .8, bty = 'n', las = 1)
  eff <- dbeta(maf, shape1 = shape[1], shape2 = shape[2])
  sd <- factor[1] * eff[abs(eff) < Inf] + factor[2]
  y_lim <- c(0 - 4 * max(sd), 2 * max(eff[eff < Inf]) + 4 * max(sd))
  ## Generate data
  eff_AA <- matrix(NA, nrow = N, ncol = length(eff))
  eff_aA <- matrix(NA, nrow = N, ncol = length(eff))
  eff_aa <- matrix(NA, nrow = N, ncol = length(eff))
  for (i in 1:N) {
    eff_AA[i,] <- eff * 0 + rnorm(n = length(eff), mean, sd)
    eff_aA[i,] <- eff * 1 + rnorm(n = length(eff), mean , sd)
    eff_aa[i,] <- eff * 2 + rnorm(n = length(eff), mean , sd)
  }
  eff_AA[,length(eff)] <- eff_AA[,length(eff)-1]
  eff_aA[,length(eff)] <- eff_aA[,length(eff)-1]
  eff_aa[,length(eff)] <- eff_aa[,length(eff)-1]
  eff_gt_lst <- list(eff_AA, eff_aA, eff_aa)
  plot(maf, eff, xlab='maf', ylab='phenotypic values', type = 'n', ylim = y_lim)
  grid(col = 'lightgrey', lty = 3)
  abline(v = thr, col = 'grey')
  mtext(text = paste('maf =', thr), 3, line = 0, at = thr, cex = .5, col = 'grey')
  points(maf, eff, cex = .5, col = 'darkgrey', pch = 19, type = 'b')
  gt_cols <- c('red', 'green', 'blue')
  gt_cols_alpha <- c(rgb(1,0,0,.3), rgb(0,1,0,.3), rgb(0,0,1,.3))
  for (i in 1:3) {
    eff_gt <- eff_gt_lst[[i]]
    gt_col <- gt_cols[i]
    gt_col_alpha <- gt_cols_alpha[i]
    mean_gt <- colMeans(eff_gt)
    sd_gt <- apply(eff_gt, 2, sd)
    upp_gt <- mean_gt + 3 * sd_gt
    low_gt <- mean_gt - 3 * sd_gt
    polygon(c(maf, rev(maf)), c(upp_gt, rev(low_gt)), col = gt_col_alpha, border = NA)
    points(maf, mean_gt, pch = 19, cex = .5, type='l', col = gt_col)

  }
  legend("topright", legend = c('aa - 2', 'aA - 1', 'AA - 0', 'Beta'), col = c('blue', 'green', 'red','grey'),
         lty = 1, bty = 'n', pch = c(NA, NA, NA, 19))
  suppressWarnings(par(old_par))
}

plot_phenos(shape = c(.1, .1), factor = c(.1,0))
plot_phenos(shape = c(1, 25), factor = c(.1,0))
