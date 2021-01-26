##' @title Plot density of markers in a region
##' @author Marcin Kierczak
##' @description Plot density of markers in non-overlaping intervals of given size
##' @param coords - coordinates of all markers
##' @param w - window (interval) size, e.g. 1000 for 1kb
##' @param plot - a logical indicating whether to plot density
##' @return a vector of densities
plot_mrk_density <- function(coords, w = 1000, plot = F) {
  # Compute marker density along the region
  # w <- 1000 # 1kb
  starts <- seq(min(coords), max(coords), by = (w+1))
  stops <- starts + w
  if (rev(stops)[1] < max(coords)) {
    starts <- c(starts, rev(starts)[1] + 1)
    stops <- c(stops, max(coords))
  } else if (rev(stops)[1] > max(coords)) {
    stops[length(stops)] <- max(coords)
  }
  lengths <- stops - starts
  dens <- rep(NA, times = length(starts))
  for (i in 1:length(starts)) {
    dens[i] <- sum(which(coords >= starts[i] & coords <= stops[i])) / lengths[i]
  }
  if (plot) {
    plot(starts/w, dens, type='h', las=1, cex=.5, pch=19, cex.axis=.7, ylab='Marker density')
    points(starts/w, dens, pch=19, cex=.5)
  }
  return(dens)
}
