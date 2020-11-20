#' @title Plots effect sizes distribution dependent on minor allele frequency.
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @param shapes_rare - params for Beta distribution for the rare variants
#' @param shapes_common - params for Beta distribution for the common variants
#' @param thr - threshold, markers with maf below or equal to it are considered rare
#' @return NULL
#' @export
#'
plot_effects <- function(shapes_rare, shapes_common, thr) {
  x <- seq(from = 0, to = .5, by = 0.001)
  effect_sizes_rare <- dbeta(x, shapes_rare[1], shapes_rare[2])
  effect_sizes_common <- dbeta(x, shapes_common[1], shapes_common[2])
  effect_sizes <- c(effect_sizes_rare[which(x <= thr)], effect_sizes_common[which(x > thr)])
  data <- tibble(maf = x, rare = effect_sizes_rare, common = effect_sizes_common, joint = effect_sizes)
  effects <- plot_effects(beta_params_rare, beta_params_common, thr_common_rare)
  effects %>%
    pivot_longer(c(rare, common, joint), names_to = 'type', values_to = 'effect size') %>%
    ggplot(mapping = aes(x = maf, y = `effect size`)) +
    geom_line() +
    geom_vline(xintercept = 0.01, col='red') +
    facet_grid(. ~ type)
}
