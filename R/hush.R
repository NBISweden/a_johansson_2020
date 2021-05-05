#' @title Silence too verbose code
#' @author Marcin Kierczak <marcin.kierczak__INSERT_AT__scilifelab.se>
#' @param  code the code to be silenced
#' @return silenced result of code execution
#'
hush <- function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}
