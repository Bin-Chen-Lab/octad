#' @export
#' @importFrom ExperimentHub ExperimentHub

.eh <- NULL

.onLoad <- function(libname, pkgname) {
  .eh <<- ExperimentHub::ExperimentHub()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to octad package. If you want to run the pipeline on the webserver, please, refer to octad.org")
}