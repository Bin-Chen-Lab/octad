#' @export
#' @import ExperimentHub
#.onLoad <- function(libname, pkgname) {
#  fl <- system.file("extdata", "metadata.csv", package=pkgname)
#  titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
#  ExperimentHub::createHubAccessors(pkgname, titles)
#}