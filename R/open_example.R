#' Open Example R Script
#'
#' Opens an example R script in the editor.
#'
#' @param name File name in `inst/examples/`
#'
#' @export
#'
#' @examples \donttest{
#' open_example("MTD-Ind.R")
#'
#' open_example("OBD-Ind.R")
#'
#' open_example("OBD-Pool.R")
#' }
#'
open_example <- function(name = "example_code.R") {
  path <- system.file("examples", name, package = "COCA")
  if (path == "") stop("File not found.")
  file.edit(path)
}

