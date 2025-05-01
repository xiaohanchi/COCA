#' Open Example R Script
#'
#' Opens an example R script in the editor.
#'
#' @param name File name in `inst/examples/`
#'
#' @export
#'
#' @examples \donttest{
#' open_example("alt_designs.R")
#'
#' }
#'
open_example <- function(name = "example_code.R") {
  path <- system.file("examples", name, package = "COCA")
  if (path == "") stop("File not found.")
  file.edit(path)
}

