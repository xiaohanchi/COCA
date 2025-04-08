#' Standardize Drug Dosage
#'
#' @param d Dose vectors for all considered doses
#'
#' @return Standardized doses
#' @export
#'
#' @examples
#' dosage_level <- c(300, 200, 0)
#' DoseStandardize(dosage_level)
DoseStandardize <- function(d) {
  return(d / max(d))
}
