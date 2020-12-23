#' Extract the second anchor from GenomicInteractions object
#' with metadata
#'
#' Creates a GRanges object from the second anchor of a GenomicInteractions
#' object, but extracts all metadata columns.
#' @keywords GenomicRanges GenomicInteractions metadata
#' @importFrom GenomicInteractions anchorTwo
#' @importFrom S4Vectors mcols
#' @importFrom magrittr %>%
#' @param gint A GenomicInteractions object
#' @return gr a GRanges object with metadata 
#' @examples
#'  load(system.file(
#'  "extdata","gs3dk_sample_data.rda",
#'  package = "GS3DKViz"))
#'  anchorTwoWithMetadata(gs3dk_gint)
#' @export
anchorTwoWithMetadata=function(gint) {
  gr=anchorTwo(gint)
  S4Vectors::mcols(gr)=S4Vectors::mcols(gint)  %>% tibble::as_tibble() %>% janitor::clean_names()
  return(gr)
}