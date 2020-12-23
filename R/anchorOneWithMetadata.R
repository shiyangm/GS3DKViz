#' Extract the first anchor from GenomicInteractions object
#' with metadata
#'
#' Creates a GRanges object from the first anchor of a GenomicInteractions
#' object, but extracts all metadata columns.
#' @keywords GenomicRanges GenomicInteractions metadata
#' @importFrom GenomicInteractions anchorOne
#' @importFrom S4Vectors mcols
#' @importFrom magrittr %>%
#' @param gint A GenomicInteractions object
#' @return gr a GRanges object with metadata 
#' @examples
#'  load(system.file(
#'  "extdata","gs3dk_sample_data.rda",
#'  package = "GS3DKViz"))
#'  anchorOneWithMetadata(gs3dk_gint)
#' @export
anchorOneWithMetadata=function(gint) {
  gr=anchorOne(gint)
  S4Vectors::mcols(gr)=S4Vectors::mcols(gint) %>% tibble::as_tibble() %>% janitor::clean_names()
  return(gr)
}