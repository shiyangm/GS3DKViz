#' Create a GenomicInteractions object
#' with metadata
#'
#' Creates a GenomicInteractions object from a CSV file
#' 
#' @keywords GenomicRanges GenomicInteractions metadata
#' @keywords GenomicRanges GenomicInteractions metadata
#' @importFrom magrittr %>%
#' @importFrom S4Vectors mcols
#' @importFrom rlang .data
#' @param csvfile Filename to be read
#' @return gint a GenomicInteractions object with metadata 
#' @examples
#' readGint(system.file(
#'   "extdata","signficant_genes_genoscan_3D_knock.csv",
#'   package = "GS3DKViz"))
#' @export
readGint=function(csvfile) {
  gint=data.table::fread(csvfile)  %>% 
    #automatically clean column names, converting to standard format.
    janitor::clean_names() %>%
    #split the gene position into start and end
    tidyr::separate(.data$gene_position,c("genestart","geneend"),  "-") %>% 
    #split the promoter position into start and end
    tidyr::separate(.data$promoter_position,c("promoterstart","promoterend"),  "-") %>% 
    #split the enhancer position into start and end
    tidyr::separate(.data$best_enhancer_position,c("enhstart","enhend"),  "-") %>%
    #convert the promoter starts to numeric, assign replacements for missing data.
    dplyr::mutate(promoterstart=as.numeric(.data$promoterstart),chr1=.data$chr,chr2=.data$chr,
                  #replace missing promoter starts with TSS
                  promoterstart=ifelse(is.na(.data$promoterstart),.data$genestart,.data$promoterstart),
                  #replace missing promoter ends with TSS
                  promoterend=ifelse(is.na(.data$promoterend),.data$genestart,.data$promoterend)) %>%  
    #create metadata column that determines whether an interaction is from
    #the ABC model or from GeneHancer.
    dplyr::mutate(enhancer_type=ifelse(grepl(pattern="GH",x=.data$enhancer_id),"GH","ABC")) %>% 
    #rearrange data columns for GenomicInteractions Objects
    dplyr::select(.data$chr1,.data$promoterstart,.data$promoterend,.data$chr2,.data$enhstart,.data$enhend,dplyr::everything()) %>% 
    #rename columns to standard GenomicInteractions format.
    dplyr::rename(start1=.data$promoterstart,end1=.data$promoterend,start2=.data$enhstart,end2=.data$enhend) %>%
    #remove chromosomes from second range (ranges should only be positions).
    dplyr::mutate(start2=gsub("chr19: ","",.data$start2)) %>% 
    #create the GI.
    makeGenomicInteractionsFromDataFrame() %>%
    #remove arcs where the promoter and enhancers overlap.
    .data$.[which(GenomicInteractions::calculateDistances(.data$.)!=0),]
  #with the GenomicInteractions object created, we'll normalize the interaction score
  #and store it in the counts column for plotting.
  S4Vectors::mcols(gint)=gint %>% tibble::as_tibble() %>% janitor::clean_names() %>% 
    dplyr::group_by(.data$enhancer_type) %>% dplyr::mutate(interaction_score=.data$interaction_score/mean(.data$interaction_score,na.rm=T))
  S4Vectors::mcols(gint)$counts<-ifelse(!is.na(S4Vectors::mcols(gint)$interaction_score),
                                        S4Vectors::mcols(gint)$interaction_score,1) 
  return(gint)
}