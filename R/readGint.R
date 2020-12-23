#' Create a GenomicInteractions object
#' with metadata
#'
#' Creates a GenomicInteractions object from a CSV file
#' 
#' @keywords GenomicRanges GenomicInteractions metadata
#' @keywords GenomicRanges GenomicInteractions metadata
#' @importFrom magrittr `%>%`
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
    tidyr::separate(gene_position,c("genestart","geneend"),  "-") %>% 
    #split the promoter position into start and end
    tidyr::separate(promoter_position,c("promoterstart","promoterend"),  "-") %>% 
    #split the enhancer position into start and end
    tidyr::separate(best_enhancer_position,c("enhstart","enhend"),  "-") %>%
    #convert the promoter starts to numeric, assign replacements for missing data.
    dplyr::mutate(promoterstart=as.numeric(promoterstart),chr1=chr,chr2=chr,
                  #replace missing promoter starts with TSS
                  promoterstart=ifelse(is.na(promoterstart),genestart,promoterstart),
                  #replace missing promoter ends with TSS
                  promoterend=ifelse(is.na(promoterend),genestart,promoterend)) %>%  
    #create metadata column that determines whether an interaction is from
    #the ABC model or from GeneHancer.
    dplyr::mutate(enhancer_type=ifelse(grepl(pattern="GH",x=enhancer_id),"GH","ABC")) %>% 
    #rearrange data columns for GenomicInteractions Objects
    dplyr::select(chr1,promoterstart,promoterend,chr2,enhstart,enhend,dplyr::everything()) %>% 
    #rename columns to standard GenomicInteractions format.
    dplyr::rename(start1=promoterstart,end1=promoterend,start2=enhstart,end2=enhend) %>%
    #remove chromosomes from second range (ranges should only be positions).
    dplyr::mutate(start2=gsub("chr19: ","",start2)) %>% 
    #create the GI.
    makeGenomicInteractionsFromDataFrame() %>%
    #remove arcs where the promoter and enhancers overlap.
    .[which(GenomicInteractions::calculateDistances(.)!=0),]
  #with the GenomicInteractions object created, we'll normalize the interaction score
  #and store it in the counts column for plotting.
  mcols(gint)=gint %>% tibble::as_tibble() %>% janitor::clean_names() %>% 
    dplyr::group_by(enhancer_type) %>% dplyr::mutate(interaction_score=interaction_score/mean(interaction_score,na.rm=T))
  mcols(gint)$counts<-ifelse(!is.na(mcols(gint)$interaction_score),
                             mcols(gint)$interaction_score,1) 
  return(gint)
}