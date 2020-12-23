#' Creates a GenomicInteractions object from a DataFrame in 
#' the appropriate format 
#' (having columns chr1,start1,end1,chr2,start2,end2,....).
#'
#' Creates a GenomicInteractions object from a DataFrame in
#' the appropriate format 
#' (having columns chr1,start1,end1,chr2,start2,end2,....).
#' @keywords GenomicRanges GenomicInteractions data.frame  tibble data.table
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @param inputdf Dataframe with six cols, with chr1,start1,end1,chr2,start2,end2
#' @param n the number of rows to use from the dataframe.
#' @param includemetadata include metadata (the columns after the first six)
#'  in output?
#' @return output_gint 
#' @examples
#' library(magrittr)
#'  gs3dk_gint <- data.table::fread(system.file(
#'  "extdata","signficant_genes_genoscan_3D_knock.csv",
#'  package = "GS3DKViz")) %>% 
#'  #automatically clean column names, converting to standard format.
#'  janitor::clean_names() %>%
#'    #split the gene position into start and end
#'    tidyr::separate(gene_position,c("genestart","geneend"),  "-") %>% 
#'    #split the promoter position into start and end
#'    tidyr::separate(promoter_position,c("promoterstart","promoterend"),  "-") %>% 
#'    #split the enhancer position into start and end
#'    tidyr::separate(best_enhancer_position,c("enhstart","enhend"),  "-") %>%
#'    #convert the promoter starts to numeric, assign replacements for missing data.
#'    dplyr::mutate(promoterstart=as.numeric(promoterstart),chr1=chr,chr2=chr,
#'                  #replace missing promoter starts with TSS
#'                  promoterstart=ifelse(is.na(promoterstart),genestart,
#'                                       promoterstart),
#'                  #replace missing promoter ends with TSS
#'                  promoterend=ifelse(is.na(promoterend),genestart
#'                                     ,promoterend)) %>%  
#'    #create metadata column that determines whether an interaction is from
#'    #the ABC model or from GeneHancer.
#'    dplyr::mutate(enhancer_type=ifelse(grepl(pattern="GH",x=enhancer_id),
#'                                       "GH","ABC")) %>% 
#'    #rearrange data columns for GenomicInteractions Objects
#'    dplyr::select(chr1,promoterstart,promoterend,chr2,enhstart,enhend,
#'                  dplyr::everything()) %>% 
#'    #rename columns to standard GenomicInteractions format.
#'    dplyr::rename(start1=promoterstart,end1=promoterend,start2=enhstart,
#'                end2=enhend) %>%
#'    #remove chromosomes from second range (ranges should only be positions).
#'    dplyr::mutate(start2=gsub("chr19: ","",start2)) %>% 
#'    #create the GI.
#' makeGenomicInteractionsFromDataFrame() %>% 
#'  .[which(GenomicInteractions::calculateDistances(.)!=0),]
#'    #remove arcs where the promoter and enhancers overlap.
#'    mcols(gs3dk_gint)=gs3dk_gint %>% tibble::as_tibble() %>% janitor::clean_names() %>% 
#'    dplyr::group_by(enhancer_type) %>%
#'    dplyr::mutate(interaction_score=interaction_score/mean(interaction_score,na.rm=T)) 
#'    mcols(gs3dk_gint)$counts<-ifelse(!is.na(mcols(gs3dk_gint)$interaction_score),
#'    mcols(gs3dk_gint)$interaction_score,1) 
#' @export

makeGenomicInteractionsFromDataFrame<-function(inputdf,n=nrow(inputdf),includemetadata=T)
{
 #appropriately renames seqnames column
  if("seqnames1" %in% colnames(inputdf) ){
    colnames(inputdf)[which(colnames(inputdf)=="seqnames1")]<-"chr1"
  }
  if("seqnames2" %in% colnames(inputdf) ){
    colnames(inputdf)[which(colnames(inputdf)=="seqnames2")]<-"chr2"
  }
  if(includemetadata==TRUE)
  {
    #creates genomic ranges objects, links them together in a genomic interactions object, with all metadata columns.
    output_gint<-GenomicInteractions( GenomicRanges::GRanges(inputdf$chr1[1:n],
                                         IRanges(as.numeric(inputdf$start1)[1:n], as.numeric(inputdf$end1)[1:n])),
                                      GenomicRanges::GRanges(inputdf$chr2[1:n],
                                         IRanges(as.numeric(inputdf$start2)[1:n], as.numeric(inputdf$end2)[1:n])),...=as.data.frame(inputdf[,7:ncol(inputdf)]))
  }
  if (includemetadata==FALSE)
  {                                 output_gint<-GenomicInteractions( GenomicRanges::GRanges(inputdf$chr1[1:n],
                                                                         IRanges::IRanges(as.numeric(inputdf$start1)[1:n], as.numeric(inputdf$end1)[1:n])),
                                                                      GenomicRanges::GRanges(inputdf$chr2[1:n],
                                                                         IRanges::IRanges(as.numeric(inputdf$start2)[1:n], as.numeric(inputdf$end2)[1:n])))
  }
  return(output_gint)
}
