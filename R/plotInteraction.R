#' Plot a GenomicInteractions object
#' 
#' Creates a Gviz promoter-enhancer plot
#' 
#' @name plotInteraction
#' @keywords GenomicRanges GenomicInteractions metadata Gviz
#' @param gint genomic ranges object on a single chromosome
#' @param chr chromsome
#' @param bounds_offset Beyond the largest and smallest elements, how much extra space in bp should be plotted?
#' @importFrom Gviz IdeogramTrack AnnotationTrack GenomeAxisTrack
#' @importFrom S4Vectors mcols
#' @importFrom GenomicInteractions InteractionTrack
#' @examples
#' #See example at the start of the vignette
#' @export
plotInteraction=function(gint,chr,bounds_offset=1.5e4,main=NULL) {
  itrack <- Gviz::IdeogramTrack(genome = "hg38", chromosome=chr)
  gtrack <- GenomeAxisTrack()
  promoter_gr=anchorOneWithMetadata(gs3dk_gint)
  enhancer_gr=anchorTwoWithMetadata(gs3dk_gint)
  promoterTrack <- AnnotationTrack(promoter_gr, genome="hg38", name="Promoters",
                                   id=S4Vectors::mcols(promoter_gr)$gene_id,  featureAnnotation="id")
  enhancerTrack <- AnnotationTrack(enhancer_gr, genome="hg38", name="Enhancers",
                                   id=S4Vectors::mcols(enhancer_gr)$gene_id,  featureAnnotation="id")
  feature(enhancerTrack)<-S4Vectors::mcols(gint)$enhancer_type
  
  displayPars(promoterTrack) <- list(fill = "olivedrab1", col = NA, 
                                     fontcolor.feature = "black", fontsize=8,
                                     just.group="below",rotation=90,rotation.group=90,rotation.item=90)
  displayPars(enhancerTrack) <- list(fill = "mediumpurple1", col = NA, 
                                     fontcolor.feature = "black", fontsize=10,
                                     just.group="below",rotation.item=90,
                                     collapse=T,mergeGroups=T,showOverplotting=T,groupAnnotation="group",group=S4Vectors::mcols(gint)$enhancer_type)
  displayPars(interaction_track) <- list(fill = "deepskyblue", col = NA, 
                                         fontcolor.feature = "black", fontsize=8,
                                         just.group="below",plot.anchors=T,plot.outside=T,col.outside="lightblue",                                   interaction.measure="counts",
                                         interaction.dimension="height",
                                         col.interactions="black",
                                         plot.trans=T,
                                         fontsize.legend=200
  )
  
  interaction_track <- GenomicInteractions::InteractionTrack(gint, name = "Interaction", chromosome = chr)
  displayPars(interaction_track)=list(col.interactions="black")
  bounds=c(gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(start1) %>% min(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(end1) %>% max(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(start2) %>% min(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(end2) %>% max())
  
  plotTracks(list(itrack,gtrack,interaction_track,promoterTrack,enhancerTrack),
             chromosome=chr,
             from = (min(bounds))-bounds_offset,
             to = (max(bounds))+bounds_offset,
             type = c("b"),
             showSampleNames = TRUE, 
             cex.sampleNames = 0.6,
             cex.main=3,
             fontsize=18,fontsize.item=12,collapse=T,min.width=5,mergeGroups=T
             ,stacking="dense" #comment out to restore stacking
             ,main=main,
             background.title = "black"
  )
  
 return(NULL) 
}