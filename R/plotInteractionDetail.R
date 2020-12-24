#' Plot a GenomicInteractions object with details for Enhancers
#' 
#' Creates a Gviz promoter-enhancer plot
#' 
#' @name plotInteraction
#' @keywords GenomicRanges GenomicInteractions metadata Gviz
#' @param gint genomic ranges object on a single chromosome
#' @param extra_prom_gr A Genomic Ranges
#' object containing extra promoters that are disconnected from interactions
#' @param chr chromsome
#' @param bounds_offset Beyond the largest and smallest elements, how much extra space in bp should be plotted?
#' @param main character string for the title of the plot
#' @importFrom Gviz IdeogramTrack AnnotationTrack GenomeAxisTrack displayPars
#' @importFrom Gviz feature feature<- plotTracks
#' @importFrom S4Vectors mcols
#' @importFrom GenomicInteractions InteractionTrack
#' @importFrom rlang .data
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter pull
#' @examples
#' #See example at the start of the vignette
#' @export
plotInteractionDetail=function(gint,chr,bounds_offset=1.5e4,main=NULL) {
  promoter_gr=anchorOneWithMetadata(gint)
  enhancer_gr=anchorTwoWithMetadata(gint)
  annotateInteractions(GIObject = gint,annotations =  list(promoter=promoter_gr,
                                    enhancer=enhancer_gr),id.col="gene_id")
  knockoff_genes=S4Vectors::mcols(gint) %>% tibble::as_tibble() %>%
    dplyr::filter(identified_by_knockoff_or_not=="Yes") %>% dplyr::pull(gene_id)
  gs3dk_gint@regions@elementMetadata$node.class=(S4Vectors::mcols(regions(
    gs3dk_gint)) %>% 
  tibble::as_tibble() %>%
  dplyr::mutate(knockoff=ifelse((promoter.id %in% knockoff_genes | 
  enhancer.id %in% knockoff_genes),"knockoff_detected","knockoff_removed"),
  node.class=knockoff) %>%
    dplyr::pull(knockoff))
  interaction_track <- InteractionTrack(gint,
  name = "Interaction",
  chromosome = as.character(gs3dk_gint %>% as.data.frame() %>%
      janitor::clean_names() %>% dplyr::filter(trait=="AD") %>%
        dplyr::pull(seqnames1)))
  displayPars(interaction_track)=list(col.interactions="black",
    col.interaction.types=c('knockoff_detected-knockoff_detected'='blue',
        'knockoff_detected-knockoff_removed'='blue',
        'knockoff_removed-knockoff_removed'='black'))
  itrack <- Gviz::IdeogramTrack(genome = "hg38", chromosome=chr)
  gtrack <- GenomeAxisTrack()
  promoterTrack <- AnnotationTrack(promoter_gr, genome="hg38", name="Promoters",
                                   id=S4Vectors::mcols(promoter_gr)$gene_id,
                                   featureAnnotation="id")
  enhancerTrack <- AnnotationTrack(enhancer_gr, genome="hg38", name="Enhancers",
                                   id=S4Vectors::mcols(enhancer_gr)$gene_id,
                                   featureAnnotation="id")
  feature(enhancerTrack)<-S4Vectors::mcols(gint)$enhancer_type
  
  Gviz::displayPars(promoterTrack) <- list(fill = "olivedrab1", col = NA, 
                                     fontcolor.feature = "black", fontsize=8,
                                     just.group="below",rotation=90,rotation.group=90,rotation.item=90)
  Gviz::displayPars(enhancerTrack) <- list(fill = "mediumpurple1", col = NA, 
                                     fontcolor.feature = "black", fontsize=10,
                                     just.group="below",rotation.item=90,
                                     collapse=T,mergeGroups=T,showOverplotting=T,groupAnnotation="group",group=S4Vectors::mcols(gint)$enhancer_type)
  selFun <- function(identifier, start, end, track, GdObject, ...){
    gcount <- table(group(GdObject))
    pxRange <- Gviz:::.pxResolution(min.width = 50, coord = "x")
    return((end - start) < pxRange && gcount[identifier] == 1 && (distanceToNearest(GdObject[group(GdObject) == identifier]@range,GdObject[group(GdObject) != identifier]@range)@elementMetadata$distance)<50
    )
  }
  detFun <- function(identifier, GdObject.original, ...){
    plotTracks(list(GenomeAxisTrack(scale = 0.3, size = 0.2, cex = 0.7), 
                    GdObject.original[group(GdObject.original) == identifier]),         add = TRUE, showTitle = FALSE)
  }
  deTrackEnh2 <- AnnotationTrack(name = "Enhancers",enhancer_gr, fun = detFun, 
                                 selectFun = selFun,
                                 groupDetails = TRUE, details.size = 0.5, 
                                 detailsConnector.cex = 0.5, 
                                 detailsConnector.lty = "dotted",
                                 shape = c("smallArrow", "arrow"), 
                                 groupAnnotation = "group",
                                 id=mcols(enhancer_gr)$gene_id, 
                                 group=subjectHits(findOverlaps(enhancer_gr, reduce(enhancer_gr))),featureAnnotation="id",stacking="hide",
  )
  displayPars(deTrackEnh2) <- list(fill = "mediumpurple1", col = NA, 
                                   fontcolor.feature = "black", fontsize=10,
                                   just.group="below",rotation.item=90,
                                   collapse=T,mergeGroups=T,showOverplotting=T,groupAnnotation="group",    group=mcols(gs3dk_gint)$enhancer_type,feature=mcols(gs3dk_gint)$enhancer_type,
                                   groupDetails=T,detailsConnector.lty="solid",
                                   detailsConnector.col="black",
                                   detailsBorder.col="black",
                                   detailsBorder.lty="solid",showId=F
                                   ,min.distance=5,min.width=5,title="Enhancers") 
  combined_prom_gr=c(promoter_gr,extra_prom_gr)
  CombinedPromoterTrack <- AnnotationTrack(combined_prom_gr, genome="hg38", name="Promoters",
                                           id=mcols(combined_prom_gr)$gene_id,  featureAnnotation="id",groupAnnotation="feature")
  displayPars(CombinedPromoterTrack)=list(fill = "olivedrab1", col = NA, 
                                          fontcolor.feature = "black", fontsize=8,                 just.group="below",rotation=90,rotation.group=90,rotation.item=90,
                                          min.width=1,min.distance=5)
  feature(CombinedPromoterTrack)=c(rep("connected",length(promoter_gr)),rep("unconnected",length(extra_prom_gr)))
  
  # interaction_track <- GenomicInteractions::InteractionTrack(gint, name = "Interaction", chromosome = chr)
  # Gviz::displayPars(interaction_track) <- list(fill = "deepskyblue", col = NA, 
  #                                        fontcolor.feature = "black", fontsize=8,
  #                                        just.group="below",plot.anchors=T,plot.outside=T,col.outside="lightblue",interaction.measure="counts",
  #                                        interaction.dimension="height",
  #                                        col.interactions="black",
  #                                        plot.trans=T,
  #                                        fontsize.legend=200
  # )
  # 
  # Gviz::displayPars(interaction_track)=list(col.interactions="black")
  bounds=c(gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(.data$start1) %>% min(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(.data$end1) %>% max(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(.data$start2) %>% min(),
           gint %>% as.data.frame() %>% janitor::clean_names()  %>% dplyr::pull(.data$end2) %>% max())
  
  Gviz::plotTracks(list(itrack,gtrack,interaction_track,CombinedPromoterTrack,
                        deTrackEnh2),
             chromosome=chr,
             from = (min(bounds))-bounds_offset,
             to = (max(bounds))+bounds_offset,
             type = c("b"),
             showSampleNames = TRUE, 
             cex.sampleNames = 0.6,
             cex.main=3,
             fontsize=18,fontsize.item=7,collapse=T,min.width=5,mergeGroups=T
             ,stacking="dense" #comment out to restore stacking
             ,main=main,
             background.title = "black"
  )
  
 return(NULL) 
}