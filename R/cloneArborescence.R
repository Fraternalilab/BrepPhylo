#' Function to calculate the arborescence for a fully connected graph of B cell clone members.
#'
#' @param clone Matrix / data frame holding the clone information.
#' @param cloneID Identifier of clone.
#' @param columnSeqID Name of the column holding the ID of the observation.
#' @param columnSubclass Name of the column holding the isotype subclass (see manual page of function \link{getArborescence} for a list of allowed levels in this column).
#' @param tree Character, filename of the phylogenetic tree.
#' @param outputFolder Folder, where the arborescence data and plot are to be stored.
#' @param plotFormat Either 'png' or 'pdf'. (default: 'png')
#' @param colourSpecification List of colours, a meaningful default is set.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{graph}{\code{igraph} object holding the arborescence tree.}
#'   \item{table}{A data.frame object holding the edgelist of the arborescence tree.}
#'   \item{lineage_plot}{\code{ggplot} object of a lineage plot where sequences are ordered by their distance from the germline and coloured by isotypes. See vignette for more details.}
#'   \item{csr_events}{data.frame showing edges in the tree which switches isotype, and an estimated distance-from-germline of this switching event. This is defaulted as the mean distance-from-germline of the two nodes connected by the said edge.}
#' }
#'
#' @importFrom utils write.csv
#' @importFrom grDevices dev.off pdf png
#' @importFrom graphics legend
#' @importFrom igraph graph_from_edgelist V layout.fruchterman.reingold plot.igraph
#' @importFrom alakazam getPathLengths
#' @importFrom ggmuller get_Muller_df Muller_plot
#'
#' @export cloneArborescence
cloneArborescence <- function( clone,
                               cloneID,
                               columnSeqID,
                               columnSubclass,
                               tree,
                               outputFolder,
                               plotFormat = "png",
                               colourSpecification = list( colours = list( "M" = "#FF0000",
                                                                           "IgA1" = "#9999CF",
                                                                           "IgA2" = "#000088",
                                                                           "IgG1" = "#326600",
                                                                           "IgG2" = "#58B200",
                                                                           "IgG3" = "#8BFF19",
                                                                           "IgG4" = "#CBFF99" ),
                                                           germlineColour = "#AAAAAA" ) )
{
  # check input
  if( !dir.exists( outputFolder ) )
    stop( "Output folder apparently does not exist." )
  if( !all( c( columnSeqID, columnSubclass ) %in% colnames( clone ) ) )
    stop( "Column names specified were not found." )
  if( !file.exists( tree ) )
    stop( "tree file does not exist." )
  if( !plotFormat %in% c("png", "pdf"))
    stop( "plotFormat must be either 'png' or 'pdf'." )
  
  # read in tree
  tree <- phytools::read.newick( path.expand( tree ), quiet = TRUE )
  distances <- getDistancesFromTree( tree = tree )
  germline <- colnames( distances )[ which( grepl( "IG.V", colnames( distances ) ) ) ]
  
  # get observation (sequence) IDs and remove all possible white-spaces
  seq <- trimws( x = as.character( clone[ , columnSeqID ] ), which = "both" )

  # order the clone to match the distance matrix
  # note, that observations that are not matched are automatically removed
  order <- base::match( colnames( distances ), seq )
  order <- order[ !is.na( order ) ]
  clone <- clone[ order, , drop = FALSE ]

  # check, if there are observations left
  if( nrow( clone ) < 1 )
    stop( "After ordering and matching the clone to the distance matrix, no observations are left - please check level match." )

  # prepare the data frame
  # note: it is necessary, to attach the germline here, as it is not part of the original clone
  dfSubclass <- data.frame( germline, "germline", stringsAsFactors = FALSE )
  colnames( dfSubclass ) <- c( columnSeqID, columnSubclass )
  dfSubclass <- rbind( dfSubclass, clone[ , c( columnSeqID, columnSubclass ) ] )
  colnames( dfSubclass ) <- c( "SeqID", "Class" )
  #dfSubclass <- data.frame( SeqID = colnames( distances ),
  #                          Class = as.character( c( as.character( clone[ , columnSubclass ] ), "germline" ) ),
  #                          stringsAsFactors = FALSE )

  # calculate the arborescence
  dfArborescence <- getArborescence( distances = distances, dfSubclass = dfSubclass,
                                     germline = germline)

  if( nrow( dfArborescence ) == 0 ){
    return( list( "graph" = NULL, "table" = NULL,
                  "lineage_plot" = NULL, "csr_events" = NULL) )
  }
  
  # write out the arborescence data
  write.csv( x = dfArborescence, file = paste0( outputFolder, "/clone_", cloneID, ".arbo" ), row.names = FALSE )

  # make a plot
  graph <- igraph::graph_from_edgelist( el = as.matrix( dfArborescence[ , c( "startLabel", "endLabel" ) ] ),
                                        directed = TRUE )
  layout <- igraph::layout.fruchterman.reingold( graph )
  vecVertices <- igraph::V( graph )$name
  dfSubclass <- dfSubclass[ base::match( vecVertices, dfSubclass[ , "SeqID" ] ), ]
  cols <- getColourVectorFromLevels( levels = dfSubclass[ , "Class" ], colourSpecification = colourSpecification )
  if( plotFormat == "png"){
    png( filename = paste0( outputFolder, "/clone_", cloneID, "_arbo.png" ), width = 3600, height = 3600 )
  } else {
    pdf( file= paste0( outputFolder, "/clone_", cloneID, "_arbo.pdf" ), width = 36, height = 36 )
  }
  igraph::plot.igraph( graph, layout = layout, vertex.color = cols, vertex.size = 3, edge.arrow.size = 0.75, vertex.label.cex = 1.75 )
  legend( "top",
          legend = c( "germline", unique( names( colourSpecification$colours ) ) ),
          col = c( colourSpecification$germlineColour, unlist( unique( colourSpecification$colours ) ) ),
          pch = 19,
          cex = 4.5,
          pt.cex = 4.25,
          bty = 'n' )
  dev.off()
  
  # add edge weight to graph
  igraph::E(graph)$weight <- dfArborescence$weight
  
  # JN: plot clone composition by generation
  path_lengths <- alakazam::getPathLengths(graph, root = germline)
  steps <- seq(0, max(path_lengths$distance), by = 0.01)
  steps <- c(steps, max(steps) + 0.01)
  lineage_size <- do.call( "rbind", 
                           lapply( steps , function(x){
     o <- path_lengths[path_lengths$distance <= x, ]
     o$Population <- 1
     o$Generation <- x
     o[ , c( "name", "Population", "Generation" ) ]
  }))
  names( lineage_size )[1] <- "Identity"
  lineage_edgelist <- dfArborescence[ , c( "startLabel", "endLabel" ) ]
  names( lineage_edgelist ) <- c( "Parent", "Identity" )
  muller_df <- ggmuller::get_Muller_df( lineage_edgelist, lineage_size )
  muller_df <- merge( muller_df, clone[, c( columnSeqID, columnSubclass )],
                      by.x = "Identity", by.y = columnSeqID )
  muller_plot <- ggmuller::Muller_plot( muller_df, colour_by = columnSubclass,
                                        conceal_edges = TRUE, xlab = "Distance from root",
                                        palette = unlist( colourSpecification[[ "colours" ]] ) )
  # parse class switching events
  csevents <- getCSEvents(clone, columnSeqID, columnSubclass, 
                          dfArborescence, germline = germline, path_stats = path_lengths)

  # return the arborescence
  return( list( "graph" = graph, "table" = dfArborescence,
                "lineage_plot" = muller_plot, "csr_events" = csevents) )
}
