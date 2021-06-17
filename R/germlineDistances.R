#' Function to calculate distance from germline for all observations in a given lineage tree
#'
#' @param tree_file Character, filename of the phylogenetic tree in newick format.
#'
#' @description Read a newick tree to calculate distance-from-germline. Assume germline is in the tree, with name beginning with 'IG.V' where '.' can be any character. 
#' 
#' @return A data.frame with two columns with the following names:
#' \describe{
#'   \item{SeqID}{Character, Identifier of sequence given in the tree file.}
#'   \item{distFromGermline}{Numeric, distance from the identified germline (i.e. total lengths of edges traversed in order to get from the germline to the given sequence.}
#' }
#'
#' @importFrom ape read.tree
#' @importFrom igraph distances
#' @importFrom alakazam phyloToGraph
#'
#' @export getGermlineDistance
getGermlineDistance <- function( tree_file ){
  tree <- ape::read.tree( tree_file )
  if( length(tree$tip.label) < 7 ) return(NULL)
  root <- tree$tip.label[ which(grepl("^IG.V", tree$tip.label)) ]
  dists <- igraph::distances(alakazam::phyloToGraph(tree, germline = root),
                             v = root,
                             to = tree$tip.label[which(tree$tip.label != root)])
  o <- data.frame(t(dists))
  colnames(o) <- "distFromGermline"
  o$SeqID <- rownames(o)
  rownames(o) <- NULL
  o[, c("SeqID", "distFromGermline")]
}

#' Function to calculate percentile of distance-from-germline over all clones in the selection
#'
#' @param dists numeric vector, median distance-from-germline of a set of clones.
#' @param percentiles numeric vector, a list of quantiles to consider. For percentiles, they should be defined between 0 and 1 with step 0.01.
#'
#' @description Internal function used by \code{summariseGermlineDistance}. 
#' 
#' @return numeric vector, percentile in the distance-from-germline distribution over all clones, for each corresponding entry in \code{dist}.
#' 
#' @examples
#' \dontrun{
#' # distFromGermline_table is the output of 'getGermlineDistance'
#' # dist_median is a vector of median distance-from-germline calculated per clone
#' dist_percentiles <- getDistPercentile( dist_median, 
#'   percentiles = quantile(distFromGermline_table$distFromGermline, 
#'                          probs = seq(0, 1, by = 0.01)) )
#' }
#'
getDistPercentile <- function(dists, percentiles)
{
  sapply(dists, function(d){
    min(which(d < percentiles))/100
  })
}

#' Function to summarise distance-from-germline by clones
#'
#' @param distFromGermline data.frame, output of \code{getGermlineDistance} and annotated with extra columns detailing metadata of samples.
#' @param dist_column Character, column name in \code{distFromGermline} which holds the distance-from-germline metrics.
#' @param cloneID_column Character, column name in \code{distFromGermline} which holds the clone ID to which the given sequence belongs.
#' @param summarise_variables Character vector, names of columns in \code{distFromGermline} used to partition the sequences into subsets for which summary is sought
#'
#' @description The function first summarises for each clone (optionally with other metadata columns present in the input \code{distFromGermline} data frame) the median distance-from-germline of all sequences in the clone.
#' It then used this median distance to order the clones (from lowest - i.e. fewest mutations - to highest) and treat this as an 'expected' order of the clone in terms of mutational level, and calculate for each clone its 'actual' percentile of its median distance over all other clones.
#' The discrepancy between the 'actual' and 'expected' percentiles returned by this function can subsequently be used to plot curves showing overall mutational level of a repertoire (see Example), or to seek for a quantification using the Germline Likeness metric (see \code{?getGermlineLikeness}).
#' 
#' @return A data.frame with each row corresponding to one clone, containing, in addition to metadata columns given in \code{summarise_variables}, the following columns:
#' \describe{
#'   \item{dist_median}{Numeric, between 0 to 1, the percentile of the actual median distance-from-germline of the given clone.}
#'   \item{clone_order}{Numeric, between 0 to 1, the percentile of the given clone in the distribution where all clones are ordered by their median distance.}
#' }
#'
#' @import plyr
#' @importFrom stats median quantile
#'
#' @export summariseGermlineDistance
summariseGermlineDistance <- function(distFromGermline, dist_column = "dist", 
                                      cloneID_column = "CloneID",
                                      summarise_variables){
  if( !cloneID_column %in% summarise_variables ){
    stop("'cloneID_column' should be present in 'summarise_variables' otherwise the germline distance summary will not make sense.")
  }
  distFromGermline_median <- ddply(distFromGermline,
                                   .variables = summarise_variables,
                                   colwise(median, as.quoted(dist_column), na.rm = TRUE))
  colnames( distFromGermline_median )[ ncol( distFromGermline_median ) ] <- "dist_median"
  distFromGermline_median$dist_median <- getDistPercentile(
    distFromGermline_median$dist_median,
    quantile(distFromGermline[, dist_column], probs = seq(0.01, 1, by = 0.01))
  )
  split_vars <- lapply(summarise_variables[ -which(summarise_variables == cloneID_column) ],
                       function(x) distFromGermline_median[, x])
  distFromGermline_median <- split(distFromGermline_median,
                                   f = split_vars, drop = TRUE)
  distFromGermline_median <- lapply(distFromGermline_median, function(tb){
    tb <- tb[order(tb$dist_median), ]
    tb$clone_order <- (1:nrow(tb)) / nrow(tb)
    tb
  })
  distFromGermline_median <- do.call("rbind", distFromGermline_median)
  distFromGermline_median
}

#' Function to calculate germline likeness of a pre-selected subset of clones
#'
#' @param tb data.frame, subset of the output from \code{summariseGermlineDistance} corresponding to a set of clones of interest.
#'
#' @description Internal function used by \code{getGermlineLikeness}. 
#' 
#' @return numeric, Germline Likeness metric (between 0 and 1) summarising the given set of clones.
#' 
#'
GermlineLikeness <- function(tb)
{
  dDist <- c(diff(tb$dist_median), 0)
  dClone <- c(diff(tb$clone_order), 0)
  # classical empirical area-under-curve (AUC) calculation
  sum(tb$clone_order * dDist) + sum(dDist * dClone)/2
}

#' Function to summarise Germline Likeness for repertoires
#'
#' @param distFromGermline_summary data.frame, output of \code{summariseGermlineDistance}. With metadata annotation including columns given in \code{summarise_variables}.
#' @param metadata_table data.frame with meta data, 1 row per sample/donor as appropriate (see Example).
#' @param summarise_variables Character vector, names of columns in \code{distFromGermline} used to partition the sequences into subsets for which summary is sought
#'
#' @description The function loops through each row in \code{metadata_table} and select the subset of lineages from \code{distFromGermline_summary} relevant to the given metadata combination,
#' and calculate the Germline Likeness metric for this given set of lineages.
#' 
#' @return A data.frame identical to \code{metadata_table} with the additional column "GermlineLikeness" holding the Germline Likeness metric (numeric) of lineages matching the given combination of metadata attributes
#'
#' @export getGermlineLikeness
getGermlineLikeness <- function(distFromGermline_summary, metadata_table, 
                                summarise_variables){
  if( !is.data.frame ( distFromGermline_summary ) ){
    stop("'distFromGermline_summary' should be a data.frame returned by the function 'summariseGermlineDistance'.")
  }
  if( !is.data.frame ( metadata_table ) ){
    stop("'metadata_table' should be a data.frame containing metadata, 1 row per donor/sample.")
  }
  if( any( !summarise_variables %in% colnames(metadata_table)) ){
    stop("Some of the variables in 'summarise_variables' are absent from the column names of 'metadata_table'. Check whether you have indicated the correct columns in 'summarise_variables'?")
  }
  if( any( !summarise_variables %in% colnames(distFromGermline_summary)) ){
    stop("Some of the variables in 'summarise_variables' are absent from the column names of 'distFromGermline_summary'. Check whether you have indicated the correct columns in 'summarise_variables'?")
  }
  summarise_variables_colID <- sapply(summarise_variables, function(x){
    return( which( colnames( distFromGermline_summary ) == x ) ) 
  })
  metadata_colID <- sapply(summarise_variables, function(x){
    return( which( colnames( metadata_table ) == x ) ) 
  })
  # replace NA and NULL with empty string in each relevant distFromGermline_summary column
  for(y in summarise_variables_colID ){
    distFromGermline_summary[, y] <- replace(distFromGermline_summary[, y],
                                             which(is.na(distFromGermline_summary[, y])),
                                             "") 
    distFromGermline_summary[, y] <- replace(distFromGermline_summary[, y],
                                             which(is.null(distFromGermline_summary[, y])),
                                             "") 
  }
  metadata_table$GermlineLikeness <- apply(
    metadata_table, MARGIN = 1, function(x){
      tb <- distFromGermline_summary
      # filter tb for relevant entry specified by the given row in metadata_table
      for( i in 1:length(summarise_variables_colID) ){
        metadata_col <- metadata_colID[ i ]
        distSummary_col <- summarise_variables_colID[ i ]
        # if entry in metadata_table is NA or NULL replace with empty string
        if( is.na(x[metadata_col]) ) x[metadata_col] <- ""
        if( is.null(x[metadata_col]) ) x[metadata_col] <- ""
        tb <- tb[ which( tb[, distSummary_col] == x[metadata_col] ), ]
      }
      GermlineLikeness(tb)
    }
  )
  metadata_table
}
