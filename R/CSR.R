#' Function to summarise details of different class-switch recombination events
#'
#' @param tree_stats data.frame, output of \code{getSummaryFromBatch} (specifically the element \code{csr_events} from the output list), and annotated with extra columns detailing metadata of samples.
#' @param dist_column Character, column name in \code{getSummaryFromBatch} which holds the distance-from-germline metrics. (default: 'dist')
#' @param cloneID_column Character, column name in \code{getSummaryFromBatch} which holds the clone ID to which the given sequence belongs. (Default: 'CloneID')
#' @param summarise_variables Character vector, names of columns in \code{getSummaryFromBatch} used to partition the sequences into subsets for which summary is sought.
#'
#' @description The function summarise the CSR details from \code{getSummaryFromBatch}.
#' 
#' @return A data.frame with each row corresponding to one specific CSR switch (i.e. distinct combination of start- and end-points of the switch) in different subset of the data defined by the \code{summarise_variables} column(s). The following summaries were sought and stored in separate columns:
#' \describe{
#'   \item{nEvent}{Numeric, total number of switches involving the two named isotypes in the dataset subset in question.}
#'   \item{meanDist_from_germline}{Numeric, mean distance-from-germline of CSR events involving the two named isotypes.}
#'   \item{nClone}{Numeric, total number of clonotypes in the dataset subset (should be the same across all rows for the same subset of data).}
#'   \item{nClone_events}{Numeric, between 0 and 1, the proportion of clones with the named CSR switch observed.}
#' }
#'
#' @examples 
#' \dontrun{
#' # we provide doBatchCloneAnalysis results for data from the 
#' # two individuals included here in this vignette
#' cv233 <- system.file( "extdata/CV233_batchAnalysis.rds", package = "BrepPhylo")
#' cv325 <- system.file( "extdata/CV325_batchAnalysis.rds", package = "BrepPhylo")
#' cv233 <- readRDS( cv233 )
#' cv325 <- readRDS( cv325 )
#' 
#' # the getSummaryFromBatch function flattens the output to two
#' # elements: a table containing  distance-from-germline (demonstrated 
#' # above) and another table of all CSR events observed in lineages
#' cv233 <- getSummaryFromBatch( cv233 )
#' cv233 <- cv233$csr_events
#' cv325 <- getSummaryFromBatch( cv325 )
#' cv325 <- cv325$csr_events
#' # just add an extra column with PatientID
#' cv233$PatientID <- "CV233"
#' cv325$PatientID <- "CV325"
#' # we can combine the two tables
#' csr_details <- rbind( cv233, cv325 )
#' csr_summary <- summariseCSR( csr_details, 
#'                              dist_column = "distFromGermline", 
#'                              cloneID_column = "CloneID",
#'                              summarise_variables = "PatientID")
#' }
#'
#' @import plyr
#'
#' @export summariseCSR
summariseCSR <- function(tree_stats, cloneID_column = "CloneID", 
                         dist_column = "dist", 
                         summarise_variables){
  if( !is.data.frame( tree_stats ) )
    stop( "'tree_stats' should be a data.frame.")
  if( any( !c("startIsotype", "endIsotype") %in% colnames( tree_stats ) ) )
    stop( "'tree_stats' appear not to be in the format of csr_events table as expected. Are you sure this is taken from getSummaryFromBatch()?" )
  if( any( !c(cloneID_column, dist_column) %in% colnames( tree_stats ) ) )
    stop( "Check whether both of 'cloneID_column' and 'dist_column' is found in the column names of 'tree_stats'.")
  if( any( !summarise_variables %in% colnames( tree_stats ) ) )
    stop( "Check whether all variables indicated in 'summarise_variables' are found in the column names of 'tree_stats'.")
  summarise_variables <- summarise_variables[ which( !summarise_variables %in%
                                                       c("startIsotype", "endIsotype")) ]
  csr_summary <- do.call("ddply",
                         list(tree_stats, 
                              c(summarise_variables, "startIsotype", "endIsotype"), 
                              summarize, 
                              nEvent = call("length", 
                                            call("unique", as.symbol( cloneID_column ) )),
                              meanDist_from_germline = call("mean", as.symbol( dist_column ), na.rm = TRUE )))
  # colwise(median, as.quoted(dist_column), na.rm = TRUE)
  csr_summary <- merge(csr_summary,
                       do.call("ddply",
                               list(tree_stats, 
                                    c(summarise_variables), 
                                    summarize, 
                                    nClone = call("length", 
                                                  call("unique", as.symbol( cloneID_column ) )))),
                           by = summarise_variables )
  csr_summary$csr_frequency <- csr_summary$nEvent / csr_summary$nClone
  csr_summary
}

#' Function to plot class-switch recombination summaries of repertoires
#'
#' @param csr_summary data.frame, output of \code{summariseCSR}.
#'
#' @description The function plot the CSR summary from \code{summariseCSR}. Remember to set factors for columns represented in the plots to alter the order of these labels on axes of the plot (see vignette & example below).
#' 
#' @return A \code{ggplot2} plot object. You can further manipulate it, e.g. plot data from different subsets in separate panels (with \code{facet_wrap()} in \code{ggplot2}) etc.
#'
#' @import ggplot2
#'
#' @examples 
#' \dontrun{
#' # first set the order of the isotypes so that they reflect
#' # the actual physical order on the genome (hence the 
#' # order possible in CSR)
#' csr_summary$startIsotype <- factor(csr_summary$startIsotype,
#'                                    levels = c("IgM", "IgG3", "IgG1", "IgA1",
#'                                               "IgG2", "IgG4", "IgE", "IgA2"),
#'                                    labels = c("M", "G3", "G1", "A1", "G2",
#'                                               "G4", "E", "A2"))
#' csr_summary$endIsotype <- factor(csr_summary$endIsotype,
#'                                  levels = c("IgM", "IgG3", "IgG1", "IgA1",
#'                                             "IgG2", "IgG4", "IgE", "IgA2"),
#'                                  labels = c("M", "G3", "G1", "A1", "G2",
#'                                             "G4", "E", "A2"))
#' 
#' # plot it; it uses ggplot2 so you can extend this with ggplot2 functions
#' # here separate into different panels using the column 'PatientID'
#' library(ggplot2)
#' plotCSRsummary( csr_summary ) + facet_wrap(~ PatientID)
#' 
#' }
#'
#' @export plotCSRsummary
plotCSRsummary <- function(csr_summary){
  ggplot(csr_summary, 
         aes(x = endIsotype, y = startIsotype,
             fill = meanDist_from_germline, size = csr_frequency)) +
    geom_point(pch = 21) + cowplot::theme_cowplot() +
    scale_size_continuous(range = c(0, 10), labels = scales::percent,
                          name = "% clones with\nCSR events") +
    scale_fill_viridis_c(limits = c(0, 0.2),
                         name = "mean distance\nfrom germline") +
    xlab("To") + ylab("From")
}