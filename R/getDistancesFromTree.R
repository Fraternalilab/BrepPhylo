#' Function to calculate the distances for a B cell clone lineage from a tree object.
#'
#' @param tree Input tree, type \code{phylo}.
#' @param digits Number of digits to consider (default: 4).
#'
#'
#' @export getDistancesFromTree
#' @importFrom adephylo distTips
getDistancesFromTree <- function( tree,
                               digits = 4 )
{
  return( round( x = as.matrix( adephylo::distTips( x = tree ) ), digits = digits ) )
}