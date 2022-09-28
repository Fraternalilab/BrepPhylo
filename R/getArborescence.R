#' Function to transform a distance matrix of observations into a list of "Arc"s, a tuple defined in the python script "edmond_arborescence.py"
#'
#' @param distances Either a "dist" or a matrix with the observations.
#' @param dfSubclass A two-column data frame holding:
#' \describe{
#'   \item{SeqID}{The sequence ID as used in the distance object.}
#'   \item{Class}{The isotype subclass of the observations.}
#' }
#' @param germline character, the sequence to be taken as the germline
#' Note, that the germline has subclass type "germline" and observations without a valid subclass must be labelled 'NA' (will be ignored).
#'
#' @return A data frame with the columns "startID" (number of start), "endID" (number of end point), "weight" (weight of the arc), "startLabel" (the labels for the start IDs) and "endID" (the labels for the end IDs) for an arborescence minimum spanning tree.
#'
#' @importFrom reticulate source_python
#' 
#' @export getArborescence
#' 
getArborescence <- function( distances,
                             dfSubclass,
                             germline )
{
  Arc <- NULL
  min_spanning_arborescence <- NULL
  # source the python function
  reticulate::source_python( file = paste0( system.file( package = "BrepPhylo" ), "/python/edmond_arborescence.py" ) )

  # check "distances" and transform into a full matrix, if necessary
  if( length( class( distances ) ) == 1 && class( distances ) == "dist" )
    distances <- as.matrix( distances )

  # check the input
  if( nrow( distances ) != nrow( dfSubclass ) )
    stop( "The number of observations does not match the number of rows in matrix \"matSubclass\"." )

  # add internal index to matrix
  dfSubclass[ , "InternalIndex" ] <- 1:nrow( dfSubclass )

  # initialize list of "Arc"s and fill it
  listAllArcs <- list()
  for( iRow in 1:nrow( distances ) )
  {
    curParent <- dfSubclass[ dfSubclass$SeqID == rownames( distances )[ iRow ], , drop = FALSE ]
    for( iCol in 1:ncol( distances ) )
    {
      curChild <- dfSubclass[ dfSubclass$SeqID == colnames( distances )[ iCol ], , drop = FALSE ]
      relationPossible <- is.valid.IsotypeSwitch( subclassParent = curParent[ , "Class" ],
                                                  subclassChild = curChild[ , "Class" ] )
      if( iRow != iCol )
        if( !is.na( relationPossible ) &&
            relationPossible )
        {
          # if an Arc is possible taking the genetic order into account, add it to the graph
          # note, that the internal numbering must be used in order to ensure unique integer indices
          # the weight of the Arc is given by the distance
          listAllArcs[[ length( listAllArcs ) + 1 ]] <- Arc( curParent$InternalIndex, distances[ iRow, iCol ], curChild$InternalIndex )
        }
    }
  }

  # calculate the arborescence with Edmond's algorithm
  # note, that this assumes the "source" (the germline) is the last observation
  listArborescence <- min_spanning_arborescence( listAllArcs, which( dfSubclass$SeqID == germline ) )

  # reformat the "Arc"s into a data frame
  dfArborescence <- data.frame( startID = unlist( lapply( X = listArborescence, FUN = function( x ) x$head ) ),
                                endID = unlist( lapply( X = listArborescence, FUN = function( x ) x$tail ) ),
                                weight = unlist( lapply( X = listArborescence, FUN = function( x ) x$weight ) ) )

  # replace the internal indices by the appropriate labels
  dfArborescence$startLabel <- dfSubclass[ base::match( dfArborescence$start, dfSubclass$InternalIndex ), "SeqID" ]
  dfArborescence$endLabel <- dfSubclass[ base::match( dfArborescence$end, dfSubclass$InternalIndex ), "SeqID" ]

  return( dfArborescence )
}

# the following tree should work (for tests etc.)
# input <- list( Arc( 0, 17, 1 ), Arc( 0,  1, 2 ), Arc( 0, 19, 3 ),
#                Arc( 0, 16, 4 ), Arc( 0, 16, 5 ), Arc( 0, 12, 6 ),
#                Arc( 1,  3, 2 ), Arc( 1,  3, 3 ), Arc( 1, 11, 4 ),
#                Arc( 1, 10, 4 ), Arc( 1, 12, 6 ), Arc( 2,  3, 1 ),
#                Arc( 2,  4, 3 ), Arc( 2,  8, 4 ), Arc( 2,  8, 5 ),
#                Arc( 2, 11, 6 ), Arc( 3,  3, 1 ), Arc( 3,  4, 2 ),
#                Arc( 3, 12, 4 ), Arc( 3, 18, 5 ), Arc( 3, 14, 6 ),
#                Arc( 4, 11, 1 ), Arc( 4,  8, 2 ), Arc( 4, 12, 3 ),
#                Arc( 4,  6, 5 ), Arc( 4, 10, 6 ), Arc( 5, 10, 1 ),
#                Arc( 5,  8, 2 ), Arc( 5, 11, 3 ), Arc( 5,  6, 4 ),
#                Arc( 5,  4, 6 ), Arc( 6, 12, 1 ), Arc( 6, 11, 2 ),
#                Arc( 6, 14, 3 ), Arc( 6, 10, 4 ), Arc( 6, 11, 5 ) )
# arb <- min_spanning_arborescence( input, 0 )
