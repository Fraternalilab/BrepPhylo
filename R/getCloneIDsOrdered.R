#' Function to analyse B cell clone lineages, by producing an alignment and infering a phylogenetic tree.
#'
#' @param input Matrix / data frame with values.
#' @param cloneIdColumn Name of the column holding the clone IDs.
#' @param decreasing If \code{TRUE}, the clones will be ordered according to ever decreasing sizes.
#' 
#' @examples
#' \dontrun{cloneIDdf <- getCloneIDsOrdered( input = inputDF, cloneIdColumn = "Unique_CloneID" )}
#'
#' @export
getCloneIDsOrdered <- function( input,
                                cloneIdColumn,
                                decreasing = TRUE )
{
  # check input
  if( !( cloneIdColumn %in% colnames( input ) ) )
    stop( "The column specified by \"cloneIdColumn\" was not found." )

  # get all clone IDs and initialize return object
  unique_cloneIDs <- unique( input[ , cloneIdColumn ] )
  dfReturn <- data.frame( cloneIdColumn = rep( NA, times = length( unique_cloneIDs ) ),
                          "size" = rep( NA, times = length( unique_cloneIDs ) ), stringsAsFactors = FALSE )

  # fill the matrix
  buf <- as.matrix( table( input[ , cloneIdColumn ] ) )
  dfReturn[ , 1 ] <- rownames( buf )
  dfReturn[ , 2 ] <- as.numeric( buf[ , 1 ] )

  # order the result
  dfReturn <- dfReturn[ order( dfReturn[ , 2 ], decreasing = decreasing ), ]

  # return the object
  return( dfReturn )
}