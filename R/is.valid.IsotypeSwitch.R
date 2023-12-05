#' Function that checks, whether a phylogenetic relationship in B cells is possible or not according to the genetic arrangement.
#'
#' @param subclassParent The isotype subclass of the possible parent.
#' @param subclassChild The isotype subclass of the possible child.
#' 
#' @description The genetic order in humans is "germline", 'M', 'D', "G3", "G1", "A1", "G2", "G4", 'E' and "A2". Note, that the letters "Ig" are cut off the input, in case they are present.
#'
#' @return Boolean, indicating whether this parent - child relationship is possible (\code{TRUE}) or not (\code{FALSE}). If both subclasses are the same, the result is \code{TRUE}. If the input subclasses are not valid, 'NA' is returned.
#'
#' @examples
#' # FALSE
#' is.valid.IsotypeSwitch( "IgG3", "IgM" )
#' 
#' # TRUE
#' is.valid.IsotypeSwitch( "A1", "A2" )
#' 
is.valid.IsotypeSwitch <- function( subclassParent, subclassChild )
{
  if( grepl( "L|K", subclassParent ) | grepl( "L|K", subclassChild ) ){
    # it is a light chain, accept anyway
    return(TRUE)
  }
  
  # prepare the order vector
  orderSubclasses <- c( "germline", 'M', 'D', "G3", "G1", "A1", "G2", "G4", 'E', "A2" )

  # cut off an "Ig/IGH" substring, if necessary
  subclassParent <- gsub( pattern = "Ig", replacement = "", x = subclassParent )
  subclassChild <- gsub( pattern = "Ig", replacement = "", x = subclassChild )
  subclassParent <- gsub( pattern = "IGH", replacement = "", x = subclassParent )
  subclassChild <- gsub( pattern = "IGH", replacement = "", x = subclassChild )
  
  # check input
  if( !( subclassParent %in% orderSubclasses ) ||
      !( subclassChild %in% orderSubclasses ) )
    return( NA )

  # check relationship and return the result
  if( which( subclassParent == orderSubclasses ) > which( subclassChild == orderSubclasses ) ) {
    return( FALSE )
  } else {
    return( TRUE )
  }
    
}