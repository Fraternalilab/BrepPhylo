#' Function to ensure, that a given sequence alignment has valid codons, i.e. a number of bases divisible by 3 without rest.
#'
#' @param inputAlignment Input alignment (FASTA format in \link[bio3d]{bio3d}) with nucleic acid sequences.
#' @param cutEnd If \code{TRUE}, surplus bases will be cut at the end otherwise at the beginning of the nucleotide sequences (default: \code{TRUE}).
#' @param removeSTOPcodons If \code{TRUE}, the first position in STOP codons in the sequences is replaced by a dash (default: \code{TRUE}).
#' @param maxSTOPperSequence If "removeSTOPcodons" is \code{TRUE} and the number of STOP codons for a sequence exceeds this very value then this particular sequence is removed from the alignment (default: 1).
#' @param do.warn If \code{TRUE}, a warning will issued if the aligned sequences do not have a length which is a clean multiple of 3 (default: \code{FALSE}).
#' 
#' @return A list with the two elements "alignment" (the modified alignment) and "report", the latter being a named list with the elements "noCutPositions" (number of positions, that have been cut) and the list "STOPreplaced". List "STOPreplaced" holds identifiers for the sequences and "noSTOPreplaced", the number of STOP codons which have been replaced (only for those where "noSTOPreplaced" is greater than zero).
#' 
#' @examples
#' \dontrun{alignmentValidCodons <- function( inputAlignment = ali )}
#' 
#' @export
alignmentValidCodons <- function( inputAlignment,
                                  cutEnd = TRUE,
                                  removeSTOPcodons = TRUE,
                                  maxSTOPperSequence = 1,
                                  do.warn = FALSE )
{
  codonsNoStop <- NULL
  # Note: The alignment part in a bio3d::fasta object is called "ali" and is a matrix.
  #       The other elements ("id" and "call") will not be affected by a mere sequence cutting.
  
  # check input
  if( ncol( inputAlignment[[ "ali" ]] ) < 3 )
    stop( "The sequences must at least hold a full codon (3 positions)." )
  
  # prepare return object
  listReturn <- list( "alignment" = NA, "report" = list( "noCutPositions" = 0,
                                                         "STOPreplaced" = list(),
                                                         "STOPremoved" = list() ) )
  
  # check, if cutting is necessary
  if( ( ncol( inputAlignment[[ "ali" ]] ) %% 3 ) > 0 )
  {
    noCutPositions <- ncol( inputAlignment[[ "ali" ]] ) %% 3
    if( cutEnd ) {
      inputAlignment[[ "ali" ]] <- inputAlignment[[ "ali" ]][ , 1:( ncol( inputAlignment[[ "ali" ]] ) - noCutPositions ), drop = FALSE ]
    } else {
      inputAlignment[[ "ali" ]] <- inputAlignment[[ "ali" ]][ , ( 1 + noCutPositions ):ncol( inputAlignment[[ "ali" ]] ), drop = FALSE ]
    }
    if( do.warn )
      warning( paste0( c( "Note: The sequence alignment had a length, which was not divisible by 3 without rest, so ", noCutPositions, " were cut off." ) ) )
    listReturn[[ "report" ]][[ "noCutPositions" ]] <- noCutPositions
  }
  
  # check, if STOP codons need to be removed
  if( removeSTOPcodons )
    for( iRow in nrow( inputAlignment[[ "ali" ]] ):1 )
    {
      # STOP codons: TAA, TAG and TGA
      curRow <- inputAlignment[[ "ali" ]][ iRow, , drop = FALSE ]
      curRowString <- paste( as.character( curRow ), collapse = "" )
      codons <- base::substring( text = curRowString,
                                 first = seq( from = 1, to = nchar( curRowString ) - 2, by = 3 ),
                                 last = seq( from = 3, to = nchar( curRowString ), by = 3 ) )
      
      # replace all first positions in STOP codons by '-'
      codonsNoSTOP <- codons
      positions <- which( codonsNoSTOP == "TAA" | codonsNoSTOP == "TAG" | codonsNoSTOP == "TGA" )
      codonsNoSTOP[ positions ] <- paste( '-', substr( x = codonsNoSTOP[ positions ], start = 2, stop = 3 ), sep = "" )
      
      # do update of row and print a warning, in case something had to be replaced
      if( sum( !( codons == codonsNoSTOP ) ) > 0 )
      {
        # check, if there were too many STOP codons in the sequence and remove it in case
        if( sum( !( codons == codonsNoSTOP ) ) > maxSTOPperSequence )
        {
          listReturn[[ "report" ]][[ "STOPremoved" ]][[ rownames( inputAlignment[[ "ali" ]] )[ iRow ] ]] <- sum( !codons == codonsNoSTOP )
          inputAlignment[[ "id" ]] <- inputAlignment[[ "id" ]][ inputAlignment[[ "id" ]] != rownames( inputAlignment[[ "ali" ]] )[ iRow ] ]
          inputAlignment[[ "ali" ]] <- inputAlignment[[ "ali" ]][ -iRow, , drop = FALSE ]
          if( do.warn )
            warning( "Note: Removed sequences as there were too many STOP codons." )
        } else {
          listReturn[[ "report" ]][[ "STOPreplaced" ]][[ rownames( inputAlignment[[ "ali" ]] )[ iRow ] ]] <- sum( !codons == codonsNoSTOP )
          inputAlignment[[ "ali" ]][ iRow, ] <- unlist( strsplit( x = paste( codonsNoSTOP, collapse = "" ), split = "", fixed = TRUE ) )
          if( do.warn )
            warning( paste0( c( "Note: Replaced ", sum( !( codons == codonsNoStop ) ), " first position(s) by '-', as they were (a) STOP codon(s)." ) ) )
        }
      }
    }
  listReturn[[ "alignment" ]] <- inputAlignment
  
  # return the result
  return( listReturn )
}