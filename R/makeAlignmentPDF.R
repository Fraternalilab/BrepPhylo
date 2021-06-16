#' Function to generate PDF files holding an alignment using the function \link{msaPrettyPrint_extended} and \code{LaTeX}.
#'
#' @param clone Clone dataframe.
#' @param alignment Alignment for the clone (including the germline sequence at the first position).
#' @param cloneIdColumn Name of the column holding the clone IDs.
#' @param outputFolder Path to the output folder.
#' @param colours A named vector of HTML colours for the sequences (including the '#' sign).
#' 
#' @examples
#' \dontrun{makeAlignmentPDF <- function( clone = clone,
#'                               alignment = alignment,
#'                               cloneIdColumn = "Unique_CloneID",
#'                               outputFolder = "~/Desktop",
#'                               colours = colours )}
#'
#' @export
makeAlignmentPDF <- function( clone,
                              alignment,
                              cloneIdColumn,
                              outputFolder,
                              colours )
{
  # use capital letters for LaTeX compatibility
  colours <- toupper( colours )

  # make directory (and remove it, if already present)
  cloneID <- clone[ 1, cloneIdColumn ]
  tempDir = paste0( tempdir(), "/makeAlignmentPDF_clone_", cloneID )
  if( dir.exists( paths = tempDir ) )
    unlink( x = tempDir, recursive = TRUE )
  dir.create( path = tempDir )
  setwd( tempDir )

  # do the definitions for the LaTeX call thoroughly
  stringCodeBase <- "\\shadingmode{diverse}\n\\showruler{1}{top}\n\\showsequencelogo{top}"
  levelsColour <- unique( names( colours ) )
  levelsColour <- levelsColour[ !is.na( levelsColour ) ]
  stringCodeBase <- paste( stringCodeBase, paste0( "\n\\definecolor{colour_germline}{HTML}{", substr( x = colours[ 1 ], start = 2, stop = 7 ), "}\n" ) )
  for( iLevel in 2:length( levelsColour ) )
    stringCodeBase <- paste( stringCodeBase, paste0( "\\definecolor{colour_", iLevel - 1, "}{HTML}{", substr( x = colours[ names( colours ) == levelsColour[ iLevel ] ][ 1 ], start = 2, stop = 7 ), '}' ), sep = "\n" )

  # calculate the required number of PDF pages; do the file path preparation here
  numberParts <- ceiling( ( nrow( alignment ) - 1 ) / 85 )
  archiveFilePaths <- c()

  # do the parts individually
  for( iPart in 1:numberParts )
  {
    stringCode <- stringCodeBase
    selection <- NULL
    stringCode <- paste0( stringCode, "\\namecolor{1}{colour_germline}\n" )
    if( iPart == 1 ) {
      selection <- c( 1, 2:min( nrow( alignment ), 85 ) )
    } else { selection <- c( 1, ( 2 + 84 * ( iPart - 1 ) ):min( nrow( alignment ), ( 85 + 84 * ( iPart - 1 ) ) ) ) }

    # build vectors of indices to be used
    vecSelectedLevels <- names( colours )[ selection ]
    for( iLevel in 2:length( levelsColour ) )
      stringCode <- paste( stringCode, paste0( "\\namecolor{", paste( which( vecSelectedLevels == levelsColour[ iLevel ] ), collapse = ',' ),"}{colour_", iLevel - 1, '}' ), sep = "\n" )

    msaPrettyPrint_extended( x = alignment, output = "pdf", file = paste0( tempDir, "/clone_", cloneID, "_part_", iPart, ".pdf" ),
                             paperWidth = 8.27, paperHeight = 11.69, showLogo = "top",
                             showNames = "left", askForOverwrite = FALSE, psFonts = TRUE,
                             code = stringCode, vecAdditionalHeader = c( "\\usepackage{xcolor}\n" ), subset = selection )

    # move the PDF to the appropriate location
    file.copy( from = paste0( tempDir, "/clone_", cloneID, "_part_", iPart, ".pdf" ),
               to = paste0( outputFolder, "/clone_", cloneID, "_part_", iPart, ".pdf" ) )
  }
}