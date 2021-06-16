#' Function to add a clone ID and its size to a file.
#' @param folder The path to the folder, where the file lies (filename is always "clones_completed.csv").
#' @param cloneID The clone ID to be added.
#' @param size The size of the clone.
#' @importFrom utils read.csv write.csv
#' @export
addCloneToCompleted <- function( folder, cloneID, size )
{
  # check, if file exists
  path <- paste0( folder, "/clones_completed.csv" )
  matClones <- matrix( c( NA, NA ), nrow = 1 )[ -1, , drop = FALSE ]
  colnames( matClones ) <- c( "cloneID", "size" )
  if( file.exists( path ) )
    matClones <- read.csv( file = path, header = TRUE, stringsAsFactors = FALSE )

  # add the observation, if not present already
  if( nrow( matClones ) < 1 ||
      all( matClones[ , "cloneID" ] != cloneID ) )
    matClones <- rbind( matClones, c( cloneID, size ) )

  # write the result
  write.csv( x = matClones, file = path, row.names = FALSE )
}


#' Function to generate a vector of colours from a vector of levels and a colour specification list
#' @param levels Vector of levels.
#' @param colourSpecification List that specifies which levels are to be represented by which respective colour. Make sure, that the definition matches the levels in parameter \code{levels}.
#' @return A vector of colours in the order of the input levels.
#' @export
getColourVectorFromLevels <- function( levels,
                                       colourSpecification )
{
  # set all colours, except for the unknown ones and the germline
  listColours <- colourSpecification[[ "colours" ]][ as.character( levels ) ]
  
  # set all unknown levels to be black and to be of type "UNK"
  isUNK <- unlist( lapply( X = listColours, FUN = function( x ) is.null( x ) ) )
  listColours[ isUNK ] <- "#000000"
  vecColours <- unlist( listColours )
  names( vecColours )[ is.na( names( vecColours ) ) ] <- "UNK"
  
  # set the germline colour and level
  names( vecColours )[ which( as.character( levels ) == "germline" ) ] <- "germline"
  vecColours[ which( as.character( levels ) == "germline" ) ] <- colourSpecification[[ "germlineColour" ]]
  
  # return the result
  return( vecColours )
}

#' Function to print the current date and the time.
#' @return A string with the current datetime in a specific format.
curDateTime <- function()
{
  return( format( x = Sys.time(), format = "%Y-%m-%d %X" ) )
}

#' Function to parse class-switching events from arborescence tree analysis
#'
#' @param clone data.frame of clone information.
#' @param columnSeqID column name holding Sequence IDs.
#' @param columnSubclass column name holding subclasses information.
#' @param dfArborescence data.frame of arborescence tree edge list.
#' @param germline Root ID for germline.
#' @param path_stats data.frame of distances from germline calculated for the arborescence tree.
#' 
#' @description This function parses edges in the arborescence tree that switches isotype/subisotype, 
#' and offer a distance-from-germline estimate of this switching event, as the mean distance-from-germline 
#' of the two ends of each of such edge.
#'
#' @return A data.frame. 
#'
getCSEvents <- function(clone, 
                        columnSeqID, 
                        columnSubclass, 
                        dfArborescence, 
                        germline,
                        path_stats)
{
  dfArborescence <- merge( dfArborescence, 
                           clone[, c( columnSeqID, columnSubclass )],
                           by.x = "startLabel", by.y = columnSeqID, 
                           all.x = TRUE, all.y = FALSE, sort = FALSE)
  names( dfArborescence )[ ncol( dfArborescence ) ] <- "startIsotype"
  dfArborescence <- merge( dfArborescence, 
                           clone[, c( columnSeqID, columnSubclass )],
                           by.x = "endLabel", by.y = columnSeqID, 
                           all.x = TRUE, all.y = FALSE, sort = FALSE)
  names( dfArborescence )[ ncol( dfArborescence ) ] <- "endIsotype"
  # dfArborescence[ dfArborescence$startLabel == germline, "startIsotype"] <- "germline"
  dfArborescence <- dfArborescence[ !is.na( dfArborescence$startIsotype ) &
                                      !is.na( dfArborescence$endIsotype ) &
                                      dfArborescence$startIsotype != dfArborescence$endIsotype, ]
  if( nrow(dfArborescence) > 0 ){
    dfArborescence <- merge( dfArborescence,
                             path_stats[, c("name", "distance")],
                             by.x = "startLabel", by.y = "name",
                             all.x = TRUE, all.y = FALSE, sort = FALSE)
    names( dfArborescence )[ ncol(dfArborescence) ] <- "startDistance"
    dfArborescence <- merge( dfArborescence,
                             path_stats[, c("name", "distance")],
                             by.x = "endLabel", by.y = "name",
                             all.x = TRUE, all.y = FALSE, sort = FALSE)
    names( dfArborescence )[ ncol(dfArborescence) ] <- "endDistance"
    dfArborescence$distFromGermline <- apply( dfArborescence[, c("startDistance",
                                                                 "endDistance")],
                                              MARGIN = 1, mean, na.rm = TRUE)
    return( dfArborescence[, c( "startLabel", "endLabel", 
                                "startIsotype", "endIsotype",
                                "distFromGermline" )] )
  } else return( data.frame() )
}
#' Function to parse DNAStringSet into nested list
#'
#' @param DNAStrings object of class \code{DNAStringsSet}, holding the sequences. Can be obtained by \code{Biostrings::readDNAstringSet( fasta )}, where 'fasta' is the FASTA file holding the sequences.
#' 
#' @description This function parses user-supplied sequences to give a 'nested list' of sequences as desired by \code{cloneLineage} to
#' fetch germline sequence to perform lineage inference. NOTE: it is assumed sequences are directly fetched from IMGT and obeys the FASTA header format imposed by IMGT!
#'
#' @return A list of list, the first layer organised by species and the second layer organised by allele. 
#'
convertDNAStringSetToList <- function(DNAStrings)
{
  seqnames <- names( DNAStrings )
  tb <- do.call("rbind", lapply(seqnames, function( x ){
    seq <- as.character( DNAStrings[x] )
    x <- unlist( strsplit( x, split = "|", fixed = TRUE ) )
    species <- unlist( strsplit( x[3], split = "_" ) )[1]
    allele <- x[2]
    data.frame( species = species, allele = allele, sequence = seq,
                stringsAsFactors = FALSE)
  }))
  o <- lapply( unique( tb$species ), function( sp ){
    alleles <- tb[ which( tb$species == sp ), "allele" ]
    oo <- lapply( alleles, function( allele ) tb[ which ( tb$species == sp &
                                                            tb$allele == allele ), 
                                                  "sequence"])
    names(oo) <- alleles
    oo
  })
  names( o ) <- unique( tb$species )
  o
}
