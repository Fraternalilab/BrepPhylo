#' Function to analyse B cell clone lineages, by producing an alignment and infering a phylogenetic tree.
#'
#' @param input Matrix / data frame with values.
#' @param outputFolder Parent folder, in which the clones are to be deposited.
#' @param species Either "Bos_taurus", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Oryctolagus_cuniculus", "Sus_scrofa." (default: "Homo_sapiens").
#' @param minCloneSize Minimum clone size to be considered, clones with less members are not processed (default: 10).
#' @param whichClones If not \code{NULL} but a vector of clone IDs (must match the levels in column "cloneIdColumn"), only these clones are considered. Note, that parameter \code{minCloneSize} is still in effect.
#' @param ignoreClones This parameter can be set to one of the following:
#' \describe{
#'   \item{\code{NULL}}{No clones will be ignored and all that match the minimum size specified will be calculated. This is the default.}
#'   \item{"completed"}{The program will try to load the file "clones_completed.csv" in the directory specified by \code{outputFolder}. If found, it will ignore clones that are in the "cloneID" column of this file. If the file has not been generated yet or deleted, the program and will proceed with all clones. The file "clones_completed.csv" is generated and updated automatically by this function.}
#'   \item{A vector of clone IDs}{Any clone ID in this vector will not be considered.}
#' }
#' Note, that this parameter must be \code{NULL} if parameter \code{whichClones} is specified.
#' @param colourSpecification List with two elements:
#' \describe{
#'   \item{colours}{Named list, where the levels are the names and the HTML colour the value.}
#'   \item{coloursColumn}{Name of the column, which holds the levels for colouring.}
#'   \item{germlineColour}{The HTML colour value of the germline sequence.}
#' }
#' @param sequenceColumn Name of the column holding the sequences.
#' @param cloneIdColumn Name of the column holding the clone IDs.
#' @param labelColumn Name of the column holding the labels for the plots. If \code{NULL}, the index will be used.
#' @param germlineIdColumn Name of column holding the germline identifier (default: "V.GENE.and.allele").
#' @param germlineSet Optional, filepath pointing to a FASTA file containing germline sequences to be used.
#' @param algorithm Either "ClustalOmega" or "ClustalW" (pretty slow, though), specifying the alignment algorithm (default: "ClustalOmega").
#' @param plotType Either "phylogram", "cladogram", "fan" or "radial", specifying the layout of the phylogenetic tree (default: "phylogram").
#' @param plotFormat Either 'png' or 'pdf'. (default: 'png')
#' @param gatherPlots If \code{TRUE}, the generated plots will be copied from the individual clone folders to a folder "plots" in the output directory (default: \code{TRUE}).
#' @param useTempDir If \code{TRUE}, generate temporary directory and write results there. (default: \code{TRUE})
#' @param treeConstruction A list with the following elements:
#' \describe{
#'   \item{type}{"simple", "dnapars" or "igphyml". "simple" refers to a neighbour-joining tree; "dnapars" constructs a maximum parsimony tree using the phylip "dnapars" program (a local installation of the program is required). "igphyml" constructs a maximum likelihood tree taking into account mutational hotspot contexts in immunoglobulins, but is the most time consuming. Default "dnapars".}
#'   \item{parameters}{List of parameters. If \code{type} is "simple", always \code{NULL}. If \code{type} is "dnapars", an element \code{collapseAndReplace} can be set, which is either TRUE (default) or FALSE. This controls whether identical sequences are collapsed prior to tree construction, and subsequently be placed back into the tree as sister nodes of its duplicates. This is recommended to save computational time while preserving accuracy of tree construction. If \code{type} is "igphyml", an element \code{accuracy} can be set, which is either "basic", "high" or "extreme" (see manual of function \link{IgPhyML} for details).}
#' }
#' @param makeAlignmentPDFConstruction A list with the following elements:
#' \describe{
#'   \item{makeAlignmentPDF}{If \code{TRUE}, an alignment as PDF is printed out (default: \code{TRUE}).}
#' }
#' 
#' @examples
#' \dontrun{ cloneLineage <- function( input = input,
#'                           outputFolder = "~/Desktop/out",
#'                           species = "Homo_sapiens",
#'                           minCloneSize = 300,
#'                           colourSpecification = list( colours = list( "IgA" = "#0000FF",
#'                                                                       "IgG" = "#FF0000",
#'                                                                       "IgM" = "#00FF00" ),
#'                                                       coloursColumn = "Class",
#'                                                       germlineColour = "#AAAAAA" ),
#'                           sequenceColumn = "Seqs",
#'                           cloneIdColumn = "CloneID",
#'                           germlineIdColumn = "V.GENE.and.allele",
#'                           labelColumn = "SeqID",
#'                           algorithm = "ClustalOmega",
#'                           plotType = "fan",
#'                           gatherPlots = TRUE,
#'                           treeConstruction = list( "type" = "igphyml",
#'                                                    "parameters" = list( "accuracy" = "basic" ) ),
#'                           makeAlignmentPDFConstruction = list( "makeAlignmentPDF" = FALSE ) )}
#'
#' @importFrom utils read.csv write.csv capture.output
#' @importFrom methods new
#' @importFrom grDevices dev.off pdf png
#' @importFrom graphics legend
#' @importFrom Biostrings DNAStringSet readDNAMultipleAlignment readDNAStringSet writeXStringSet
#' @importFrom bio3d read.fasta write.fasta
#' @importFrom alakazam getPathLengths
#' @importFrom msa msa msaConvert
#' @importFrom alakazam buildPhylipLineage graphToPhylo
#' @importFrom igraph V neighbors edge.attributes get.edge.ids vertex edge
#' @importFrom phytools bind.tip
#' @importFrom ape drop.tip root nj as.phylo di2multi write.tree plot.phylo tiplabels
#' @importFrom seqinr dist.alignment read.alignment
#' 
#' @export
cloneLineage <- function( input,
                          outputFolder,
                          species = "Human",
                          minCloneSize = 10,
                          whichClones = NULL,
                          ignoreClones = NULL,
                          colourSpecification,
                          sequenceColumn,
                          cloneIdColumn,
                          germlineIdColumn = "V.GENE.and.allele",
                          labelColumn = NULL,
                          germlineSet = NULL,
                          algorithm = "ClustalOmega",
                          plotType = "phylogram",
                          plotFormat = "png",
                          gatherPlots = TRUE,
                          useTempDir = TRUE,
                          treeConstruction = list( "type" = "dnapars",
                                                   "parameters" = NULL ),
                          makeAlignmentPDFConstruction = list( "makeAlignmentPDF" = FALSE ) )
{
  # NOTE: The object "Vgermline" is part of the package and loaded automatically when the package is loaded.
  #       The data is stored in file "data/Vgermline.rda" and it might be necessary to update it at one point.
  #       In order to do that, one needs to generate a list called "Vgermline", with the "species" as named elements
  #       which in turn are named lists, where the IMGT V gene is the key and the nucleic acid sequence the value.
  #       For example: Vgermline <- list( "Human" = list( "IGHV1-1*1" = "acagatgactagctagtact", <and so on>
  #       The object can then be exported to overwrite the original one using: usethis::use_data( Vgermline, internal = FALSE, overwrite = TRUE )
  
  # check input
  if( any( !( c( colourSpecification[[ "coloursColumn" ]], sequenceColumn, cloneIdColumn, germlineIdColumn ) %in%
              colnames( input ) ) ) )
    stop( "One column or multiple columns specified missing." )
  if( !is.null( labelColumn ) && !( labelColumn %in% colnames( input ) ) )
    stop( "Specified labelColumn does not exist." )
  if( !( algorithm %in% c( "ClustalOmega", "ClustalW" ) ) )
    stop( "Parameter \"algorithm\" must be either \"ClustalOmega\" or \"ClustalW\"." )
  if( !( plotType %in% c( "fan", "phylogram", "cladogram", "radial" ) ) )
    stop( "Parameter \"plotType\" must be either \"fan\", \"phylogram\", \"cladogram\" or \"radial\"." )
  if( !file.exists( outputFolder ) )
    stop( "Specified outputFolder does not exist." )
  if( nrow( input ) < 1 )
    stop( "Parameter \"input\" must have at least one row." )
  if( !( treeConstruction[[ "type" ]] %in% c( "simple", "dnapars", "igphyml" ) ) )
    stop( "Parameter \"type\" of input list \"treeConstruction\" must be either \"simple\", \"dnapars\" or \"igphyml\"." )
  if( !is.logical( gatherPlots ) )
    stop( "Parameter \"gatherPlots\" must be a boolean value." )
  if( !is.null( whichClones ) && !is.null( ignoreClones ) )
    stop( "Parameters \"whichClones\" and \"ignoreClones\" cannot be specified simultaneously." )
  if( !is.null( germlineSet ) && !is.character( germlineSet ))
    stop( "germlineSet should be a character specifying a FASTA file containing germline sequences to be used." )
  if( !is.null( germlineSet ) && !file.exists( germlineSet ))
    stop( "germlineSet is not NULL but the supplied file does not exist." )
  if( !plotFormat %in% c("png", "pdf"))
    stop( "plotFormat must be either 'png' or 'pdf'." )
  if( !is.list( treeConstruction ) )
    stop( "'treeConstruction' should be a list. See documentation for details." )
  if( ! "type" %in% names( treeConstruction ))
    stop( "'type' must be one of the item names in the list 'treeConstruction'. See documentation for details." )
  if( ! treeConstruction$type %in% c( "simple", "dnapars", "igphyml" ) )
    stop( "'type' in 'treeConstruction' must be any one of: 'simple', 'dnapars' or 'igphyml'." )
  if( treeConstruction$type == "dnapars" ){
    if( "executable" %in% names( treeConstruction$parameters ) ){
      dnapars_path <- treeConstruction$parameters$executable
      if( !file.exists( dnapars_path )){
        stop( paste0( c( "Could not find dnapars executable. Have you indicated the correct path?\n",
                         "If you have installed the software system-wide you can omit indication of 'executable'",
                         " in treeConstruction parameters - the function will automatically look for it.\n") ) )
      }
    } else {
      dnapars_path <- as.character( Sys.which( names = "dnapars" ) )
      if( !file.exists( dnapars_path ) )
        stop( paste0( c( "Could not find dnapars executable. If installed, make sure you append the path to the",
                         " executable, e.g. in linux:\necho \"export PATH=$PATH:<path>/exe\" >> ~/.bash_profile\n",
                         "where <path> is the full path to the phylip folder with the\nsubfolder \"exe\" containing",
                         " the phylip executables.\n") ) )
    }
    if( "collapseAndReplace" %in% names( treeConstruction$parameters ) ){
      if( ! treeConstruction$parameters$collapseAndReplace %in% c(TRUE, FALSE) ){
        stop("\"collapseAndReplace\" parameter in \"treeConstruction\" must be either TRUE or FALSE.")
      } else collapseAndReplace <- treeConstruction$parameters$collapseAndReplace
    } else collapseAndReplace <- TRUE
  }
  if( treeConstruction$type == "igphyml" ){
    if( "executable" %in% names( treeConstruction$parameters ) ){
      igphyml_path <- treeConstruction$parameters$executable
      if( !file.exists( igphyml_path )){
        stop( paste0( c( "Could not find dnapars executable. Have you indicated the correct path?\n",
                         "If you have installed the software system-wide you can omit indication of 'executable'",
                         " in treeConstruction parameters - the function will automatically look for it.\n") ) )
      }
    } else {
      igphyml_path <- ""
    }
    if( "accuracy" %in% names( treeConstruction$parameters ) ){
      if( ! treeConstruction$parameters$accuracy %in% c("basic", "high", "extreme") ){
        stop("\"accuracy\" parameter in \"treeConstruction\" must be any one of 'basic', 'high' or 'extreme'.")
      } 
    } else treeConstruction$parameters$accuracy <- "basic"
  }
  
  # If germlineSet is supplied and the file exists, use this and initialise Vgermline
  # otherwise Vgermline will be taken from the lazy-load Vgermline set.
  if( !is.null( germlineSet ) && file.exists( germlineSet )){
    Vgermline <- Biostrings::readDNAStringSet( germlineSet )
    Vgermline <- convertDNAStringSetToList( Vgermline )
  }
  
  # check if Vgermline agrees with supplied sequence column in terms of IMGT gapping.
  # if one is gapped but the other is not, throw error.
  if( class( Vgermline ) == "DNAStringSet" ){
    germline_dot <- any( grepl( ".", Vgermline, fixed = TRUE ) )
  } else if( class( Vgermline ) == "list" ){
    germline_dot <- grepl( ".", Vgermline[[ 1 ]][[ 1 ]], fixed = TRUE )
  }
  sequence_dot <- any( grepl( ".", input[, sequenceColumn ], fixed = TRUE ) )
  if( sequence_dot && !germline_dot && is.null( germlineSet )){
    # load gapped version
    Vgermline <- Vgermline_gapped
    germline_dot <- grepl( ".", Vgermline[[ 1 ]][[ 1 ]], fixed = TRUE )
  }
  if( germline_dot != sequence_dot )
    stop("Only one of the germline set and the sequence column is gapped (containing '.'). If the sequences to be analysed is gapped please provide a gapped germline FASTA file using parameter 'germlineSet' (such sequences can be downloaded from IMGT). If sequences to be analysed are NOT gapped but you supplied a germlineSet with gapped sequence, please either supplied an ungapped germline set, or set germlineSet = NULL to use an in-built germline database.")

  if( !( species %in% names( Vgermline ) ) )
    stop( "Specified species is not supported." )
  
  # if a list of clones has been specified, set the minimum size of a clone to 1 and select the respective clones
  if( !is.null( whichClones ) )
    input <- input[ input[ , cloneIdColumn ] %in% whichClones, , drop = FALSE ]

  # if a list of clones has been specified to be ignored, deselect the respective clones
  if( !is.null( ignoreClones ) )
  {
    if( length( ignoreClones ) == 1 && ignoreClones == "completed" ) {
      # try to load the file "clones_completed.csv" in the output directory, which may have been made by a previous
      # run or engineered by hand; if the file is not present, proceed the usual way, otherwise remove all clones affected
      if( file.exists( paste0( outputFolder, "/clones_completed.csv" ) ) )
      {
        ignoreClones <- read.csv( file = paste0( outputFolder, "/clones_completed.csv" ) )[ , "cloneID" ]
      } else { ignoreClones <- NA }
    }

    # do not consider clone IDs in the "ignoreClones" vector
    input <- input[ !( input[ , cloneIdColumn ] %in% ignoreClones ), , drop = FALSE ]
  }

  # initialise variables
  set.seed( 3 )
  vecCloneIDs <- unique( input[ , cloneIdColumn ] )
  skippedClones <- 0

  # check, if there are any clones left
  if( length( vecCloneIDs ) < 1 )
    stop( "No clones found with this selection." )

  # prepare the colours (set all to capital letters for LaTeX afterwards)
  for( iColour in 1:length( colourSpecification[[ "colours" ]] ) )
    colourSpecification[[ "colours" ]][[ iColour ]] <- toupper( trimws( x = colourSpecification[[ "colours" ]][[ iColour ]], which = "both" ) )

  # write out colour code
  dfColours <- data.frame( names( colourSpecification[[ "colours" ]] ),
                           unlist( colourSpecification[[ "colours" ]] ) )
  names( dfColours ) <- c( colourSpecification[[ "coloursColumn" ]], "colours" )
  write.csv( x = dfColours,
             file = paste0( outputFolder, "/colour_code.csv" ),
             row.names = FALSE )
  rm( dfColours )

  # write the parameters out
  parametersLogPath <- paste0( outputFolder, "/parameters.log" )
  if( file.exists( parametersLogPath ) )
    unlink( x = parametersLogPath )
  parametersLog <- file( description = parametersLogPath, open = 'w' )
  asChar_labelColumn <- labelColumn
  if( is.null( asChar_labelColumn  ) )
    asChar_labelColumn <- "NULL"
  cat( paste0( "# ", curDateTime(), " Parameters of function call \"cloneLineage()\"",
               "\nspecies: \"", species, "\"",
               "\nminCloneSize: \"", minCloneSize, "\"",
               "\nsequenceColumn: \"", sequenceColumn, "\"",
               "\ncloneIdColumn: \"", cloneIdColumn, "\"",
               "\ngermlineIdColumn: \"", germlineIdColumn, "\"",
               "\nlabelColumn: \"", asChar_labelColumn, "\"",
               "\nalgorithm: \"", algorithm, "\"" ),
       file = parametersLog )
  close( parametersLog )

  # prepare an plot output folder, if required
  if( gatherPlots && !file.exists( paste0( outputFolder, "/plots" ) ) )
    dir.create( paste0( outputFolder, "/plots" ), showWarnings = FALSE )

  # do the clone analysis
  for( iClone in 1:length( vecCloneIDs ) )
  {
    result <- tryCatch({

      # get the current clone and skip it, if the size threshold is not reached
      cloneMembers <- which( input[ cloneIdColumn ] == vecCloneIDs[ iClone ] )
      cloneMembers <- cloneMembers[ !is.na( cloneMembers ) ]
      cloneSize <- length( cloneMembers )
      if( cloneSize < minCloneSize )
      {
        skippedClones <- skippedClones + 1
        next
      }
      curClone <- input[ cloneMembers, , drop = FALSE ]

      # make a directory to store the results for this specific clone
      curDir <- paste0( outputFolder, "/clone_", vecCloneIDs[ iClone ] )
      dir.create( curDir, showWarnings = FALSE )

      # start the log file
      logFilePath <- paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".log" )
      if( file.exists( logFilePath ) )
        unlink( x = logFilePath )
      curClone_log <- file( description = logFilePath, open = 'a' )
      cat( paste0( "# ", curDateTime(), " Calculation commenced" ), file = curClone_log )

      # store the current clone
      write.csv( x = curClone,
                 file = paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".csv" ),
                 row.names = FALSE )
      cat( paste0( "\n# ", curDateTime(), " Wrote clone ", vecCloneIDs[ iClone ], " out (dimensions: ", nrow( curClone ), 'x', ncol( curClone ), ')' ), file = curClone_log )

      # get the germline sequence for a representative observation
      # JN: change from just the first one to the germline Vgene most frequently found in the set
      representativeObs <- curClone[ , germlineIdColumn ]
      germlineToUse <- sapply(representativeObs, function(x){
        unlist( strsplit( as.character(x), split = ' ', fixed = TRUE) )[1]
      })
      germlineToUse <- names( which.max( table( germlineToUse ) ) )[1]
      germline_sequence <- Vgermline[[ species ]][[ germlineToUse ]]
      if( is.null( germline_sequence ) )
        stop( paste0( "Specified germline sequence \"", germlineToUse, "\" has not been found in internal database." ) )
      cat( paste0( "\n# ", curDateTime(), " Use ", species, " germline \"", germlineToUse, "\": ", germline_sequence ), file = curClone_log )

      # get the names as a vector; if there is no specified "labelColumn", use just the numbers
      # in the current clone
      vecNames <- c( germlineToUse, 1:nrow( curClone ) )
      if( !is.null( labelColumn ) )
      {
        vecNames <- c( germlineToUse, trimws( as.character( curClone[ , labelColumn ] ), which = "both" ) )
        if( length( vecNames ) != length( unique( vecNames ) ) )
        {
          cat( paste0( "\n# ", curDateTime(), " Skipped clone ", vecCloneIDs[ iClone ], " since the specified parameter \"labelColumn\" does not have unique identifiers." ), file = curClone_log )
          close( curClone_log )
          next
        }
      }

      # do the alignment
      # JN: changes due to switching to use IMGT gapped version. Trim germline sequence to closest multiple of 3 (e.g. 320 --> 318)
      # and then trim the gapped repertoire to the same length.
      # Next write this into a temporary FASTA and read using bio3d::read.fasta
      # and then proceed with preparing/removing sequences to remove those with excessive STOPs (see below).
      if( sequence_dot & germline_dot ){
        if( treeConstruction$type == "igphyml" ){
          round_to <- floor( nchar( germline_sequence ) / 3 ) * 3
          germline_sequence <- substr( germline_sequence, 1, round_to )
          sequences <- sapply( curClone[ , sequenceColumn ], function( x ){
            substr( x, 1, round_to )
          })
        } else {
          sequences <- curClone[ , sequenceColumn ]
        }
        sequences <- c( germline_sequence, sequences )
        DNAss <- Biostrings::DNAStringSet( trimws( sequences ) )
        names( DNAss ) <- vecNames
        alignment <- tempfile()
        Biostrings::writeXStringSet( DNAss, alignment )
        # transform the alignment into a bio3d object, prepare the FASTA path and write it out
        alignment_FASTA <- bio3d::read.fasta( alignment )
      } else {
        # if not using gapped sequences then do conventional MSA
        sequences <- c( germline_sequence, as.character( curClone[ , sequenceColumn ] ) )
        DNAss <- Biostrings::DNAStringSet( trimws( sequences ) )
        names( DNAss ) <- vecNames
        outputMSA <- capture.output( alignment <- msa::msa( inputSeqs = DNAss, method = algorithm, order = "input", type = "dna" ) )
        # transform the alignment into a bio3d object, prepare the FASTA path and write it out
        alignment_FASTA <- msa::msaConvert( x = alignment, type = "bio3d::fasta" )
      }
      cat( paste0( "\n# ", curDateTime(), " Alignment calculated" ), file = curClone_log )
      
      # BUT: As IgPhyML needs a proper alignment with a sequence length
      #      equal to a multiple of 3, i.e. full codons, one needs to write a "pruned" one
      #      if the tree is to be constructed using IgPhyML; apart from this, STOP codons must be removed
      if( treeConstruction[[ "type" ]] == "igphyml" )
      {
        # JN: sometimes the reading frame is wrong in that gaps are introduced
        # before the beginning of the germline (which should be in frame!)
        # trim those gaps from the alignment before pruning alignment
        germline_aligned <- alignment_FASTA[[ "ali" ]][ germlineToUse , ]
        germline_aligned <- paste( germline_aligned, collapse = "" )
        germline_gaps <- gregexpr( "-" , germline_aligned, fixed = TRUE )[[ 1 ]]
        if( !all( germline_gaps == -1 ) ){
          # find the beginning of the germline sequence (ie after removing gaps)
          trim_from <- which( sapply( 1:nchar(germline_aligned), 
                                      function(x) !(x %in% germline_gaps) ) )[ 1 ]
          alignment_FASTA[[ "ali" ]] <- alignment_FASTA[[ "ali" ]][, trim_from:ncol(alignment_FASTA[[ "ali" ]]) ]
        }
        valCodReturn <- alignmentValidCodons( inputAlignment = alignment_FASTA,
                                              cutEnd = TRUE,
                                              removeSTOPcodons = TRUE,
                                              do.warn = FALSE )
        alignment_FASTA <- valCodReturn[[ "alignment" ]]
        cat( paste0( "\n# ", curDateTime(), " Checked alignment for validity, i.e. if IgPhyML requirements are met: ", valCodReturn[[ "report" ]][[ "noCutPositions" ]], " positions were cut (end)" ), file = curClone_log )
        listSTOPreplaced <- valCodReturn[[ "report" ]][[ "STOPreplaced" ]]
        listSTOPremoved <- valCodReturn[[ "report" ]][[ "STOPremoved" ]]
        if( length( listSTOPreplaced ) > 0 )
          for( iSTOPreplaced in 1:length( listSTOPreplaced ) )
            cat( paste0( "\n# --- Replaced ", listSTOPreplaced[[ iSTOPreplaced ]], " first position(s) in STOP codon(s) in sequence \"", names( listSTOPreplaced )[ iSTOPreplaced ], "\"" ), file = curClone_log )
        if( length( listSTOPremoved ) > 0 )
          for( iSTOPremoved in 1:length( listSTOPremoved ) )
            cat( paste0( "\n# --- Removed sequence \"", names( listSTOPremoved )[ iSTOPremoved ], "\" as it contained ", listSTOPremoved[[ iSTOPremoved ]], " STOP codon(s)" ), file = curClone_log )
        if( length( alignment_FASTA[[ "id" ]] ) == 0 | germlineToUse %in% names( listSTOPremoved ) )
        {
          cat( paste0( "\n# ", curDateTime(), " Skipped clone ", vecCloneIDs[ iClone ], " as no sequences remain after codon validation (try to set \"treeConstruction$type\" to \"simple\")." ), file = curClone_log )
          close( curClone_log )
          next
        }
      }
      pathFASTA <- paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".fasta" )
      bio3d::write.fasta( alignment = alignment_FASTA, file = pathFASTA )
      if( nrow(alignment_FASTA$ali) < minCloneSize )
      {
        skippedClones <- skippedClones + 1
        next
      }

      # depending on the way the tree is to be constructed, either do it yourself or use the IgPhyML package
      # prepare the tree object
      phyloTree <- NULL
      if( treeConstruction[[ "type" ]] == "igphyml" )
      {
        # call function "IgPhyML()", which will construct a phylogenetic tree for the alignment
        IgPhyML <- IgPhyML( path = pathFASTA,
                            root = germlineToUse,
                            accuracy = treeConstruction[[ "parameters" ]][[ "accuracy" ]],
                            log = curClone_log, 
                            executable = igphyml_path,
                            useTempDir = useTempDir)

        if( !IgPhyML[[ "success" ]] )
        {
          cat( paste0( "\n# ", curDateTime(), " Problem occured, when invoking IgPhyML - skip clone:\n\n\n", IgPhyML[[ "error_message" ]] ), file = curClone_log )
          warning( paste0( "Problem occured, when invoking IgPhyML for clone ", vecCloneIDs[ iClone ], ", which will be skipped." ) )
          close( curClone_log )
          next
        }
        cat( paste0( "\n# ", curDateTime(), " Constructed tree (IgPhyML, accuracy: \"", treeConstruction[[ "parameters" ]][[ "accuracy" ]], "\")" ), file = curClone_log )

        # set the phylo tree
        phyloTree <- IgPhyML[[ "fit" ]][[ "tree" ]]

        # write the stats out
        writeLines( text = IgPhyML[[ "fit" ]][[ "stats" ]], con = paste0( curDir, "/clone_", vecCloneIDs[ iClone ], "_IgPhyML_fit.stats" ) )
      }
      else if( treeConstruction[[ "type" ]] == "dnapars" ){
        tb <- Biostrings::readDNAStringSet( pathFASTA )
        germline_seq <- as.character( tb[ germlineToUse ] )
        tb <- as.data.frame( tb )
        tb$collapsed <- 1
        tb$label <- rownames( tb )
        tb <- tb[, c( 3,1,2 )]
        colnames(tb) <- c("sequence_id", "sequence", "collapse_count")
        if( collapseAndReplace ){
          # collapse repeated sequences
          which_dup <- which( duplicated( tb$sequence ) )
          if( length( which_dup ) > 0 ){
            tb_unique <- tb[ -which_dup, ]
          } else tb_unique <- tb
        } else tb_unique <- tb
        # if collapsed only has 1 sequence left and it is the germline,
        # output message and skip this clone
        if( nrow( tb_unique ) == 1 & tb_unique[ 1, 1] == germlineToUse ){
          cat( paste0( "\n# ", curDateTime(), " Constructed tree (dnapars) skipped as this clone contains only copies of the germline sequence." ), file = curClone_log )
          next
        }
        changeoclone <- new( "ChangeoClone", data = tb_unique,
                   clone = as.character( vecCloneIDs[ iClone ] ), 
                   germline = gsub(".", "-", germline_seq, fixed = TRUE),
                   v_gene = germlineToUse, j_gene = "", junc_len = 1 ) # these are useless annotation just to be added later
        if( useTempDir ){
          tmppath <- NULL
          rmtmp <- TRUE
        } else {
          tmppath <- curDir
          rmtmp <- FALSE
        }
        phyloGraph <- alakazam::buildPhylipLineage( changeoclone, dnapars_path, verbose = FALSE,
                                                    rm_temp = rmtmp, 
                                                    temp_path = tmppath,
                                                    branch_length = "distance")
        if( collapseAndReplace ){
          # for each removed sequence in the collapsing, add them back to the tree as
          # sister leaves to their duplicates.
          tips_to_add <- tb[ which( !tb$sequence_id %in% names( igraph::V(phyloGraph) ) ), 
                             "sequence_id" ]
          if( length( tips_to_add ) > 0 ){
            for( dropped in tips_to_add) {
              identical_tip <- 
                tb_unique[ which( tb_unique$sequence == tb[ which(tb$sequence_id == dropped), "sequence"] ), "sequence_id" ]
              ancestor <- names( igraph::neighbors( phyloGraph, identical_tip, mode = "in") )
              edgeWeight <- igraph::edge.attributes( phyloGraph,
                                                     igraph::get.edge.ids( phyloGraph,
                                                                           vp = c(ancestor, identical_tip) ) )$weight
              phyloGraph <- phyloGraph + igraph::vertex( dropped,
                                                         sequence = tb[ which(tb$sequence_id == dropped), "sequence"],
                                                         clone = as.character( vecCloneIDs[ iClone ] ), 
                                                         v_gene = germlineToUse, j_gene = "", junc_len = 1)
              phyloGraph <- phyloGraph + igraph::edge( ancestor, dropped,
                                                       weight = edgeWeight,
                                                       label = edgeWeight)
              #            phyloTree <- phytools::bind.tip( phyloTree, 
              #                                             tip.label = dropped,
              #                                             edge.length = 0,
              #                                             where = which( phyloTree$tip.label == identical_tip ))
            }
          }
        }
        phyloTree <- alakazam::graphToPhylo( phyloGraph )
        phyloTree <- ape::drop.tip( phyloTree, "Germline" )
        phyloTree <- ape::root( phyloTree, germlineToUse )
        cat( paste0( "\n# ", curDateTime(), " Constructed tree (dnapars)" ), file = curClone_log )
      }
      else
      {
        # get the distances
        distances <- as.matrix( seqinr::dist.alignment( seqinr::read.alignment( pathFASTA, format = "fasta" ) ) )
        colnames( distances ) <- vecNames
        rownames( distances ) <- vecNames

        # build the lineage tree
        tree <- ape::nj( distances )
        phyloTree <- ape::as.phylo( x = tree )
        phyloTree <- ape::root( tree, germlineToUse )
        # collapse multichotomies (i.e. branches with very short lengths)
        phyloTree <- ape::di2multi(phyloTree, top = 0.01)
        cat( paste0( "\n# ", curDateTime(), " Constructed tree (simple)" ), file = curClone_log )
      }

      # write out the tree in Newick format
      ape::write.tree( phy = phyloTree,
                       file = paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".tree" ) )
      cat( paste0( "\n# ", curDateTime(), " Tree written" ), file = curClone_log )

      # prepare the colours for the phylogeny tree; make sure, the colours match the
      # tips and labels; note, that "match()" will automatically take care of observations
      # lost during the codon validation
      listColours <- colourSpecification[[ "colours" ]][ as.character( curClone[ , colourSpecification[[ "coloursColumn" ]] ] ) ]
      isUNK <- unlist( lapply( X = listColours, FUN = function( x ) is.null( x ) ) )
      listColours[ isUNK ] <- "#000000"
      vecColours <- c( colourSpecification[[ "germlineColour" ]],
                       unlist( listColours ) )
      names( vecColours )[ is.na( names( vecColours ) ) ] <- "UNK"
      names( vecColours )[ 1 ] <- "germline"
      vecColoursTree <- vecColours[ base::match( phyloTree$tip.label, vecNames ) ]

      # make the plot
      if( plotFormat == "png"){
        plotPath <- paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".png" )
        png( width = 1800, height = 1800, filename = plotPath )
      } else if( plotFormat == "pdf" ){
        plotPath <- paste0( curDir, "/clone_", vecCloneIDs[ iClone ], ".pdf" )
        pdf( width = 18, height = 18, file = plotPath )
      }
      ape::plot.phylo( x = phyloTree,
                       type = plotType,
                       tip.color = vecColoursTree,
                       label.offset = 0.035 )
      ape::tiplabels( pch = 19, col = vecColoursTree )
      legend( "topright",
              legend = c( "germline", unique( names( colourSpecification$colours ) ) ),
              col = c( colourSpecification$germlineColour, unlist( unique( colourSpecification$colours ) ) ),
              pch = 19,
              cex = 2.5,
              pt.cex = 2.25,
              bty = 'n' )
      dev.off()
      if( gatherPlots )
        base::file.copy( from = plotPath, to = paste0( outputFolder, "/plots/", basename( plotPath ) ) )
      cat( paste0( "\n# ", curDateTime(), " Tree plotted" ), file = curClone_log )

      # make the alignment as PDF, if requested
      if( makeAlignmentPDFConstruction[[ "makeAlignmentPDF" ]] )
      {
        alignment <- Biostrings::readDNAMultipleAlignment( pathFASTA )
        makeAlignmentPDF( clone = curClone,
                          alignment = alignment,
                          cloneIdColumn = cloneIdColumn,
                          outputFolder = curDir,
                          colours = vecColours )
        cat( paste0( "\n# ", curDateTime(), " Alignment PDF generated" ), file = curClone_log )
      }

      # update the list of completed clones
      addCloneToCompleted( folder = outputFolder, cloneID = vecCloneIDs[ iClone ], size = nrow( curClone ) )

      # show progress, store end time and close log file
      cat( paste0( "Clone ", vecCloneIDs[ iClone ], " (", cloneSize, " members) has been completed." ), sep = "\n" )
      cat( paste0( "\n# Calculation completed: ", curDateTime() ), sep = "\n", file = curClone_log )
      close( con = curClone_log )

    }, error = function( e ) {
      cat( paste0( "Calculation of clone ", vecCloneIDs[ iClone ], " failed\nError message:\n", toString( x = e ), "\nProceeding with next clone!\n\n\n" ) )
      return( NULL )
    })
  }
  cat( paste0( "A total of ", skippedClones, " clones (of ", length( vecCloneIDs ), ") have been skipped." ), sep = "\n" )
}
