#' Perform phylogenetic inference & tree analysis in a batch
#'
#' @param inputDF input data.frame containing repertoire data.
#' @param outputFolder character, path for the folder where output files are to be stored.
#' @param species characer, species. For now only "Human" is supported and accepted.
#' @param minCloneSize numeric, minimum clone size to be considered for tree inference. (default: 3)
#' @param sequence_column character, column name in \code{inputDF} which holds the sequences. (default: "Sequence")
#' @param cloneID_column character, column name in \code{inputDF} which holds the clone IDs. (default: "CloneID")
#' @param label_column character, column name in \code{inputDF} which holds the sequence identifiers. (default: "Seq_ID")
#' @param colourColumn character, column name in \code{inputDF} which is to be categorised by colours in plots. (default: "Subclass")
#' @param IGHVgeneandallele_column character, column name in \code{inputDF} which holds annotated germline V gene names. (default: "V.GENE.and.allele")
#' @param germlineSet Optional, filepath pointing to a FASTA file containing germline sequences to be used. (default: NULL)
#' @param plotFormat Either 'png' or 'pdf' (default: 'png')
#' @param phyloTreeType "simple", "dnapars" or "igphyml". "simple" refers to a neighbour-joining tree; "dnapars" constructs a maximum parsimony tree using the phylip "dnapars" program (a local installation of the program is required). "igphyml" constructs a maximum likelihood tree taking into account mutational hotspot contexts in immunoglobulins, but is the most time consuming. (default: "dnapars")
#' @param phyloTreeOptions A list to be fed as the 'parameter' entry in the \code{treeConstruction} argument of the \code{cloneLineage} function (see documentation of the \code{cloneLineage} function). If \code{NULL}, default settings will be used. (default: \code{NULL})
#' @param makeArboTree Should arborescence tree be calculated? (default: TRUE)
#' @param useTempDir If \code{TRUE}, generate temporary directory and write results there. (default: \code{TRUE})
#'
#' @description This function loops through clones and perform IgPhyML phylogenetic tree inference and arborescence tree construction, and calculates tree metrics, as detailed in vignette of this package. 
#'
#' @return A list with each element named by a clone ID, each itself a list with two elements:
#' \describe{
#'   \item{distances}{data.frame with the distance-from-germline calculated from the lineage trees for each sequence in the clone. See Vignette for details.}
#'   \item{csr_events}{data.frame listing all class-switching events and an estimate distance-from-germline at which such event takes place.}
#' }
#' 
#' @importFrom utils read.csv
#' @importFrom ape read.tree
#' @export doBatchCloneAnalysis
doBatchCloneAnalysis <- function( inputDF,
                                  outputFolder,
                                  species = "Human",
                                  minCloneSize = 3,
                                  sequence_column = "Sequence",
                                  cloneID_column = "CloneID",
                                  label_column = "Seq_ID",
                                  colourColumn = "Subclass",
                                  IGHVgeneandallele_column = "V.GENE.and.allele",
                                  germlineSet = NULL,
                                  plotFormat = "png",
                                  phyloTreeType = "simple",
                                  phyloTreeOptions = NULL,
                                  makeArboTree = TRUE,
                                  useTempDir = TRUE)
{
  # initialise variables
  set.seed( 3 )
  vecCloneIDs <- unique( inputDF[ , cloneID_column ] )
  skippedClones <- 0
  
  # generate an output folder
  outputFolder <- path.expand( outputFolder )
  dir.create( outputFolder, showWarnings = FALSE )
  arboFolder <- paste0( outputFolder, "/arborescence" )
  dir.create( arboFolder, showWarnings = FALSE )
  
  output <- list()
  
  # do the clone analysis
  for( iClone in 1:length( vecCloneIDs ) )
  {
    # get the current clone
    cloneMembers <- which( inputDF[ cloneID_column ] == vecCloneIDs[ iClone ] )
    cloneMembers <- cloneMembers[ !is.na( cloneMembers ) ]
    cloneSize <- length( cloneMembers )
    if( cloneSize < minCloneSize )
    {
      skippedClones <- skippedClones + 1
      next
    }
    curClone <- inputDF[ cloneMembers, , drop = FALSE ]

    # generate an output folder
    outputf <- path.expand( paste0( outputFolder ) )
    dir.create( outputf, showWarnings = FALSE )
    # if analysis already run then skip it and don't waste time re-generate the outputs
    if( file.exists( paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", 
                             vecCloneIDs[ iClone ], ".tree" ) ) ){
      next
    }
    output[[ as.character( vecCloneIDs[ iClone ] ) ]] <- list()
    
    # run the tree reconstruction by calling IgPhyML
    if( !is.null( phyloTreeOptions ) ){
      params = phyloTreeOptions
    }
    if( phyloTreeType == "igphyml" ){
      if( is.null( phyloTreeOptions )) params <- list( "accuracy" = "basic")
      tree_type <- list( "type" = "igphyml",
                         "parameters" = params)
    } else if( phyloTreeType == "dnapars" ){
      if( is.null( phyloTreeOptions )) params <- list( "collapseAndReplace" = TRUE )
      tree_type <- list( "type" = "dnapars", 
                         "parameters" =  params)
    } else {
      tree_type <- list( "type" = "simple" )
    }
    cloneLineage( input = curClone,
                  outputFolder = outputf,
                  species = species,
                  whichClones = vecCloneIDs[ iClone ],
                  colourSpecification = list( colours = list( "IgM" = "#FF0000",
                                                              "IgD" = "#FF00FF",
                                                              "IgA" = "#9999CF",
                                                              "IgA1" = "#9999CF",
                                                              "IgA2" = "#000088",
                                                              "IgG" = "#326600",
                                                              "IgG1" = "#326600",
                                                              "IgG2" = "#58B200",
                                                              "IgG3" = "#8BFF19",
                                                              "IgG4" = "#CBFF99" ),
                                              coloursColumn = colourColumn,
                                              germlineColour = "#AAAAAA" ),
                  sequenceColumn = sequence_column,
                  cloneIdColumn = cloneID_column,
                  germlineIdColumn = IGHVgeneandallele_column,
                  minCloneSize = minCloneSize,
                  labelColumn = label_column,
                  germlineSet = germlineSet,
                  useTempDir = useTempDir,
                  plotFormat = plotFormat,
                  treeConstruction = tree_type,
                  makeAlignmentPDFConstruction = list( "makeAlignmentPDF" = FALSE ) )
    
    if( file.exists( paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], ".csv" ) ) && 
        file.exists( paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], ".tree" )) ){
      # show progress
      cat( paste0( "Clone ", vecCloneIDs[ iClone ], " lineage tree has been constructed.\n" ) )
      
      # make arborescence tree
      if( makeArboTree ){
        clone <- read.csv( paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], ".csv" ) )
        if( all( is.na( clone[, colourColumn] ) ) ){
          # no sequences have annotated subclass. skip arborescence tree
          next
        }
        tree <- paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], ".tree" )
        arborescence <- cloneArborescence( clone = clone,
                                           cloneID = vecCloneIDs[ iClone ],
                                           columnSeqID = label_column,
                                           columnSubclass = colourColumn, 
                                           tree = tree,
                                           plotFormat = "pdf",
                                           colourSpecification = list( colours = list( "IgM" = "#FF0000",
                                                                                       "IgD" = "#FF00FF",
                                                                                       "IgA" = "#9999CF",
                                                                                       "IgA1" = "#9999CF",
                                                                                       "IgA2" = "#000088",
                                                                                       "IgG" = "#326600",
                                                                                       "IgG1" = "#326600",
                                                                                       "IgG2" = "#58B200",
                                                                                       "IgG3" = "#8BFF19",
                                                                                       "IgG4" = "#CBFF99" ),
                                                                       germlineColour = "#AAAAAA" ),
                                           outputFolder = arboFolder )
        
        # show progress
        cat( paste0( "Clone ", vecCloneIDs[ iClone ], " arborescence tree has been considered.\n" ) )
        # germline <- arborescence$table[ grepl( "^IG.V", arborescence$table$startLabel ), "startLabel" ]
      }# else {
        #tree <- paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], ".tree" )
        #germline <- phytools::read.newick( tree )$tip.label
        #germline <- germline[ grepl( "^IG.V", germline ) ]
      #}
      # get tree stats
      if( phyloTreeType == "igphyml" ){
        fit_stats <- paste0( outputf, "/clone_", vecCloneIDs[ iClone ], "/clone_", vecCloneIDs[ iClone ], "_IgPhyML_fit.stats" )
      } else {
        fit_stats <- NULL
      }
      if( makeArboTree ){
        arbo_tree <- arborescence$graph
        csr_events <- arborescence$csr_events
      } else {
        arbo_tree <- NULL
        csr_events <- NULL
      }
      # tree <- ape::read.tree( tree )
      distances <- getGermlineDistance( tree )
      output[[ as.character( vecCloneIDs[ iClone ] ) ]] <- list( distances = distances,
                                                                 csr_events = csr_events )

      # show progress
      cat( paste0( "Clone ", vecCloneIDs[ iClone ], " tree statistics have been parsed.\n" ) )
    } else {
      cat( paste0( "Clone ", vecCloneIDs[ iClone ], " tree construction skipped.\n" ) )
    }
    
  }
    
  cat( paste0( "A total of ", skippedClones, " clones (of ", length( vecCloneIDs ), ") have been skipped due to the size requirement.\n" ) )
  return( output )
}