#' Function that calls IgPhyML to calculate a phylogenetic tree and returns it including meta-data.
#'
#' @param path Input FASTA file path.
#' @param root Root ID (germline sequence).
#' @param accuracy Either "basic", "high" or "extreme", this parameter affects the internal calls of IgPhyML which will affect runtime considerably (default: "basic").
#' @param executable If not empty, this file will be used for the call of \code{IgPhyML} otherwise the function will attempt to find the executable (default: "").
#' @param log If not \code{NULL}, log entries over the execution will be written to the given connection.
#' @param useTempDir If \code{TRUE}, generate temporary directory and write results there. (default: \code{TRUE})
#'
#' @description Note, that what this function implements is a two-step optimization: First, the topology is searched for using M0/GY49 (\href{https://doi.org/10.1093/oxfordjournals.molbev.a040153}{Goldman and Yang (1994)}; \href{https://doi.org/10.1016/S0169-5347(00)01994-7}{Yang and Bielawski (2000)}), then this topology is fixed and used 
#' as starting point for the fitting by the HLP17 model (specialized for antibody lineages \href{https://doi.org/10.1534/genetics.116.196303}{Hoehn, Lunter and Pybus (2017)}). One could further optimize the topology by specifying "-o tlr" instead of "-o lr".
#'
#' @return A list holding the intermediary and final trees and additional data (elements: "success", "topology", "fit" )/
#'
#' @examples
#' \dontrun{IgPhyML( path = "/path/to/input.fasta", root = "V4-59" )}
#' 
#' @importFrom phytools read.newick
#' 
#' @export
IgPhyML <- function( path,
                     root,
                     accuracy = "basic",
                     executable = "",
                     log = NULL,
                     useTempDir = TRUE)
{
  # list of IgPhyML options:
  #########
  # -i (or --input) seq_file_name (Required) | seq_file_name is the name of the nucleotide or amino-acid sequence file in PHYLIP format.
  # -q (or --sequential) | Changes interleaved format (default) to sequential format.
  # -m (or --model) model | model: substitution model name. Codon based models: HLP17 | GY (default)
  # -f (or --frequencies) empirical, model, optimized, fT, fC, fA, fG, fT1, fC1, fA1, fG1, fT2, fC2, fA2, fG2, fT3, fC3, fA3, fG3 or fC1, fC2, ... , fC64
  #                       empirical: (default) the character frequencies are determined as follows: 
  #                                  - Codon sequences: (Empirical) the equilibrium codon frequencies are estimated by counting
  #                                                     the occurence of bases or codons in the alignment according to the frequency model that is selected.
  # optimize: the character frequencies are determined as follows: 
  #            - Nucleotide sequences: (ML) the equilibrium base frequencies are estimated using maximum likelihood 
  #            - Amino-acid sequences: (ML) the equilibrium amino-acid frequencies are estimated using maximum likelihood
  #            - Codon sequences     : (ML) the equilibrium codon frequencies are estimated using maximum likelihood
  # --fmodel frequency model | Which frequency model to use.
  #                            frequency model = F1XCODONS | F1X4 | F3X4 | CF3X4 (default)
  # -t (or --ts/tv) ts/tv_ratio | This option sets the transition/transversion ratio.
  #                               ts/tv_ratio > 0: Set transition/transversion ratio to the value.
  #                               ts/tv_ratio = e: Get maximum likelihood estimate of transition/transversion ratio (default under HLP17).
  #                               ts/tv_ratio = 1: Fix estimate (default under GY94).
  # -w (or --omega) model (Required if a codon model other than HLP17 is used) The omega parameter or nonsynonymous/synonymous rate ratio.
  #                                                                            model = DM0: (default) Single omega model.
  # -s (or --search) move | Tree topology search operation option. Can be either NNI (default, fast) or SPR (a bit slower than NNI) or BEST (best of NNI and SPR search).
  # -u (or --inputtree) user_tree_file
  #                     user_tree_file: starting tree filename. The tree must be in Newick format.
  # -o (or --optimize) params | This option focuses on specific parameter optimisation.
  #                             params = tlr: (default) tree topology (t), branch length (l) and rate parameters (r) are optimised.
  #                             params = tl: tree topology and branch length are optimised.
  #                             params = lr: branch length and rate parameters are optimised.
  #                             params = l: branch length are optimised.
  #                             params = r: rate parameters are optimised.
  #                             params = n: no parameter is optimised.
  # --run_id ID_string | Append the string ID_string at the end of each PhyML output file. This option may be useful when running simulations involving PhyML.
  # --root rootID | Sets the root.
  # --threads Maximum number of threads.
  #########

  # check input
  if( !file.exists( path ) )
    stop( paste0( "Input file \"", path, "\" does not exist." ) )
  if( !( accuracy %in% c( "basic", "high", "extreme" ) ) )
    stop( "Parameter \"accuracy\" must be either \"basic\", \"high\" or \"extreme\"." )

  # check, whether "IgPhyML" is installed and available
  IgPhyML_path <- executable
  if( IgPhyML_path == "" )
    IgPhyML_path <- as.character( Sys.which( names = "igphyml" ) )
  if( !file.exists( IgPhyML_path ) )
    stop( paste0( c( "Could not find IgPhyML executable. If installed, set the path manually ",
                     " using parameter \"executable\". Installation instructions:\n\n 1) Install libatlas (dev version):\n",
                     "apt-get install libatlas-base-dev\n2) Download IgPhyML:\nhttps://github.com/kbhoehn/IgPhyML\n",
                     "3) Execute in IgPhyML folder:\n./make_phyml_blas_omp (or the adequate version without BLAS and / or OpenMP)\n",
                     "sudo make install\n" ) ) )

  # generate a temporary directory
  if( useTempDir ){
    tempDir <- paste0( tempdir(), "/IgPhyML" )
    if( dir.exists( paths = tempDir ) )
      unlink( x = tempDir, recursive = TRUE )
    dir.create( path = tempDir )
    setwd( dir = tempDir )
    if( !is.null( log ) )
      cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Switched to temporary directory" ), file = log )
    
    # copy input file to directory
    file.copy( from = path, to = paste0( tempDir, '/', basename( path ) ), overwrite = FALSE )
    path <- paste0( tempDir, '/', basename( path ) )
  } else {
    tempDir <- dirname( path )
    curDir <- getwd()
    setwd( tempDir )
    path <- basename( path )
    tempDir <- '.'
  }
  
  # prepare the arguments
  argumentsTopology <- c( paste0( "-i ", path ),
                          "-m GY",
                          "-w M0",    # as in manual
                          "-t e",     # as in manual
                          "--threads 1",
                          "--run_id gy94" )
  argumentsFit <- c( paste0( "-i ", path ),
                     "-o lr",
                     "-m HLP17",
                     paste0( "--root ", root ),
                     "--threads 1",
                     paste0( "-u ", basename( path ), "_igphyml_tree_gy94.txt" ) )
  if( accuracy %in% c( "high", "extreme" ) ) 
  {
    # replace the fast "Nearest Neighbour Interchange" (NNI) by the slower, but more sensitive "Subtree Pruning and Regrafting" (SPR)
    argumentsTopology <- c( argumentsTopology, "-s SPR" )
    argumentsFit <- c( argumentsFit, "-s SPR" )

    if( accuracy == "extreme" )
      # apart from the branch length and rate parameters, also optimize the tree topology further using the slow HLP17 model
      argumentsFit[ 2 ] <- "-o tlr"
  }

  # prepare return list
  result <- list( "success" = TRUE, "error_message" = NA, "topology" = NA, "fit" = NA )

  # execute the program
  # this will require a two-step fashion, where the tree is first approximated and then
  # refined; the GY model is fast, but does not correct for mutation hotspots
  # technical note: make input "\n" for finishing execution in case something went wrong
  if( !is.null( log ) )
    cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Commencing GY tree construction" ), file = log )
  
  base::system2( command = IgPhyML_path,
                 args = argumentsTopology,
                 stdout = paste0( tempDir, "/std_topology.out" ),
                 stderr = paste0( tempDir, "/std_topology.err" ),
                 input = "\n" )
  if( any( grepl( pattern = "Type enter to exit.", x = readLines( con = paste0( tempDir, "/std_topology.out" ), warn = FALSE ) ) ) )
  {
    if( !is.null( log ) )
      cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Error occured during execution, check the log files." ), file = log )
    result[[ "success" ]] <- FALSE
    result[[ "error_message" ]] <- paste( c( base::readLines( con = paste0( tempDir, "/std_topology.out" ), warn = FALSE ) ), collapse = "\n" )
  }
  else
  {
    if( !is.null( log ) )
      cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Commencing HLP17 tree refinement" ), file = log )
    base::system2( command = IgPhyML_path,
                   args = argumentsFit,
                   stdout = paste0( tempDir, "/std_fit.out" ),
                   stderr = paste0( tempDir, "/std_fit.err" ),
                   input = "\n" )

    if( any( grepl( pattern = "Type enter to exit|site lhood = -nan", x = readLines( con = paste0( tempDir, "/std_fit.out" ), warn = FALSE ) ) ) )
    {
      if( !is.null( log ) )
        cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Error occured during execution, check the log files." ), file = log )
      result[[ "success" ]] <- FALSE
      result[[ "error_message" ]] <- paste( c( base::readLines( con = paste0( tempDir, "/std_fit.out" ), warn = FALSE ) ), collapse = "\n" )
      if( !useTempDir ){
        # change back directory
        setwd( curDir )
      }
      return( result )
    }
    if( !is.null( log ) )
      cat( paste0( "\n# ", curDateTime(), " --- IgPhyML: Finished construction of trees" ), file = log )

    # update the list
    if( file.exists( paste0( tempDir, '/', basename( path ), "_igphyml_tree.txt" ) ) ){
      result[[ "topology" ]] <- list( tree = phytools::read.newick( file = paste0( tempDir, '/', basename( path ), "_igphyml_tree_gy94.txt" ), quiet = TRUE ),
                                      stats = base::readLines( con = paste0( tempDir, '/', basename( path ), "_igphyml_stats_gy94.txt" ), warn = FALSE ),
                                      stdout = base::readLines( con = paste0( tempDir, "/std_topology.out" ), warn = FALSE ),
                                      stderr = base::readLines( con = paste0( tempDir, "/std_topology.err" ), warn = FALSE ) )
      result[[ "fit" ]] <- list( tree = phytools::read.newick( file = paste0( tempDir, '/', basename( path ), "_igphyml_tree.txt" ), quiet = TRUE ),
                                 stats = base::readLines( con = paste0( tempDir, '/', basename( path ), "_igphyml_stats.txt" ), warn = FALSE ),
                                 stdout = base::readLines( con = paste0( tempDir, "/std_fit.out" ), warn = FALSE ),
                                 stderr = base::readLines( con = paste0( tempDir, "/std_fit.err" ), warn = FALSE ) )
    }
  }
  if( !useTempDir ){
    # change back directory
    setwd( curDir )
  }

  # return result
  return( result )
}