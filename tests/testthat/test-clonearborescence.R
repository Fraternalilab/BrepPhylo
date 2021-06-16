context( "test-clonearborescence" )

library( BrepPhylo )
clone <- read.csv( system.file( "extdata/clone_188.csv",
                                package = "BrepPhylo" ) )
tree <- system.file( "extdata/clone_188.tree",  package = "BrepPhylo" )
outputFolder <- paste0( tempdir(), "/cloneArborescence_test" )
dir.create( outputFolder, showWarnings = FALSE )

test_that( desc = "Input checks working?",
           code = {
             expect_error( object = BrepPhylo::cloneArborescence( clone = clone,
                                                                  cloneID = 188,
                                                                  columnSeqID = "Seq_ID",
                                                                  columnSubclass = "Subclass",
                                                                  tree = tree,
                                                                  outputFolder = "/no/path/ready" ),
                           "Output folder apparently does not exist." )
             expect_error( object = BrepPhylo::cloneArborescence( clone = clone,
                                                                  cloneID = 188,
                                                                  columnSeqID = "Seq_ID___",
                                                                  columnSubclass = "Subclass",
                                                                  tree = tree,
                                                                  outputFolder = outputFolder ),
                           "Column names specified were not found." )
           })

test_that( desc = "Calculation working?",
           code = {
             arborescence <- cloneArborescence( clone = clone,
                                                cloneID = "188",
                                                columnSeqID = "Seq_ID",
                                                columnSubclass = "Subclass",
                                                tree = tree,
                                                outputFolder = outputFolder )
             
             expect_equal( object = arborescence$table[ 1:7, 1 ], c( 89, 40, 37, 43, 37, 45, 61 ) )
             expect_equal( object = arborescence$table[ 1, "startLabel" ], "IGHV3-49*04" )
           })
