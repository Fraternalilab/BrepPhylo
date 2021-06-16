context( "test-clonelineage" )

library( BrepPhylo )
input <- read.csv( system.file( "extdata/clone_188.csv",
                                package = "BrepPhylo" ) )
outputFolder <- paste0( tempdir(), "/cloneLineage_test" )
dir.create( outputFolder, showWarnings = FALSE )

test_that( desc = "Input checks working?",
           code = {
             expect_error( object = BrepPhylo::cloneLineage( input = input,
                                                             outputFolder = outputFolder,
                                                             species = "Cattus_cattus" ),
                           "Specified species is not supported." )
             expect_error( object = BrepPhylo::cloneLineage( input = input,
                                                             outputFolder = outputFolder,
                                                             colourSpecification = list( coloursColumn = "Subclass" ),
                                                             sequenceColumn = "VDJ_SORTED",
                                                             cloneIdColumn = "CloneID___",
                                                             germlineIdColumn = "V.GENE.and.allele" ),
                           "One column or multiple columns specified missing." )
             expect_error( object = BrepPhylo::cloneLineage( input = input,
                                                             outputFolder = "/none/path",
                                                             colourSpecification = list( coloursColumn = "Subclass" ),
                                                             sequenceColumn = "VDJ_SORTED",
                                                             cloneIdColumn = "CloneID",
                                                             germlineIdColumn = "V.GENE.and.allele" ),
                           "Specified outputFolder does not exist." )
             })

test_that( desc = "Calculation working?",
           code = {
             BrepPhylo::cloneLineage( input = input,
                                      outputFolder = outputFolder,
                                      species = "Human",
                                      whichClones = 188,
                                      colourSpecification = list( colours = list( "M" = "#FF0000",
                                                                                  "D" = "#FF00FF",
                                                                                  "IgA1" = "#00FFFF",
                                                                                  "IgA2" = "#000088",
                                                                                  "IgG1" = "#7FFF88",
                                                                                  "IgG2" = "#7FFF00",
                                                                                  "IgG3" = "#682288",
                                                                                  "IgG4" = "#880000" ),
                                                                  coloursColumn = "Subclass",
                                                                  germlineColour = "#AAAAAA" ),
                                      sequenceColumn = "VDJ_SORTED",
                                      cloneIdColumn = "CloneID",
                                      germlineIdColumn = "V.GENE.and.allele",
                                      labelColumn = "Seq_ID",
                                      treeConstruction = list( "type" = "simple" ) )
      
             tree <- read.newick( file = paste0( outputFolder, "/clone_188/clone_188.tree" ) )
             expect_equal( object = tree$Nnode, 88 )
             expect_equal( object = tree$tip.label[ 79:89 ], c( "E50B_A1_147682_UID141490_RD40", "E50B_A1_104319_UID103810_RD29", "E50B_A1_130569_UID126657_RD9", 
                                                                "E50B_G1_112874_UID111301_RD15", "E50B_G_127476_UID123986_RD28", "E50B_G1_152709_UID145746_RD39",
                                                                "E50B_G2_155822_UID148353_RD38", "E50B_G1_104838_UID104264_RD31", "E50B_A1_106597_UID105803_RD21",
                                                                "E50B_G_120135_UID117612_RD33", "IGHV3-49*04" ) )
           })
