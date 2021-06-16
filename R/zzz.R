# "onLoad" functionality: Set up python interface
scipy <- NULL
.onLoad <- function( libname, pkgname )
{
  # use superassignment to update global reference to scipy
  # this delays the loading of python a bit
  scipy <<- reticulate::import( "scipy", delay_load = TRUE )
}