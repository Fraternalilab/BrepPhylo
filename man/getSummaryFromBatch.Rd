% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/support.R
\name{getSummaryFromBatch}
\alias{getSummaryFromBatch}
\title{Function to summarise distance-from-germline and CSR events from batch lineage analysis}
\usage{
getSummaryFromBatch(batchOutput)
}
\arguments{
\item{batchOutput}{list, output from \code{doBatchCloneAnalysis}.}
}
\value{
A list with two elements: 
\describe{
  \item{distances}{data.frame with the distance-from-germline calculated from any sequences in all the lineage trees included in the input.}
  \item{csr_events}{data.frame listing all class-switching events across all the lineages include in the input, and an estimate distance-from-germline at which such event takes place.}
}
}
\description{
This function processes the output from \code{doBatchCloneAnalysis} and merge entries representing individual lineage trees into two data frames, one containing distance-from-germline of all sequences in all lineages, and the other containing CSR events identified in all lineages.
}
