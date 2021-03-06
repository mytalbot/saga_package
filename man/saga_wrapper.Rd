% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saga_wrapper.R
\name{saga_wrapper}
\alias{saga_wrapper}
\title{SAGA Wrapper}
\usage{
saga_wrapper(smplpath, showModel = 1, doGSEA = 0)
}
\arguments{
\item{smplpath}{sample path}

\item{showModel}{Can be 0 or 1. If set to 1 (default) PCA plot (+samples) will be shown.}

\item{doGSEA}{Can be 0 or 1. If set to 1, GSEA analysis will be added to the analysis.}
}
\value{
\code{output} optimized SVM model for classification.
}
\description{
The \code{saga_wrapper} function integrates all SAGA functions into one. With this function you are not as flexible as
with the single ones but you can run the analysis/ classification in one step.
However, this requires a complete SampleInformatin.txt file. For GSEA analysis the SampleInformation file must also fulfill
the neccessary criteria. See saga_vignette for more information.
}
\details{
\code{saga_wrapper}
}
