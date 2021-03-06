% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastliclust.R
\name{crop}
\alias{crop}
\title{Crop the result matrix}
\usage{
crop(flInput)
}
\arguments{
\item{flInput}{}
}
\value{
a flInput-format result for processing with \code{\link{toHclust}}
}
\description{
From a \code{\link{fastLiclust}}() result, crops out the actually useful part.
In addition, it checks for disconnected subgraphs and connects them together
with the highest observed distance.
}
