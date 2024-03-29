% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellenONE_metadata_map.R
\name{link_cellenONE_Raw}
\alias{link_cellenONE_Raw}
\title{Link CellenONE metadata}
\usage{
link_cellenONE_Raw(QQC, cells_file)
}
\arguments{
\item{QQC}{A QuantQC object}

\item{allCells}{A named list of file paths to cellenONE isolation files}
}
\value{
A QQC object with \code{QQC at meta.data} slot and \code{cellenONE.meta} slot linking raw data to meta data from cellenONE
}
\description{
This function takes a QuantQC (QQC) object and file paths to cellenONE isolation files
as a named list.
}
\examples{
link_cellenONE_Raw(TestSamples,SortData_3_7)
}
