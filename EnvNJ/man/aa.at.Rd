% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancillary.R
\name{aa.at}
\alias{aa.at}
\title{Residue Found at the Requested Position}
\usage{
aa.at(at, target)
}
\arguments{
\item{at}{the position in the primary structure of the protein.}

\item{target}{a character string specifying the sequence of the protein of interest.}
}
\value{
Returns a single character representing the residue found at the indicated position in the indicated protein.
}
\description{
Returns the residue found at the requested position.
}
\examples{
aa.at(2, 'MFSQQQRCVE')
}
\seealso{
aa.comp()
}
