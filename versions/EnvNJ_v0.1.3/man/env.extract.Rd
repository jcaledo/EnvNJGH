% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/envcoding.R
\name{env.extract}
\alias{env.extract}
\title{Sequence Environment Around a Given Position}
\usage{
env.extract(seq, c, r)
}
\arguments{
\item{seq}{a string protein sequence.}

\item{c}{center of the environment.}

\item{r}{radius of the environment.}
}
\value{
Returns a  a strings with the extracted environment sequence.
}
\description{
Extracts the sequence environment around a given position.
}
\examples{
env.extract('ARGGQQQCATSYV', c = 8,  r = 2)
}
\references{
Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
}
\seealso{
env.matrices(), env.sp()
}
