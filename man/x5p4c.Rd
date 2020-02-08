\name{x5p4c}
\alias{x5p4c}
\docType{data}
\title{
Synthetic Data Set contains 5 Variables and 4 Clusters 
}
\description{
A synthetically created data frame consists of five continous variables forming four clusters.
}
\usage{data(x5p4c)}
\format{A data frame with 400 rows and 5 numeric variables:
 \describe{
  \item{p1}{a continous variable with one mode}
  \item{p2}{a continous variable with four modes}
  \item{p3}{a continous variable with two modes}
  \item{p4}{a continous variable with three modes}
  \item{p5}{a continous variable with two modes}
 }
}
\note{
The data set \code{x5p4c} is recommended to use in comparing the performances of the internal validity indexes in cluster analysis.
}
\examples{
data(x5p4c)
# descriptive statistics of the variables
summary(x5p4c)
# plot the histogram of the variable p2
hist(x5p4c$p2, breaks=15)
# scatter plots of the variable pairs
pairs(x5p4c)
}
\concept{synthetic datasets}
\concept{clustering}
\concept{cluster analysis}
\keyword{Basics|datasets}
