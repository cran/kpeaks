\name{kpeaks-package}
\alias{kpeaks-package}
\docType{package}
\title{Determination of K Using Peak Counts of Features for Clustering}
\description{The input argument \var{k} that represents the number of clusters is needed to start all of the partitioning clustering algorithms. In unsupervised learning applications, an optimal value of this argument is widely determined by using the internal validity indexes. Since these indexes suggest a \var{k} value which is computed on the clustering results after several runs of a clustering algorithm, they are computationally expensive. On the contrary, '\pkg{kpeaks}' enables to estimate \var{k} before running any clustering algorithm. It is based on a simple novel technique using the descriptive statistics of peak counts of the features in a dataset.}

\details{
The package '\pkg{kpeaks}' contains five functions and one synthetically created dataset for testing purposes. In order to suggest an estimate of \var{k}, the function \code{\link{findk}} internally calls the functions \code{\link{genpolygon}} and \code{\link{findpolypeaks}}, respectively. The frequency polygons can be visually inspected by using the function \code{\link{plotpolygon}}. Using \code{\link{rmshoulders}} is recommended to flatten or remove the the shoulder peaks around the main peaks of a frequency polygon, if any.
}

\author{Zeynel Cebeci, Cagatay Cebeci}

\seealso{
\code{\link{findk}},
\code{\link{findpolypeaks}},
\code{\link{genpolygon}},
\code{\link{plotpolygon}},
\code{\link{rmshoulders}}
}

\concept{finding number of clusters}
\concept{partitioning-based clustering}

\keyword{Statistics|cluster}