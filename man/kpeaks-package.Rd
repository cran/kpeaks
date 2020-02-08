\name{kpeaks-package}
\alias{kpeaks-package}
\docType{package}
\title{Determination of K Using Peak Counts of Features for Clustering}
\description{The input argument \var{k}, represents the number of clusters is needed to start all the partitioning clustering algorithms. In unsupervised learning applications, an optimal value of this argument is widely determined by using the internal validity indexes. Since these indexes suggest a \var{k} value which is computed on the clustering results obtained with several runs of a clustering algorithm, they are computationally expensive. On the contrary, the package '\pkg{kpeaks}' enables to estimate \var{k} before running any clustering algorithm. It is based on a simple novel technique using the descriptive statistics of peak counts of the features in a dataset.}

\details{
The package '\pkg{kpeaks}' contains five functions and one synthetically created dataset for testing purposes. In order to suggest an estimate of \var{k}, the function \code{\link{findk}} internally calls the functions \code{\link{genpolygon}} and \code{\link{findpolypeaks}}, respectively. The frequency polygons can be visually inspected by using the function \code{\link{plotpolygon}}. Using the function \code{\link{rmshoulders}} is recommended to flatten or remove the the shoulder peaks around the main peaks of a frequency polygon, if any.
}

\author{Zeynel Cebeci, Cagatay Cebeci}

\references{
Cebeci, Z. & Cebeci, C. (2018). "A novel technique for fast determination of K in partitioning cluster analysis", \emph{Journal of Agricultural Informatics}, 9(2), 1-11. 
doi: \href{https://doi.org/10.17700/jai.2018.9.2.442}{10.17700/jai.2018.9.2.442}.

Cebeci, Z. & Cebeci, C. (2018). "kpeaks: An R Package for Quick Selection of K for Cluster Analysis", In \emph{2018 International Conference on Artificial Intelligence and Data Processing (IDAP)}, IEEE. 
doi: \href{https://doi.org/10.1109/IDAP.2018.8620896}{10.1109/IDAP.2018.8620896}.
}

\seealso{
\code{\link{findk}},
\code{\link{findpolypeaks}},
\code{\link{genpolygon}},
\code{\link{plotpolygon}},
\code{\link{rmshoulders}}
}

\concept{finding number of clusters}
\concept{partitioning-based clustering}

\keyword{cluster}