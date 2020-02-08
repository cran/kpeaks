\name{findpolypeaks}
\alias{findpolypeaks}
\title{Find the Peaks of a Frequency Polygon}
\description{Frequency polygons are graphics to reveal the shapes of data distributions as histograms do. The peaks of frequency polygons are required in several data mining applications. \code{findpolypeaks} finds the peaks in a frequency polygon by using the frequencies and middles values of the classes of it.}
\usage{
findpolypeaks(xm, xc, tcmethod, tc)
}
\arguments{
  \item{xm}{a numeric vector contains the middle values of the classes of the frequency polygon (or the bins of a histogram).}
  \item{xc}{an integer vector contains the frequencies of the classes of the frequency polygon.}
  \item{tcmethod}{a string represents the threshold method to discard the empty and the small bins whose frequencies are smaller than a threshold frequency value. Default method is \option{usr}. Alternatively, the methods given below can be used to compute a threshold frequency value using the descriptive statistics of the frequencies in \code{xc}.
    \itemize{
    \item \option{sd1} and \option{sd2} use the standard deviation.
    \item \option{q1} uses the first quartile (Q1).
    \item \option{iqr} uses the interquartile range (IQR).
    \item \option{avg} uses the arithmetic mean.
    \item \option{min} and \option{min2} use the minimum.
    \item \option{log2} uses the two-base logarithm of \var{n}, vector size.
    \item \option{usr} uses a user-specified value.}}
 \item{tc}{an integer which is used as the threshold frequency value for discarding the empty and small height classes in the frequency polygon. Default value is 1 if the threshold option \option{usr} is chosen. Depending on the selected methods, the value of \code{tc} equals to: 
    \itemize{
    \item one standart deviation with the method \option{sd1}, 
    \item one quarter of the standart deviation with the method \option{sd2}, 
    \item the first quartile with the method \option{q1},
    \item one quarter of the interquartile range with method \option{iqr},
    \item 10\% of the arithmetic mean with the method \option{avg}, 
    \item the minimum value with the method \option{min}, 
    \item two times of minimum with method \option{min2},
    \item two-base logarithm of the number of classes divided by ten with the method \option{log2},
    \item an arbitrary number specified with the method \option{usr}.}}
}
\details{
The peaks are determined after removing the empty and small height classes whose frequencies are below the chosen threshold frequency. Default threshold value is 1 that means that all the classes which have frequencies of 0 and 1 are removed in the input vectors \code{xm} and \code{xc}. 
}
\value{
 \item{pm}{a data frame with two columns which are \var{pvalues} and \var{pfreqs} containing the middle values and frequencies of the peaks which determined in the frequency polygon, respectively.}
 \item{np}{an integer representing the number of peaks in the frequency polygon.}
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
  \code{\link{genpolygon}},
  \code{\link{rmshoulders}}
}
\examples{
data(x5p4c)
# Using a user-specified number of bins, build the frequency polygon of p2 in the data set x5p4c
hvals <- genpolygon(x5p4c$p2, binrule="usr", nbins=20)
plotpolygon(x5p4c$p2, nbins=hvals$nbins, ptype="ph")

# Find the peaks in the frequency polygon by using the threshold method min
resfpp1 <- findpolypeaks(hvals$mids, hvals$freqs, tcmethod="min")
print(resfpp1)

# Find the peaks in the frequency polygon by using the threshold equals to 5
resfpp2 <- findpolypeaks(hvals$mids, hvals$freqs, tcmethod="usr", tc=5)
print(resfpp2)

data(iris)
# By using Doane rule, build the frequency polygon of the 4th feature in the data set iris
hvals <- genpolygon(iris[,4], binrule="doane")
plotpolygon(iris[,4], nbins=hvals$nbins, ptype="p")

#Find the peaks in the frequency polygon by using the threshold method avg
resfpp3 <- findpolypeaks(hvals$mids, hvals$freqs, tcmethod="avg")
print(resfpp3)
}
\concept{frequency polygons}
\concept{frequency distribution}
\concept{peaks in a distribution}

\keyword{cluster}
