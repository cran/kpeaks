\name{rmshoulders}
\alias{rmshoulders}
\title{
Shoulders Removal in Frequency Polygons
}
\description{
Removes the shoulders around the main peaks in a frequency polygon.
}

\usage{
rmshoulders(xm, xc, trmethod, tv)
}

\arguments{
\item{xm}{
a numeric vector containing the middle values of peaks of a frequency polygon.
}
\item{xc}{
an integer vector containing the frequencies of peaks of a frequency polygon.
}
\item{trmethod}{a string representing the type of shoulders removal option for computing a threshold value. Default method is \option{usr}. The alternatives are \option{sd}, \option{q1}, \option{iqr}, \option{avg} and \option{med}. These methods compute the threshold distance value using some statistics of the distances between the middle values of two successive peaks in the vector \code{xm}.

  \itemize{
    \item \option{sd} uses the standard deviation.
    \item \option{q1} uses the first quartile (Q1).
    \item \option{q3} uses the third quartile (Q3).
    \item \option{iqr} uses the interquartile range (IQR).
    \item \option{avg} uses the arithmetic mean.
    \item \option{med} uses the median.
    \item \option{usr} uses a user-specified number.
  }
}
\item{tv}{a numeric value to be used as the threshold distance for deciding the shoulders. Default threshold is 1 if the removal method \option{usr} is chosen. Depending on the selected removal method \code{tv} equals to:
  \itemize{
    \item one standart deviation if \code{trmethod} is \option{sd},
    \item the first quartile if \code{trmethod} is \option{q1}, 
    \item the third quartile if \code{trmethod} is \option{q3}, 
    \item one quarter of the interquartile range if \code{trmethod} is \option{iqr}, 
    \item the arithmetic mean if \code{trmethod} is \option{avg},
    \item the median if \code{trmethod} is \option{med},
    \item a user-specified number if \code{trmethod} is \option{usr}. 
  }
 }
}

\details{
Literally speaking, a \dfn{shoulder peak} or shortly \dfn{shoulder} is a secondary peak in a close location before or after the main peak of a mountain. In a frequency polygon, a shoulder is a smaller peak that is quite close to a higher peak resulting a non-obvious valley between them. Shoulders may occur randomly due to some reasons such as random noises or selecting higher number of classes in histogram building etc. Usually, it is desired to remove them from the peaks vector of a frequency polygon. In '\pkg{kpeaks}', a peak considered as a shoulder when its height is smaller than the height of its neighbor peak and its distance to its neighbor is also lower than a threshold distance value. In order to compute a threshold distance value, here, we propose to use seven options as listed in the section \sQuote{arguments}. The options \code{q1} and \code{iqr} can be applied to remove the minor shoulders that are very near to the main peaks while \code{q3} is recommended to eliminate the substantial shoulders in the processed frequency polygon. The remaining options may be more efficient for removing the moderate shoulders.
}

\value{
 \item{pm}{a data frame with two columns whose names are \var{pvalues} and \var{pfreqs} for the middle values and the frequencies of the peaks after removal process, respectively.}
 \item{np}{an integer representing the number of peaks after removal of the shoulders.}
}

\note{
The function \code{rmshoulders} normally should be called with the input values that are returned by the function \code{\link{findpolypeaks}}.
}

\author{
Zeynel Cebeci, Cagatay Cebeci
}

\references{
Cebeci, Z. & Cebeci, C. (2018). "A novel technique for fast determination of K in partitioning cluster analysis", \emph{Journal of Agricultural Informatics}, 9(2), 1-11. 
doi: \href{https://doi.org/10.17700/jai.2018.9.2.442}{10.17700/jai.2018.9.2.442}.

Cebeci, Z. & Cebeci, C. (2018). "kpeaks: An R Package for Quick Selection of K for Cluster Analysis", In \emph{2018 International Conference on Artificial Intelligence and Data Processing (IDAP)}, IEEE. 
doi: \href{https://doi.org/10.1109/IDAP.2018.8620896}{10.1109/IDAP.2018.8620896}.
}

\seealso{
\code{\link{findpolypeaks}},
\code{\link{plotpolygon}},
\code{\link{genpolygon}}
}

\examples{
# Build a data vector with three peaks
x1 <-rnorm(100, mean=20, sd=5)
x2 <-rnorm(50, mean=50, sd=5)
x3 <-rnorm(150, mean=90, sd=10)
x <- c(x1,x3,x2)

# generate the frequency polygon and histogram of x by using Doane rule
hvals <- genpolygon(x, binrule="doane")
plotpolygon(x, nbins=hvals$nbins, ptype="p")

# find the peaks in frequency polygon of x by using the default threshold frequency
resfpp <- findpolypeaks(xm=hvals$mids, xc=hvals$freqs)
print(resfpp)

# remove the shoulders with the threshold distance option 'avg'
resrs <- rmshoulders(resfpp$pm[,1], resfpp$pm[,2], trmethod = "avg")
print(resrs)

# remove the shoulders with the threshold distance option 'iqr'
resrs <- rmshoulders(resfpp$pm[,1], resfpp$pm[,2], trmethod = "iqr")
print(resrs)

data(x5p4c)
# plot the frequnecy polygon and histogram of p2 in x5p4c data set 
hvals <- genpolygon(x5p4c$p2, binrule="usr", nbins=30)
plotpolygon(x5p4c$p2, nbins=hvals$nbins, ptype="ph")

# find the peaks in frequency polygon of p2 
resfpp <- findpolypeaks(xm=hvals$mids, xc=hvals$freqs, tcmethod = "min")
print(resfpp)

# remove the shoulders with threshold distance option 'q1'
resrs <- rmshoulders(resfpp$pm[,1], resfpp$pm[,2], trmethod = "q1")
print(resrs)

\dontrun{
data(iris)
# plot the frequency polygon and histogram of Petal.Length in iris data set 
# by using a user-defined class number 
hvals <- genpolygon(iris$Petal.Length, binrule="usr", nbins=30)
plotpolygon(iris$Petal.Length, nbins=hvals$nbins, ptype="p")

# find the peaks in frequency polygon of Petal.Length with default 
# threshold frequency value
resfpp <- findpolypeaks(xm=hvals$mids, xc=hvals$freqs)
print(resfpp)

# remove the shoulders with threshold option 'med'
resrs <- rmshoulders(resfpp$pm[,1], resfpp$pm[,2], trmethod = "med")
print(resrs)
}
}
\concept{frequency polygons}
\concept{peaks in a frequency polygon}
