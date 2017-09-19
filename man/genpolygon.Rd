\name{genpolygon}
\alias{genpolygon}
\title{Generate the Classes to Build a Frequency Polygon}
\description{Constructs the histogram of a feature by using a selected binning rule, returns the middle values and frequencies of classes for further works on the frequency polygon.}
\usage{
genpolygon(x, binrule, nbins, disp = FALSE)
}
\arguments{
  \item{x}{a numeric vector containing the observations for a feature.}
  \item{binrule}{name of the rule in order to compute the number of bins to build the histogram.}
  \item{nbins}{an integer representing the number of bins which is computed by using the selected binning rule. Default rule is \option{sturges}. Depending on the selected rule \code{nbins} equals to (In the formulae, \code{n} is the number of observations.):
    \itemize{
      \item \code{floor(sqrt(n))} if the rule is \option{sqr},
      \item \code{ceiling(1+log(n, 2))} if the rule is \option{sturges},
      \item \code{ceiling(1+3.332*log(n, 10))} if the rule is \option{huntsberger},
      \item \code{ceiling(5*log(n, 10))} if the rule is \option{bc},
      \item \code{ceiling(n^(1/3))} if the rule is \option{cencov},
      \item \code{ceiling(2*n^(1/3))} if the rule is \option{rice},
      \item \code{ceiling((2*n)^(1/3))} if the rule is \option{ts},
      \item \code{ceiling(((max(x)-min(x))/(3.5*sqrt(var(x))*n^(-1/3))))} if the rule is \option{scott},
      \item \code{ceiling(((max(x)-min(x))/(2*IQR(x)*n^(-1/3)))} if the rule is \option{fd},
      \item \code{ceiling(1+log(n,2)+log(1+abs(skewness(x))/(6*(n-2)/((n+1)*(n+3))^0.5),2))} if the rule is \option{doane},
      \item \code{ceiling(log(n)/2*pi)} if the rule is \option{cebeci},
      \item a user-specified integer if the rule is \option{usr}.
      }
  }
  \item{disp}{a logical value should be set to \code{TRUE} to display the histogram.}
}
\details{
According to Hyndman (1995), Sturges's rule was the first rule to calculate \var{k}, the number of classes to build a histogram. Most of the statistical packages use this simple rule for determining the number of classes in constructing histograms. Brooks & Carruthers (1953) proposed a rule using \eqn{log_{10}}{log10} instead of \eqn{log_{2}}{log2} giving always larger \var{k} when compared to Sturges's rule.  The rule by Huntsberger (1962) yields nearly equal result to those of Sturges's rule. These two rules work well if \var{n} is less than 200. Scott (1992) argued that Sturges's rule leads to generate oversmoothed histograms in case of large number of \var{n}. In his rule, Cencov (1962) used the cube root of \var{n} simply. This rule was followed by its extensions, i.e., Rice rule and Terrell & Scott (1985) rule.  When compared to the others, the square root rule produces larger \var{k} (Davies & Goldsmith, 1980).\cr\cr
Most of the rules simply include only \var{n} as the input argument. On the other hand, the rules using variation and shape of data distributions can provide more optimal \var{k} values. For instance, Doane (1976) extended the Sturges's rule by adding the standardized skewness in order to overcome the problem with non-normal distributions need more classes. In order to estimate optimal \var{k} values, Scott (1979) added the standard deviation to his formula. Freedman and Diaconis (1981) proposed to use the interquartile range (IQR) statistic which is less sensitive to outliers than the standard deviation. In a study on unsupervised discretization methods, Cebeci & Yildiz (2017) tested a binning rule formula based on the ten-base logarithm of \var{n} divided by \code{2*pi}. They also argued that the rules Freedman-Diaconis and Doane were slightly performed better than the other rules based on the training model accuracies on a chicken egg quality traits dataset. Therefore, using the above mentioned rules may be more effective in determining the peaks of a frequency polygon. 
}
\value{
 \item{xm}{a numeric vector containing the middle values of bins.}
 \item{xc}{an integer vector containing the frequencies of the bins.}
 \item{nbins}{an integer containing the number of bins to build the histogram.}
}

\references{
Brooks C E P & Carruthers N (1953). Handbook of statistical methods in meteorology. H M Stationary Office, London.

Cebeci Z & Yildiz F (2017). Unsupervised discretization of continuous variables in a chicken egg quality traits dataset. \emph{Turk. J Agriculture-Food Sci. & Tech.} 5(4): 315-320.

Cencov N N (1962). Evaluation of an unknown distribution density from observations. \emph{Soviet Mathematics} 3: 1559-1562.

Davies O L & Goldsmith P L (1980). Statistical methods in research and production. 4th edn, Longman: London.

Doane D P (1976). Aesthetic frequency classification. \emph{American Statistician} 30(4):181-183.

Freedman D & Diaconis P (1981). On the histogram as a density estimator: L2 Theory. \emph{Zeit. Wahr. ver. Geb.} 57(4):453-476.

Hyndman R J (1995). The problem with Sturges rule for constructing histograms. \url{http://robjhyndman.com/papers/sturges.pdf}. 

Huntsberger D V (1962). Elements of statistical inference. London: Prentice-Hall.

Scott D W (1992). Multivariate density estimation: Theory, Practice and Visualization. John Wiley & Sons: New York.

Sturges H (1926). The choice of a class-interval. \emph{J Amer. Statist. Assoc.} 21(153):65-66.

Terrell G R & Scott D W (1985). Oversmoothed nonparametric density estimates. \emph{J Amer. Statist. Assoc.} 80(389):209-214. 
}

\author{Zeynel Cebeci, Cagatay Cebeci}

\seealso{
 \code{\link{findk}}, 
 \code{\link{findpolypeaks}}, 
 \code{\link{plotpolygon}}
}

\examples{
x <- rnorm(n=100, mean=5, sd=0.5)
# Construct the histogram of x according to the Sturges rule with no display
hvals <- genpolygon(x, binrule = "sturges")
print(hvals)

# Plot the histogram of x by using the user-specified number of classes
hvals <- genpolygon(x, binrule = "usr", nbins = 20, disp = TRUE)
print(hvals)

# Plot the histogram of the second feature in iris dataset 
# by using the Freedman-Diaconis (fd) rule
data(iris)
hvals <- genpolygon(iris[,2], binrule = "fd", disp = TRUE)
print(hvals)
}
\concept{histograms}
\concept{frequency polygons}