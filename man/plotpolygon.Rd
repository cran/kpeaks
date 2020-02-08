\name{plotpolygon}
\alias{plotpolygon}
\title{
Plot Frequency Polygons
}
\description{
Plots the frequency polygon and histogram of a feature with some options.
}
\usage{
plotpolygon(x, nbins, ptype, bcol = "gray", pcol = "blue")
}
\arguments{
\item{x}{
a numeric vector containing the observations of a feature, or a numeric matrix when \code{ptype} is \option{sp}.
}
\item{nbins}{
an integer for the number of classes in the frequency polygon.
}

\item{bcol}{
a string for the color of bins. Default is \option{gray}.
}

\item{pcol}{
a string for the color of polygon lines. Default is \option{blue}.
}

\item{ptype}{
a string specifying the type of plot. Use \option{p} for plotting the polygon only or \option{ph} for plotting the polygon with the histogram. Default is \option{sp} for the scatterplots between the pairs of features and the polygons on the diagonal panel.
}

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
\code{\link{genpolygon}}
}
\examples{
# plot the frequency polygon of the 2nd feature in x5p4c data set
data(x5p4c)
hvals <- genpolygon(x5p4c[,2], binrule="usr", nbins=20)

# plot the frequency polygon of the 2nd feature in x5p4c data set
plotpolygon(x5p4c[,2], nbins=hvals$nbins, ptype="p")

# plot the histogram and frequency polygon of the 2nd feature in x5p4c data set
plotpolygon(x5p4c[,2], nbins=hvals$nbins, ptype="ph", bcol="orange", pcol="blue")

# plot the pairwise scatter plots of the features in x5p4c data set
pairs(x5p4c, diag.panel=plotpolygon, upper.panel=NULL, cex.labels=1.5)

# plot the histogram and frequency polygon of Petal.Width in iris data set 
data(iris)
hvals <- genpolygon(iris$Petal.Width, binrule="doane")
plotpolygon(iris$Petal.Width, nbins=hvals$nbins, ptype="ph")
}
\concept{frequency polygons}
\concept{histograms}
\concept{scatter plots}
