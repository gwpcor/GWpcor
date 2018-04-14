\name{gwpcor}
\alias{gwpcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographically weighted partial correlation
}
\description{
This function calculates the geographically weighted partial correlations between two variables given others. The function is designed from gwss function from the GWmodel package, and the cor2pcor function from the corpcor package.
}
\usage{
gwpcor(sdata, summary.locat, vars, method = c("pearson", "spearman"), kernel = "bisquare", adaptive = FALSE, bw, p = 2, theta = 0, longlat = F, dMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdata}{
a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp

}

  \item{summary.locat}{
a Spatial*DataFrame object for providing summary locations, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp

}

  \item{vars}{
a vector of variable names to be summarized

}
  \item{method}{
A character string indicating which partial correlation coefficient to compute. Currently, only "pearson" or "spearman" are accepted.

}
  \item{kernel}{
function chosen as follows:
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise

}
  \item{adaptive}{
if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)

}
  \item{p}{
the power of the Minkowski distance, default is 2, i.e. the Euclidean distance

}
  \item{theta}{
an angle in radians to rotate the coordinate system, default is 0

}
  \item{longlat}{
if TRUE, great circle distances will be calculated

}
  \item{dMat}{
a pre-specified distance matrix, it can be calculated by the function gw.dist

}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{SDF }{
a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package “sp”) with local covariances, local correlations (Pearson's), local correlations (Spearman's),  p-values of local correlations (Pearson's), p-values of local correlations (Spearman's), local partial correlations (Pearson's), local partial correlations (Spearman's), p-values of local partial correlations (Pearson's), and p-values of local partial correlations (Spearman's).
}
  \item{vars}{
  Names of variables used for the calculation.
}
   \item{kernel}{
   The name of kernel used for the calculation.
}
   \item{adaptive}{
   TRUE/FALSE
   }

   \item{bw}{
   The bandwidth size used for the calculation.
   }

}
\references{
Percival J. and Tsutsumida N. (2017) Geographically weighted partial correlation for spatial analysis, GI_forum, Issue 1, 36-43, URL http://dx.doi.org/10.1553/giscience2017_01_s36

Isabella Gollini, Binbin Lu, Martin Charlton, Christopher Brunsdon, Paul Harris (2015). GWmodel: An R Package for Exploring Spatial Heterogeneity Using Geographically Weighted Models. Journal of Statistical Software, 63(17), 1-50. URL http://www.jstatsoft.org/v63/i17/.

Binbin Lu, Paul Harris, Martin Charlton, Christopher Brunsdon (2014). The GWmodel R package: further topics for exploring spatial heterogeneity using geographically weighted models. Geo-spatial Information Science, 17(2), 85-101. URL http://dx.doi.org/10.1080/10095020.2014.917453
}
\author{
N. Tsutsumida and J. Percival.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#NOTE: This example only shows how to implement gwpcor using sample data (meuse) in sp package. Results do not suggest any meanings.

#import data from sp package, and create a SpatialPointsDataFrame
data("meuse")
coordinates(meuse) = ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")

#implement gw partial correlation as an example
#the bandwidth is arbitrary.
res <- gwpcor(sdata = meuse, vars = c("cadmium","copper", "zinc"), method = "pearson", kernel = "bisquare",adaptive = TRUE,bw = 50,longlat = FALSE)

#show results
spplot(res$SDF, c("corr_cadmium.copper","corr_cadmium.zinc","corr_copper.zinc")) #GW-cor
spplot(res$SDF, c("pcorr_cadmium.copper","pcorr_cadmium.zinc","pcorr_copper.zinc")) #GW-pcor
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
