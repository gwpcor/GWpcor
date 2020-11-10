\name{gwpcor}
\alias{gwpcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographically Weighted Correlation and Partial Correlation
}
\description{
This function calculates the geographically weighted correlation and partial correlation between two variables given others. The function is designed by the gwss function from the GWmodel package, and the cor2pcor function from the corpcor package.
}
\usage{
gwpcor(sdata, res_dp, vars, method = c("pearson", "spearman"), 
kernel = "bisquare", adaptive = FALSE, bw, dMat, foreach = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdata}{
a Spatial*DataFrame (i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp),
or a sf object. 

}

  \item{res_dp}{
A Spatial*DataFrame object for providing summary locations, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp, or a sf object. 

}

  \item{vars}{
A vector of variable names to be used for the analysis.

}
  \item{method}{
A character string indicating which correlation and partial correlation coefficients to compute. "pearson" or "spearman" are accepted.

}
  \item{kernel}{
function chosen as follows:
gaussian: wgt = exp( -0.5 * (vdist / bw)^2);
exponential: wgt = exp( -vdist / bw);
bisquare: wgt = ( 1 -(vdist / bw)^2)^2 if vdist < bw, wgt = 0 otherwise;
tricube: wgt = (1 -(vdist / bw)^3)^3 if vdist < bw, wgt = 0 otherwise;
boxcar: wgt = 1 if dist < bw, wgt = 0 otherwise

}
  \item{adaptive}{
if TRUE, an adaptive kernel where the bandwidth (bw) corresponds to the proportion of the number of nearest neighbours (i.e. adaptive distance) is employed. The default is FALSE, where a fixed kernel is employed (bandwidth is a fixed distance).

}
  \item{bw}{
Bandwidth size. If adaptive kernel, bw should be the proportion of the number of nearest neighbours (0 < bw <=1). For fixed kernel, the Euclid distance.

}

  \item{dMat}{
A pre-specified distance matrix, it can be calculated by the function st_distance().

}
  \item{foreach}{
Whether parallel computation is implemented or not.

}

}

\value{
  \item{SDF }{
  A SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package “sp”) when the input is Spatial*DataFrame or a sf class object when input is sf, with local covariances, local correlations (Pearson's), local correlations (Spearman's),  p-values of local correlations (Pearson's), p-values of local correlations (Spearman's), local partial correlations (Pearson's), local partial correlations (Spearman's), p-values of local partial correlations (Pearson's), and p-values of local partial correlations (Spearman's).
}
  \item{vars}{
  Names of variables used for the calculation.
}
   \item{kernel}{
   The name of kernel used for the calculation.
}
   \item{adaptive}{
   Whether aadaptive kernel is employed or not (TRUE/FALSE),
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
Tsutsumida N. and Percival J.
}

\examples{
#NOTE: This example only shows how to implement gwpcor using sample data (meuse) in sp package. 
#Results do not suggest any meanings.

#import data from sp package
library(sp)
library(sf)
data(meuse, package = "sp") 
meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992)

#implement gwpcor as an example
#the bandwidth is arbitrary.
res <- gwpcor(sdata = meuse_sf, vars = c("cadmium","copper", "zinc"),
method = "pearson", kernel = "bisquare",adaptive = TRUE, 
bw = 0.25, foreach = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "SPATIAL" }% use one of  RShowDoc("KEYWORDS")
\keyword{ "MULTIVALIATE" }% __ONLY ONE__ keyword per line
