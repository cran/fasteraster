\name{fasteraster}
\alias{fasteraster}
\docType{package}
\title{
  Raster Image Processing and Vector Recognition
}
\description{
  If there is a need to recognise edges on a raster image or a bitmap or any kind of a matrix, one can find packages
  that does only 90 degrees vectorization. Typically the nature of artefact images is linear and can be vectorized in much more
  efficient way than draw a series of 90 degrees lines. The fasteraster package does recognition of lines using only one pass.
  It also allows to calculate weight and the mass centers for the recognized zones or polygons.
}
\details{
  Use raster2vector to vectorize input matrix.
  Use rasterZoneAnalyzer to break the input matrix by zones and get their weight and mass center.
}
\author{
Andy Bosyi <andy@bosyi.com>.

Maintainer: Andy Bosyi <andy@bosyi.com>
}
\keyword{ package }
\seealso{
  http://bosyi.com/craft/
}
\examples{
  \dontrun{
  library(fasteraster);
  library(datasets);

  polygons <- raster2vector(volcano, 120, 200, 20, 1);
  image(volcano, col = rev(grey.colors(100)), useRaster = TRUE)
  plot(0, type = "l", xlim = c(0, nrow(volcano)), ylim = c(0, ncol(volcano)))
  a <- lapply(polygons, function(x) lines(rbind(x, x[1,])))

  zones <- rasterZoneAnalyzer(volcano, 120, 200, 20);
  a <- text(zones[ , 3], zones[ , 4], labels = zones[ , 2]);
  }
}
