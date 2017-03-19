library(fasteraster);
library(datasets);

polygons <- raster2vector(volcano, 120, 200, 20, 1);
image(volcano, col = rev(grey.colors(100)), useRaster = TRUE)
plot(0, type = "l", xlim = c(0, nrow(volcano)), ylim = c(0, ncol(volcano)))
a <- lapply(polygons, function(x) lines(rbind(x, x[1,])))

zones <- rasterZoneAnalyzer(volcano, 120, 200, 20);
a <- text(zones[ , 3], zones[ , 4], labels = zones[ , 2]);

