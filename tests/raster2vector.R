library(fasteraster);
library(datasets);

#load("dem5");

#inp <- opacities;

#inp <- inp[480:520, 90:180];

#res <- raster2vector(inp, 0, 1, 0.1, 1, FALSE);

#image(inp, col = rev(grey.colors(100)), useRaster = TRUE)
#plot(0, type = "l", xlim = c(0, nrow(inp)), ylim = c(0, ncol(inp)))
#a <- lapply(res, function(x) lines(rbind(x, x[1,])))

res <- raster2vector(volcano, 120, 200, 20);

image(volcano, col = rev(grey.colors(100)), useRaster = TRUE)
plot(0, type = "l", xlim = c(0, nrow(volcano)), ylim = c(0, ncol(volcano)))
a <- lapply(res, function(x) lines(rbind(x, x[1,])))
