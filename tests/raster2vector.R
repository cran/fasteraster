library(fasteraster);
library(datasets);

load("dem4");

inp <- dem;
inp[220, 400] <- 0.3;
inp[221, 400] <- 0.35;
inp[222, 400] <- 0.4;
inp[220, 401] <- 0.45;
inp[221, 401] <- 0.45;
inp[222, 401] <- 0.5;
inp[220, 402] <- 0.55;
inp[221, 402] <- 0.5;
inp[222, 402] <- 0.6;

inp <- inp[150:300, 350:500];

res <- raster2vector(inp, 0, 1, 0.3, 1, FALSE);

image(inp, col = rev(grey.colors(100)), useRaster = TRUE)
plot(0, type = "l", xlim = c(0, nrow(inp)), ylim = c(0, ncol(inp)))
a <- lapply(res, function(x) lines(rbind(x, x[1,])))

inp <- volcano;

res <- raster2vector(volcano, 120, 200, 20);

image(inp, col = rev(grey.colors(100)), useRaster = TRUE)
plot(0, type = "l", xlim = c(0, nrow(inp)), ylim = c(0, ncol(inp)))
a <- lapply(res, function(x) lines(rbind(x, x[1,])))
