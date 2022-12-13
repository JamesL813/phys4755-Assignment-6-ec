library("scatterplot3d") # load

tetra <- read.csv("data/xyz4.data", header = TRUE, sep = " ")
png(filename = "img/tetra.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(tetra$x, tetra$y, tetra$z,
    type = "h", pch = 16,
    highlight.3d = TRUE, box = TRUE,
    main="shape", xlab = "x", ylab = "y", zlab = "z")

octa <- read.csv("data/xyz6.data", header = TRUE, sep = " ")
png(filename = "img/octa.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(octa$x, octa$y, octa$z,
    type = "h", pch = 16,
    highlight.3d = TRUE, box = TRUE,
    main="shape", xlab = "x", ylab = "y", zlab = "z")

cube <- read.csv("data/xyz8.data", header = TRUE, sep = " ")
png(filename = "img/cube.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(cube$x, cube$y, cube$z,
    type = "h", pch = 16,
    highlight.3d = TRUE, box = TRUE,
    main="shape", xlab = "x", ylab = "y", zlab = "z")

icosa <- read.csv("data/xyz12.data", header = TRUE, sep = " ")
png(filename = "img/icosa.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(icosa$x, icosa$y, icosa$z,
    type = "h", pch = 16,
    highlight.3d = TRUE, box = TRUE,
    main="shape", xlab = "x", ylab = "y", zlab = "z")

dodeca <- read.csv("data/xyz20.data", header = TRUE, sep = " ")
png(filename = "img/dodeca.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(dodeca$x, dodeca$y, dodeca$z,
    type = "h", pch = 16,
    highlight.3d = TRUE, box = TRUE,
    main="shape", xlab = "x", ylab = "y", zlab = "z")
