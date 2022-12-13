library("scatterplot3d") # load

data <- read.csv("xyz.data", header = TRUE, sep = " ")
png(filename = "shape.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(data$x, data$y, data$z,
    type = "h",
    pch = 16,
    highlight.3d = TRUE,
    box = TRUE,
    main="shape",
    xlab = "x",
    ylab = "y",
    zlab = "z")
