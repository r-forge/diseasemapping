library("geostatsp")

data("loaloa")
loaloa = unwrap(loaloa)
elevationLoa = unwrap(elevationLoa)
eviLoa = unwrap(eviLoa)

pts = geostatData(
    y ~ elev + evi,
    data = loaloa,
    covariates = list(elev = elevationLoa, evi = eviLoa),
    grid = squareRaster(loaloa, cells = 20, buffer = 1e4)
)

stopifnot(identical(names(pts), c("data", "grid", "covariates")))
stopifnot(inherits(pts$data, "SpatVector"))
stopifnot(inherits(pts$grid, "SpatRaster"))
stopifnot(all(c("elev", "evi") %in% names(pts$data)))
stopifnot(all(c("elev", "evi") %in% names(pts$covariates)))

myCrs = crs("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
dataR = rast(matrix(1:100, 10, 10), extent = ext(0, 1000, 0, 1000), crs = myCrs)
names(dataR) = "y"
covR = rast(matrix(seq_len(100), 10, 10), extent = ext(0, 1000, 0, 1000), crs = myCrs)
names(covR) = "x"

ras = geostatData(
    y ~ x,
    data = dataR,
    covariates = covR
)

stopifnot(identical(names(ras), c("data", "grid", "covariates", "formula")))
stopifnot(inherits(ras$grid, "SpatRaster"))
stopifnot("x" %in% names(ras$data))
stopifnot("x" %in% names(ras$covariates))

ok = try(
    geostatData(y ~ x, data = data.frame(y = 1, x = 2)),
    silent = TRUE
)
stopifnot(inherits(ok, "try-error"))
