

devtools::install_github("wrathematics/thrust")
devtools::install_github("gpuRcore/gpuRcuda")


ORDER = 100
A <- matrix(rnorm(ORDER^2), nrow=ORDER)
B <- matrix(rnorm(ORDER^2), nrow=ORDER)

library("gpuRcuda")
gpuA <- cudaMatrix(A, type="double")
gpuB <- cudaMatrix(B, type="double")

gpuC <- gpuA %*% gpuB

out = matrix(0, nrow(gpuA), ncol(gpuA))

gpuRcuda:::cudaMatToSEXP(gpuA@address, out, 8L)

out[1,1]
A[1,1]


# create a matrix
