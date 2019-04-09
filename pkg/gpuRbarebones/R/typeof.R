#' @title Get gpuR object type
#' @description \code{typeof} determines the type (i.e. storage mode) of a 
#' gpuR object
#' @param x A gpuR object
#' @rdname typeof-gpuR-methods
#' @author Charles Determan Jr.
#' @rdname typeof-gpuR-methods
#' @export
setMethod('typeof', signature(x="vclMatrix"),
          function(x) {
              switch(class(x),
                     "ivclMatrix" = "integer",
                     "fvclMatrix" = "float",
                     "dvclMatrix" = "double",
                     "cvclMatrix" = "fcomplex",
                     "zvclMatrix" = "dcomplex",
                     "ivclMatrixBlock" = "integer",
                     "fvclMatrixBlock" = "float",
                     "dvclMatrixBlock" = "double",
                     stop("unrecognized class"))
          })


#' @rdname typeof-gpuR-methods
#' @export
setMethod('typeof', signature(x="vclVector"),
          function(x) {
              switch(class(x),
                     "ivclVector" = "integer",
                     "fvclVector" = "float",
                     "dvclVector" = "double",
                     "ivclVectorSlice" = "integer",
                     "fvclVectorSlice" = "float",
                     "dvclVectorSlice" = "double",
                     stop("unrecognized vclVector class"))
          })


