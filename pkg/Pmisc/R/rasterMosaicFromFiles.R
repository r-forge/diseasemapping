#' Raster mosaic from file names
#' 
#' @description Mosaic raster files into a new raster with a larger spatial extent
#'
#' @param x a string specifying raster files
#' @param cropOut area to crop the output raster to
#' @param fun Function, see \code{\link[raster]{mosaic}}
#' @param tolerance, numeric, see \code{\link[raster]{mosaic}}
#' @param filename Output file, see \code{\link[raster]{writeRaster}}
#' @param overwrite Logical, see \code{\link[raster]{writeRaster}}
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#' @return A \code{\link[raster]{raster}}
#' @seealso \code{\link[raster]{mosaic}}, \code{\link[raster]{writeRaster}}
#' @export 
mosaicFromFiles = function(x, cropOut=NULL, fun='mean',
  tolerance = 0.05,
  filename = paste(tempfile(), '.grd', sep=''),
  overwrite =  file.exists(filename),
  ...
) {
  
  if(!is.null(cropOut)) {
    raster1 = raster::raster(x[1])
    
    cropOut = raster::extend(
      raster::extent(
        raster::projectExtent(
          cropOut, 
          raster1)
      ), 
      raster::res(raster1)/10
    )
    
    rasterList = mapply(
      function(Dx) {
        res = try(raster::crop(raster::raster(Dx), cropOut))      
        if(class(res) =='try-error')
          res = NULL
        res
      },
      Dx=x)
  } else {
    rasterList = mapply(raster::raster, x)
  }
  
  rasterList = rasterList[!unlist(lapply(rasterList, is.null))]
  
  names(rasterList) = c('x','y', 
    paste(rep('r', length(rasterList)-2), 
      seq(from=3, 
        len=pmax(0,length(rasterList)-2), 
        by=1), 
      sep=''))
  
  rasterList$fun = fun
  rasterList$tolerance = tolerance
  
  res = do.call(
    raster::mosaic,
    rasterList  
  )
  
  writeRaster(res, 
    filename = filename, ..., 
    overwrite =  overwrite)
}
