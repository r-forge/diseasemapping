#' @title Create streams on GPU
#'
#' @description Create streams in R.
#' 
#' @param seedR an R vector, which is the initial seed of streams.
#' @param Nstreams number of streams to create.
#' @param keepinitial logical, whether to keep the initial seed in the created stream.
#' @return A stream object on GPU.
#' @useDynLib gpuRandom
#' @export




CreateStreamsGpu = function(seedR=c(12345, 12345, 12345, 12345, 12345, 12345), 
                            Nstreams, 
                            keepInitial=1) {


  
   # if(typeof(seedR)=="integer"){
   #  
   #  myseedR = packBits( intToBits(as.vector(rbind(0L, seedR))), type = 'double')
   #  
   #  }else if (typeof(seedR)=="double"){
   #  myseedR = packBits( numToBits(seedR), type = 'double')  
   #  }
  
  
    myseed <- gpuR::vclVector(as.integer(seedR), type="integer")  
    streamsR<-vclMatrix(0L, nrow=Nstreams, ncol=18, type="integer")
    
    gpuRandom:::CreateStreamsGpuBackend(myseed, streamsR, keepInitial)
  
  
    streamsR
  
  

}










