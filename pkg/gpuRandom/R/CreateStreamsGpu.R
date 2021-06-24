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




CreateStreamsGpu = function(seedR, Nstreams, keepinitial=1) {
  
  streamsR<-vclMatrix(0L, nrow=Nstreams, ncol=18, type="double")
  
  gpuR::colnames(streamsR) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
                               "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
                               "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")
  
  
 # if(typeof(seedR)=="double"){
    
    myseedR = packBits( intToBits(as.vector(rbind(0L, seedR))), type = 'double')
    
    myseed <- gpuR::vclVector(myseedR, type="double")
 # }
  
  
  gpuRandom:::CreateStreamsGpuBackend(myseed, streamsR, keepinitial)
  
  
  streamsR
  
  
  
  
}

