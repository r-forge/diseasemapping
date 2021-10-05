#' @title Create streams on GPU
#'
#' @description Create streams in R.
#' 
#' @param n number of streams to create.
#' @param seed the initial seed of streams.
#' @return A stream object on GPU.
#' @useDynLib gpuRandom
#' @export




createStreamsGpu = function(n, 
                            seed=12345){

    seed = rep_len(seed, 6)
    seedVec <- gpuR::vclVector(as.integer(seed), type="integer")  
    streamsR<-gpuR::vclMatrix(0L, nrow=as.integer(n), ncol=12, type="integer")

    gpuRandom:::CreateStreamsGpuBackend(seedVec, streamsR, keepInitial=1)
  
    streamsR
  
}










