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




createStreams = function(n, 
                        seed=12345) {

    seed = rep_len(seed, 6)
    seedVec <- gpuR::vclVector(as.integer(seed), type="integer")  
    streamsR<-gpuR::vclMatrix(0L, nrow=as.integer(n), ncol=12, type="integer")

    keepInitial=1
    gpuRandom:::CreateStreamsGpuBackend(seedVec, streamsR, keepInitial)
  
    streamsR
  
}










