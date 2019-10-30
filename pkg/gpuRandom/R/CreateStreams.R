#' @title Create streams
#'
#' @description Create streams in R.
#' 
#' @param x A vector indicating the number of work-items.
#' @return A R-streams.
#' @useDynLib gpuRandom
#' @export

CreateStreams = function(
  numWorkItems) {
  
  
  cpp_mrg31k3pCreateStreams(numWorkItems) 
  

}



