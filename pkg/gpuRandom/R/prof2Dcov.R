#' @title 2D profile Log-likelihood for covariance parameters and lambda
#' @import ggplot2
#' @import data.table
#' @useDynLib gpuRandom
#' @export





        prof2Dcov <- function(LogLik,  # cpu matrix
                              paramToEstimate = c("range","nugget"),
                              cilevel,  # decimal
                              ParamList, 
                              params, # cpu matrix, 
                              Ndata,
                              Nobs,
                              Ncov,
                              verbose=FALSE){
            
          
          if (length(paramToEstimate) > 2){
            stop('should fix the rest params for a 2D plot')
          }
          maximum <- max(LogLik)
          breaks = maximum - qchisq(cilevel,  df = 2)/2
         
          result = as.data.table(cbind(LogLik, params[ ,paramToEstimate]))


          profileLogLik2D <- result[, .(profile=max(.SD)), by=paramToEstimate]
          profileLogLik2D <- as.matrix(profileLogLik2D)
          
          lMatrix = matrix(profileLogLik2D[,3], length(ParamList[[1]]), length(ParamList[[2]]))
          contour(ParamList[[1]]/100, ParamList[[2]], lMatrix,
                  col = par("fg"), lty = par("lty"), lwd = par("lwd"),
                  add = FALSE, levels = breaks)
          points(x = params[,paramToEstimate[1]]/100, y = params[,paramToEstimate[2]], pch=20, col='grey', cex=0.5)
          
        }
        
        
        
    
        
        
        
        
        
        
        # index <- which(LogLik == max(LogLik, na.rm = TRUE), arr.ind = TRUE)
        # 
        # makedata <- data.frame(cbind(params[,paramToEstimate], LogLik[,index[2]]))
        # colnames(makedata) <- c(paramToEstimate, 'LogLik')
        # 
        # lowx = min(params[,paramToEstimate[1]])
        # uppx = max(params[,paramToEstimate[1]])
        # lowy = min(params[,paramToEstimate[2]])
        # uppy = max(params[,paramToEstimate[2]])
        # 
        # 
        # contourplot <- ggplot(makedata, aes_string(x = paramToEstimate[1], y = paramToEstimate[2], z='LogLik')) + theme_bw() +
        #   scale_y_continuous(limits=c(lowy, uppy))+
        #   scale_x_continuous(limits=c(lowx, uppx) )+
        #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        #   geom_point(size = 0.5, colour='grey') +
        #   stat_contour(breaks=breaks) 
        # 
        # 
        # 
        # 
        # contourplot        