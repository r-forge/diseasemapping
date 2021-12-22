#' @title 1D profile Log-likelihood for covariance parameters and lambda
#' @import data.table
#' @importFrom rootSolve uniroot.all
#' @useDynLib gpuRandom
#' @export
  




   prof1dCov <- function(LogLik,  # cpu matrix
                         XVYXVX,  # cpu matrix
                         ssqResidual,  # cpu matrix
                         paramToEstimate = c("range","shape","nugget","anisoAngleDegrees", "anisoRatio", "boxcox"),
                         cilevel=0.95,  # decimal
                         params, # cpu matrix, 
                         boxcox,  # boxcox vallues, consistent with other functions
                         Ndata,
                         Nobs,
                         Ncov,
                         reml, 
                         verbose=FALSE){
  
         
     
     
          
         maximum <- max(LogLik)
         breaks = maximum - qchisq(cilevel,  df = 1)/2
         
         ############## output matrix ####################
         Table <- matrix(NA, nrow=length(union(paramToEstimate, 'boxcox'))+Ncov+1, ncol=3)
         rownames(Table) <-  c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = ''), "sdSpatial", union(paramToEstimate, 'boxcox'))
         colnames(Table) <-  c("estimate", paste(c('lower', 'upper'), cilevel*100, 'ci', sep = ''))
         
         
         # myplots <- vector("list", length(paramToEstimate))
         # names(myplots) <- paramToEstimate
         ############### profile for covariance parameters #####################
         if(is.element('range',paramToEstimate)){
           # get profile log-lik for range values
           result = data.table::as.data.table(cbind(LogLik, params[,"range"]))
           colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "range")
           profileLogLik <- result[, .(profile=max(.SD)), by=range]
           f1 <- approxfun(profileLogLik$range, profileLogLik$profile-breaks)  
           plot(profileLogLik$range, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='range')
           # plot(profileLogLik$range, profileLogLik$profile)
           curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
           abline(h =0, lty = 2)
           #myplots[['range']] <- plot.range
           # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(profileLogLik$range)
           upper = max(profileLogLik$range)
           MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
           # maxvalue <- rangeresults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci for range")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci for range")}
           }
           
           if(length(ci)==0 | length(ci)>2){
             warning("require a better param matrix")
             ci <- c(NA, NA)
           }
           Table["range",] <- c(MLE, ci)
           
         }
         
         
         
         if(is.element('shape',paramToEstimate)){
           result = as.data.table(cbind(LogLik, params[,"shape"]))
           colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "shape")
           profileLogLik <- result[, .(profile=max(.SD)), by=shape]
           plot(profileLogLik$shape, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab="shape")
           f1 <- approxfun(profileLogLik$shape, profileLogLik$profile-breaks)  
           #f1 <- splinefun(profileLogLik$shape, profileLogLik$profile-breaks, method = "monoH.FC")
           curve(f1, add = TRUE, col = 2) 
           abline(h =0, lty = 2)
           #myplots[['shape']] <- plot.shape
           # shaperesults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(profileLogLik$shape)
           upper = max(profileLogLik$shape)
           MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
           # maxvalue <- shaperesults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci for shape")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci for shape")}
           }
           
           if(length(ci)==0 | length(ci)>2){
             warning("require a better param matrix")
             ci <- c(NA, NA)
           }
           Table["shape",] <- c(MLE, ci)

         }       
         
         
         
         
         if(is.element('nugget',paramToEstimate)){
           result = as.data.table(cbind(LogLik, params[,"nugget"]))
           colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "nugget")
           profileLogLik <-result[, .(profile=max(.SD)), by=nugget]
           plot(profileLogLik$nugget, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='nugget')
           f1 <- approxfun(profileLogLik$nugget, profileLogLik$profile-breaks)  
           # f1 <- splinefun(profileLogLik$nugget, profileLogLik$profile-breaks, method = "monoH.FC")
           curve(f1, add = TRUE, col = 2, n=1001) 
           abline(h =0, lty = 2)
           #myplots[['nugget']] <- plot.nugget
           # nuggetresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(profileLogLik$nugget)
           upper = max(profileLogLik$nugget)
           MLE <-  optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
           # maxvalue <- nuggetresults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci for nugget")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci for nugget")}
           }
           
           if(length(ci)==0 | length(ci)>2){
             warning("require a better param matrix")
             ci <- c(NA, NA)
           }
           Table["nugget",] <- c(MLE, ci)

         }
         
         
         
         
         if(is.element('anisoRatio',paramToEstimate)){
           result = as.data.table(cbind(LogLik, params[,"anisoRatio"]))
           colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoRatio")
           profileLogLik <- result[, .(profile=max(.SD)), by=anisoRatio]
           plot(profileLogLik$anisoRatio, profileLogLik$profile-breaks,main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='anisoRatio')
           f1 <- approxfun(profileLogLik$anisoRatio, profileLogLik$profile-breaks)  
           # f1 <- splinefun(profileLogLik$anisoRatio, profileLogLik$profile-breaks, method = "monoH.FC")
           curve(f1, add = TRUE, col = 2, n=1001) 
           abline(h =0, lty = 2)
           #myplots[['anisoRatio']] <- plot.anisoRatio
           # anisoRatioresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(profileLogLik$anisoRatio)
           upper = max(profileLogLik$anisoRatio)
           MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
           # maxvalue <- anisoRatioresults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci for anisoRatio")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci for anisoRatio")}
           }
           
           if(length(ci)==0 | length(ci)>2){
             warning("require a better param matrix")
             ci <- c(NA, NA)
           }
           Table["anisoRatio",] <- c(MLE, ci)
           
         }
         
         
         
         if(is.element('anisoAngleRadians',paramToEstimate)){
            result = as.data.table(cbind(LogLik, params[,"anisoAngleRadians"]))
            colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleRadians")
            profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleRadians]
            plot(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='anisoAngleRadians')
            f1 <- approxfun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks)
            # f1 <- splinefun(profileLogLik$anisoAngleRadians, profileLogLik$profile-breaks, method = "monoH.FC")
            curve(f1, add = TRUE, col = 2, n=1001) 
            abline(h =0, lty = 2)
            #myplots[['anisoAngleRadians']] <- plot.Degrees
            # anisoAngleRadiansresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
            lower = min(profileLogLik$anisoAngleRadians)
            upper = max(profileLogLik$anisoAngleRadians)
            MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
            # maxvalue <- anisoAngleRadiansresults$objective
            # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
            ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
            if(length(ci)==1){
               if( ci > MLE){
                  ci <- c(lower, ci)
                  message("did not find lower ci")
               }else{
                  ci <- c(ci, upper)
                  message("did not find upper ci")}
            }
            
            if(length(ci)==0 | length(ci)>2){
               warning("require a better param matrix")
               ci <- c(NA, NA)
            }
            Table["anisoAngleRadians",] <- c(MLE, ci)
            
            
         }
         
         
         
         
         
         
         
         if(is.element('anisoAngleDegrees',paramToEstimate)){
           result = as.data.table(cbind(LogLik, params[,"anisoAngleDegrees"]))
           colnames(result) <- c(paste(c('boxcox'), round(boxcox, digits = 2) ,sep = ''), "anisoAngleDegrees")
           profileLogLik <- result[, .(profile=max(.SD)), by=anisoAngleDegrees]
           plot(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='anisoAngleDegrees')
           f1 <- approxfun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks)
           # f1 <- splinefun(profileLogLik$anisoAngleDegrees, profileLogLik$profile-breaks, method = "monoH.FC")
           curve(f1, add = TRUE, col = 2, n=1001) 
           abline(h =0, lty = 2)
           #myplots[['anisoAngleDegrees']] <- plot.Degrees
           # anisoAngleDegreesresults <- optim(0.1, f1, method = "L-BFGS-B",lower = 0.1, upper = 1.5, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(profileLogLik$anisoAngleDegrees)
           upper = max(profileLogLik$anisoAngleDegrees)
           MLE <- round(optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum,digits=4)
           # maxvalue <- anisoAngleDegreesresults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci")}
           }
           
           if(length(ci)==0 | length(ci)>2){
             warning("require a better param matrix")
             ci <- c(NA, NA)
           }
           Table["anisoAngleDegrees",] <- c(MLE, ci)

           
         }
         
         
         
         
         index <- which(LogLik == max(LogLik, na.rm = TRUE), arr.ind = TRUE)
         ###############lambda hat#####################
         if(is.element('boxcox',paramToEstimate)  & length(boxcox)>5 ){
           likForboxcox = cbind(boxcox, apply(LogLik, 2,  max) )
           f1 <- approxfun(likForboxcox[,1], likForboxcox[,2]-breaks)
           # f1 <- splinefun(likForboxcox[,1], likForboxcox[,2]-breaks, method = "monoH.FC")
           plot(likForboxcox[,1], likForboxcox[,2]-breaks, main="Profile LogL, y axis adjusted", ylab= "proLogL-breaks", xlab='boxcox')
           curve(f1(x), add = TRUE, col = 2, n = 1001)   #the number of x values at which to evaluate
           abline(h =0, lty = 2)
           #myplots[['boxcox']] <- plot.boxcox
           # rangeresults <- optim(26000, f1, method = "L-BFGS-B",lower = 20000, upper = 240000, hessian = FALSE, control=list(fnscale=-1) )
           lower = min(boxcox)
           upper = max(boxcox)
           MLE <- optimize(f1, c(lower, upper), maximum = TRUE, tol = 0.0001)$maximum
           # maxvalue <- rangeresults$objective
           # breaks = maxvalue - qchisq(cilevel,  df = 1)/2
           ci<-rootSolve::uniroot.all(f1, lower = lower, upper = upper)
           if(length(ci)==1){
             if( ci > MLE){
               ci <- c(lower, ci)
               message("did not find lower ci for boxcox")
             }else{
               ci <- c(ci, upper)
               message("did not find upper ci for boxcox")}
           }else if(length(ci)>2){
             warning("error in param matrix")
             ci <- c(NA, NA)
           }
           Table["boxcox",] <- c(MLE, ci)
         }else if(is.element('boxcox',paramToEstimate)  & length(boxcox) <= 5){
           message("boxcox: not enough values for interpolation!")
         }else{
           Table["boxcox",1] <- boxcox[index[2]]
         }
         
         NcolTotal = Ndata + Ncov
         
         ###############betahat#####################
         Betahat <- matrix(0, nrow=Ncov, ncol=Ndata)
         a<-c( ((index[1]-1)*Ncov+1) : (index[1]*Ncov) )
         mat <- XVYXVX[a,((Ndata+1):NcolTotal)]
         mat[upper.tri(mat)] <- mat[lower.tri(mat)]
         Betahat <- solve(mat) %*% XVYXVX[a,index[2]]
         
         Table[c("intercept", paste(c('betahat'), seq_len(Ncov-1),sep = '')),1] <- Betahat
         
         
         
         
         #################sigma hat#########################
         if(reml==FALSE)  {
           Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/Nobs)
         }else{         
           Table["sdSpatial",1] <- sqrt(ssqResidual[index[1],index[2]]/(Nobs - Ncov))
         }
         
         
         
         Output <- list(summary = Table,
                        breaks = breaks
                        )
         
         
         
         Output
         
         
         
         
   }      
         
         
        