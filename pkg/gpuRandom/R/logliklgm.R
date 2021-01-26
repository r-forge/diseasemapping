#' @title pre_function for likfitGpu
#'
#' @useDynLib gpuRandom
#' @export


#before start, we have SpatialPointsDataFrame and spatial model

lgmGpuObjectes <- function(modelname, mydat, type=c("double", "float")){
        
        covariates = model.matrix(modelname$model$formula, data=modelname$data)
        
        temp = model.frame(modelname$model$formula, data=modelname$data)
        
        response=temp[,as.character(attributes(terms(temp))$variables)[2]]
        
        yX=vclMatrix(cbind(response,covariates),type=type)
        
        coordsGpu<-vclMatrix(mydat@coords,type=type)
        
        n = length(response)
        p = ncol(covariates)
        
        output<-list(yX=yX, coordsGpu=coordsGpu, n=n, p=p)
        
        output
}


#' @title Estimate Log-likelihood for Gaussian random fields
#'
#' @useDynLib gpuRandom
#' @export  

likfitGpu <- function( modelname, mydat, type=c("double", "float"),
                       paramsBatch, #vclMatrix of parameter sets,Vbatch, # matern correlation vclmatrix,diagMat, # D of cholesky decomposition
                       betas=NULL, #a vclmatrix  #given by the user or provided from formula
                       form = c("loglik", "ml", "mlFixSigma", "mlFixBeta", "reml", "remlPro"),
                       workgroupSize,
                       localSize,
                       NlocalCache,
                       verbose=FALSE){
        
        if((is.null(betas) & form == "loglik") | (is.null(betas) & form == "mlFixSigma") ) stop("need betas")
        
        yX <- lgmGpuObjectes(modelname, mydat, type)$yX
        
        coordsGpu <- lgmGpuObjectes(modelname, mydat, type)$coordsGpu   
        
        n <- lgmGpuObjectes(modelname, mydat, type)$n 
        p <- lgmGpuObjectes(modelname, mydat, type)$p  
        
        y <-  vclMatrix(yX[,1], nrow=n, ncol=1, type = type) 
        X <-  vclMatrix(yX[,c(2:(1+p))],type = type                                                                                    )
        
        form = c(loglik=1, ml=2, mlFixSigma=3, mlFixBeta=4, reml=5, remlPro=6)[form]
        
        rowbatch = nrow(paramsBatch)
        colbatch = 1     #colbatch = ncol(y)
        
        localSizechol<-localSize
        localSizechol[2]<-workgroupSize[2]
        
        
        Vbatch = vclMatrix(0, nrow(paramsBatch)*n, n, type = type)
        diagMat = vclMatrix(0, nrow(paramsBatch), n, type = type) 
        ab <- vclMatrix(0, nrow(Vbatch), colbatch+p, type = gpuR::typeof(Vbatch))
        temp2 <- vclMatrix(0, ncol(ab)*rowbatch, ncol(ab), type = gpuR::typeof(Vbatch))
        diagP <- vclMatrix(0, rowbatch, p, type = gpuR::typeof(Vbatch))
        temp3 <- vclMatrix(0, rowbatch*p, colbatch, type = gpuR::typeof(Vbatch))
        nine0 <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))
        aTDa <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))
        # ########################## 1, loglik or ml(beta,sigma), given beta and sigma #############################
        
        # Vbatch=LDL^T, cholesky decomposition
        gpuRandom:::maternBatchBackend(Vbatch, coordsGpu, paramsBatch,  workgroupSize, localSize)
        gpuRandom::cholBatch(Vbatch, diagMat, numbatchD=rowbatch, Nglobal=workgroupSize, Nlocal=localSizechol, NlocalCache=NlocalCache)


        logD <- apply(log(diagMat),1,sum)
        variances<-vclMatrix(paramsBatch[,3],nrow=rowbatch, ncol=1,type = type)
        
        if (form == 1 | form == 3) { # to get one, see the notes
                # temp = y-X*beta
                temp <- y - gpuRandom::gemmBatch(X, betas, 1L, 1L, colbatch, need_transpose = FALSE, workgroupSize)
                
                # L * A = temp, backsolve for A
                A <- vclMatrix(0, nrow(Vbatch), ncol(temp), type = gpuR::typeof(Vbatch))
                gpuRandom::backsolveBatch(A, Vbatch, temp, numbatchB=1L,  diagIsOne=TRUE,workgroupSize,  localSize,  NlocalCache)
                
                
                # one0 = A^T * D^(-1) * A = (y-X*betas)^T * V^(-1) * (y-X*betas)
                one0 <- vclMatrix(0, rowbatch*colbatch, colbatch, type = gpuR::typeof(Vbatch))
               
                #won't need if there is just one y batch (if colbatch=1)
                # one<- vclMatrix(0, nrow=rowbatch, ncol=colbatch, type = gpuR::typeof(Vbatch))
                
                gpuRandom:::crossprodBatchBackend(one0, A, diagMat,  invertD=TRUE,  workgroupSize, localSize, NlocalCache)
                
                #
                # for (j in 1:colbatch){
                #         for(i in 1:rowbatch){
                #                 one[i,j] <- one0[colbatch*(i-1)+j,j]
                #         }
                # }
                
        }else { # form == 2,4,5,6
                #profile = nlog hat_sigma^2 + log|D|
                # to get hat_sigma^2,  #L(a b) = (y X)
                gpuRandom::backsolveBatch(ab, Vbatch, yX, numbatchB=1L, diagIsOne=TRUE, Nglobal=workgroupSize, Nlocal=localSize, 
                                          NlocalCache)
                

                # temp2 = (ab)^T * D^(-1) *ab
                gpuRandom:::crossprodBatchBackend(temp2, ab, diagMat, invertD=TRUE, workgroupSize, localSize, NlocalCache)
                
                # b^T * D^(-1) * b = Q * P * Q^T, cholesky of a subset (right bottom) of temp2
                gpuRandom:::cholBatchBackend(temp2, diagP, c(colbatch, p, colbatch, p), c(0, rowbatch, 0, p), rowbatch, workgroupSize, localSizechol, NlocalCache) 
                
                # Q * temp3 = (b^T * D^(-1) *a), backsolve for temp3    2 by 1
                gpuRandom:::backsolveBatchBackend(temp3, temp2, temp2,
                                                  c(0,p,0,colbatch), c(colbatch, p, colbatch, p), c(colbatch, p, 0, colbatch),rowbatch,
                                                  diagIsOne=TRUE, workgroupSize, localSize, NlocalCache)
                

                # nine0 = temp3^T * P^(-1) * temp3,  four 1 by 1 matrices
                gpuRandom:::crossprodBatchBackend(nine0, temp3, diagP, invertD=TRUE,  workgroupSize, localSize, NlocalCache) ##doesn't need selecting row/col
                
                # won't need if there is just one y batch (if colbatch=1)
                #nine <- vclMatrix(0, colbatch*rowbatch, colbatch, type = gpuR::typeof(Vbatch))  
                
                #extract a^TDa from temp2
                for (j in 1:colbatch){
                        for (i in 1:rowbatch){
                                aTDa[i,j]= temp2[(i-1)*ncol(ab)+j, j]
                        }
                }
                #extract needed cells from nine0 # won't need if there is just one y batch (if colbatch=1)
                # for (j in 1:colbatch){
                #         for (i in 1:rowbatch){
                #                 nine[i,j]= nine0[(i-1)*colbatch+j, j]
                #         }
                # }
                two = aTDa - nine0
        }
        
        if (form == 1 | form == 4){
                #part1 = n*log(sigma^2)+log |D|
                part1 <- n*log(paramsBatch[,3]) + logD
                #replicate the part1 to do the plus operation
                #part1 <- matrix(part1, nrow=length(part1), ncol=colbatch, byrow=F)
                part1 <- vclMatrix(part1,nrow=length(part1), ncol=1, type = type)
        }else if (form == 5 | form==6){
                logP <- apply(log(diagP),1,sum)
        }
        
        
       
        
        
        if (form == 1 ) { #loglik
                # n*log(sigma^2) + log |D| + one/variances
                result <- part1 + one0/variances + n*log(2*pi)
        } else if (form == 2) {#ml
                result = n*log(two) + logD + n*log(2*pi) + n
        } else if (form == 3){ # mlFixSigma/ or ml(beta,hatsigma)
                result = n*log(one0/n)+logD + n*log(2*pi) + n
        } else if (form == 4){ # mlFixBeta / or ml(hatbeta,sigma)
                result = part1 + two/variances + n*log(2*pi) 
        } else if (form == 5){ #reml
                first_part <- (n-p)*log(variances) + logD + logP
                result <- first_part + two/variances + n*log(2*pi) 
        } else if (form == 6) { #remlPro  #(n-p)*log two + log|D| + log|P|
                result <- (n-p)*log(two/(n-p)) + logD+logP + n*log(2*pi) + n-p
        }
        
        result
}


     
     
     
     
     
     
     
     
     
     
     