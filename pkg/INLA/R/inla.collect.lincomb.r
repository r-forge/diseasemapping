`inla.collect.lincomb` =
    function(results.dir,
             return.marginals.random=T,
             effect.string="lincomb",
             debug = FALSE)
{
    alldir = dir(results.dir)
    random = alldir[grep(paste("^",effect.string, sep="") ,alldir)]
    n.random = length(random)
    if(debug) print("collect random effects")

    ##read the names and model of the random effects
    if(n.random>0)
    {
        names.random = inla.namefix(character(n.random))
        model.random = inla.trim(character(n.random))
        for(i in 1:n.random)
        {
            tag = paste(results.dir,.Platform$file.sep,random[i],.Platform$file.sep,"TAG",sep="")
            names.random[i] = inla.namefix(readLines(tag,n=1))
            #modelname = inla.trim(paste(results.dir,.Platform$file.sep,random[i],.Platform$file.sep,"MODEL",sep=""))
            #model.random[i] = inla.trim(readLines(modelname,n=1))
        }

        summary.random = list()
        marginals.random = list()

        for(i in 1:n.random)
        {
            ##read the summary
            file= paste(results.dir,.Platform$file.sep,random[i],sep="")
            dir.random = dir(file)

            dd = matrix(inla.read.binary.file(file=paste(file,.Platform$file.sep,"summary.dat",sep="")),ncol=3,byrow=TRUE)
            col.nam = c("ID","mean","sd")

            ##read quantiles if existing
            if (debug) cat("...quantiles.dat if any\n")
            if(length(grep("^quantiles.dat$",dir.random))==1)
            {
                xx = inla.read.binary.file(paste(file,.Platform$file.sep,"quantiles.dat",sep=""))
                xx = t(inla.interpret.vector(xx))
                len = dim(xx)[1]
                qq = xx[seq(2,len,2),]
                if (is.null(dim(qq)))
                    qq = matrix(qq,1, length(qq))
                col.nam = c(col.nam,paste(as.character(xx[1,]),"quant",sep=""))
                dd = cbind(dd,qq)
            }

            ##read cdf if existing
            if (debug) cat("...cdf.dat if any\n")
            if(length(grep("^cdf.dat$",dir.random))==1)
            {
                xx = inla.read.binary.file(paste(file,.Platform$file.sep,"cdf.dat",sep=""))
                xx = t(inla.interpret.vector(xx))
                len = dim(xx)[1]
                qq = xx[seq(2,len,2),]
                if (is.null(dim(qq)))
                    qq = matrix(qq,1, length(qq))
                col.nam = c(col.nam,paste(as.character(xx[1,])," perc",sep=""))
                dd = cbind(dd,qq)
            }

            ##read kld
            if (debug) cat("...kld\n")
            kld1 = matrix(inla.read.binary.file(file=paste(file,.Platform$file.sep,"symmetric-kld.dat",sep="")),
                ncol=2,byrow=TRUE)

            qq = kld1[,2]
            if (is.null(dim(qq)))
                qq = matrix(qq,length(qq),1)
            dd = cbind(dd,qq)
            if (debug) cat("...kld done\n")

            col.nam = c(col.nam, "kld")
            colnames(dd) = inla.namefix(col.nam)
            summary.random[[i]] = as.data.frame(dd)

            if(return.marginals.random)
            {
                xx = inla.read.binary.file(paste(file,.Platform$file.sep,"marginal-densities.dat",sep=""))
                rr=inla.interpret.vector2(xx)
                rm(xx)
                nd = length(rr)
                names(rr) = inla.namefix(paste("index.", as.character(1:nd), sep=""))
                for(j in 1:nd)
                    colnames(rr[[j]]) = inla.namefix(c("x", "y"))
                marginals.random[[i]] = rr
            }
            else
                marginals.random=NULL
        }
        names(summary.random) = inla.namefix(names.random)
        if(!is.null(marginals.random))
            names(marginals.random) = inla.namefix(names.random)
    }
    else
    {
        if(debug) cat("No random effets")
        model.random=NULL
        summary.random=NULL
        marginals.random=NULL
    }

    res = list(summary.random=summary.random, marginals.random=marginals.random)
    return(res)
}
