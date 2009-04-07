### RCSId = "$Id: inla.R,v 1.44 2009/04/01 10:17:06 hrue Exp $"

`inla` =
    function (formula,
              family = "gaussian", 
              data = list(),
              quantiles=c(0.025,0.975),
              E = NULL,
              offset=NULL,
              scale = NULL,
              Ntrials = NULL,
              verbose = FALSE,
              control.compute = list(),
              control.predictor = list(),
              control.data = list(),
              control.inla = list(),
              control.results = list(),
              control.fixed = list(),
              control.mode = list(),
              inla.call = ifelse(.Platform$OS.type == "windows", "inla.exe", "inla"),
              num.threads = NULL,
              keep = FALSE,
              working.directory = NULL,
              only.hyperparam = FALSE,
              disable.check = FALSE,    # internal
              debug = FALSE,
              midFunction=NULL)            # internal
{
    if (!disable.check)
        inla.check.inla.call(inla.call)

    gp = inla.interpret.formula(formula)
    call = deparse(match.call())
    mf = match.call(expand.dots = FALSE)
    mf = mf[names(mf)!="midFunction"]
    
    if (gp$n.fix > 0)
        gp$model.matrix = model.matrix(gp$fixf,data=data)
    else
        gp$model.matrix = NULL
    
    ##control what should be computed
    cont.compute = inla.set.control.compute.default()
    cont.compute[(namc = names(control.compute))] = control.compute

    if(only.hyperparam)
    {
        cont.compute$hyperpar = TRUE
        cont.compute$dic = FALSE
        cont.compute$mlik = FALSE
        cont.compute$cpo = FALSE 
    } 
    
    ##control predictor section
    cont.pred = inla.set.control.predictor.default()
    cont.pred[(namc = names(control.predictor))] = control.predictor
    
    if(cont.compute$cpo || cont.compute$dic) 
        cont.pred$compute=TRUE
    if(!is.null(cont.compute$param) && length(cont.compute$param)!=2)
        stop("Two parameters for the gamma prior for the precision of predictor have to be provided")
    if(!is.null(cont.compute$initial) && length(cont.compute$initial)!=1)
        stop("One initial value for the precision of predictor have to be provided")
    if(only.hyperparam)
        cont.pred$compute = cont.pred$return.marginals = cont.pred$cdf = FALSE
        
    ##
    cont.inla =inla.set.control.inla.default(family)
    cont.inla[(namc = names(control.inla))] = control.inla

    ##control fixed
    cont.fixed = inla.set.control.fixed.default()
    cont.fixed[(namc = names(control.fixed))] = control.fixed

    ##control data section
    cont.data = inla.set.control.data.default()
    cont.data[(namc = names(control.data))] = control.data

    if(family!="T" && !is.null(cont.data$dof.max)) 
        stop("The parameter dof is defined only for family='T'")

    ##control results
    cont.res = inla.set.control.results.default()
    cont.res[(namc = names(control.results))] = control.results
    
    all.labels = character(0)
    ##if the data.dit or results.dir is specified then we keep both
    ##data and resulta
    if(!is.null(working.directory))
        keep=TRUE

    ##Create the directory where to store Model.ini and data.files and
    ##results.file
    if(keep)
    {
        ##create the directory locally or whereever specified by the user
        if(is.null(working.directory))
        {
            working.directory="inla.model"
            working.directory.start="inla.model"
        }
        else
            working.directory.start= working.directory
        ##if already exists then create one more 
        ans=file.exists(working.directory)
        kk=1
        while(ans)
        {
            working.directory=paste(working.directory.start,"-", kk, sep="")
            kk=kk+1
            ans=file.exists(working.directory)
        }
        inla.dir=working.directory
        xx=dir.create(inla.dir,showWarnings=FALSE)
        if(!xx) 
            stop(paste("You have no permission to create the directory ",inla.dir))
        else
            cat("Created working directory: ",inla.dir,"\n")
    }
    else
    {
        ##create a temporary directory
        inla.dir=tempfile()
        inla.dir=gsub("\\\\",.Platform$file.sep,inla.dir)
        dir.create(inla.dir)
    }
    ## Create a directory where to store the data files.....
    data.dir=paste(inla.dir,.Platform$file.sep,"data.files",sep="")
    dir.create(data.dir)
    ## ...and one to store the results
    results.dir = paste(inla.dir,.Platform$file.sep,"results.files",sep="")
    
    ## create the .ini file and make the problem.section
    file = paste(inla.dir,.Platform$file.sep,"Model.ini",sep="")

    if(debug) 
        print("prepare problem section")
    inla.problem.section(file,data.dir=data.dir,result.dir=results.dir,
                         hyperpar=cont.compute$hyperpar,
                         dic=cont.compute$dic, mlik=cont.compute$mlik,cpo=cont.compute$cpo,
                         quantiles=quantiles, smtp=cont.compute$smtp, q=cont.compute$q)
    
    ## PREPARE RESPONSE AND FIXED EFFECTS
    if (debug)
        cat("Prepare inla file.....")
    if(gp$n.fix!=0)
        mf$formula =gp$fixf
    else if(gp$n.random!=0)
        mf$formula =gp$randf
    else
        stop("Some covariate has to be present in the model")
    
    mf$na.action = na.pass
    mf$family = mf$control.predictor = mf$control.compute = mf$quantiles = mf$control.results = NULL
    mf$control.data = mf$control.inla = mf$control.fixed = mf$control.mode = mf$verbose = NULL
    mf$inla.call = mf$only.hyperparam = mf$keep = mf$working.directory =  mf$debug = mf$disable.check = NULL
    mf$num.threads = NULL
    mf$drop.unused.levels = TRUE
    mf[[1]] = as.name("model.frame")

    rf = mf
    rf$scale = rf$Ntrials = rf$offset = rf$E =  NULL

    ff = mf
    ff$scale = ff$Ntrials = ff$offset = ff$E = NULL

    wf = mf
    wf$scale = wf$Ntrials = wf$offset = wf$E = NULL
    
    mf = eval(mf, parent.frame())
    tot.data=nrow(mf)
    ind=seq(0,tot.data-1)

    scale = model.extract(mf, "scale")
    Ntrials = model.extract(mf, "Ntrials")
    E = model.extract(mf, "E")

    ## this takes care of the offset: `offset' is the argument, `offset.formula' is in the formula and `offset.sum'
    ## is their sum
    offset = as.vector(model.extract(mf, "offset"))
    if (!is.null(gp$offset))
    {
        ## there can be more offsets
        offset.formula = 0
        for(i in 1:length(gp$offset))
            offset.formula = offset.formula + as.vector(eval(parse(text=gp$offset[i]),data))
    }
    else
        offset.formula = NULL

    if (!is.null(offset.formula) && is.null(offset))
        offset.sum = offset.formula
    else if (is.null(offset.formula) && !is.null(offset))
        offset.sum = offset
    else if (!is.null(offset.formula) && !is.null(offset.formula))
    {
        if (length(offset.formula) == length(offset))
            offset.sum = offset.formula + offset
        else
            stop("The offset defined in formula and in argument has different length.")
    }
    else
        offset.sum = NULL

    ## cat("offset.formula ", offset.formula, "\n")
    ## cat("offset         ", offset, "\n")
    ## cat("offset.sum     ", offset.sum, "\n")
    
    ## Create a file with the response
    file.data = inla.create.data.file(mf=mf, E=E, scale=scale, Ntrials=Ntrials, event=NULL, family=family,
        data.dir=data.dir, file=file, debug=debug)

    ## add a section to the ini file
    ## first check all the arguments in control.data
    known.likelihoodmodels.with.ntheta.0 = c("poisson", "binomial")
    known.likelihoodmodels.with.ntheta.1 = c("gaussian", "stochvol_t", "zeroinflated_poisson_0",
        "zeroinflated_poisson_1", "zeroinflated_binomial_0", "zeroinflated_binomial_1")
    known.likelihoodmodels.with.ntheta.2 = c("T", "stochvol_nig")
    known.likelihoodmodels.with.ntheta.3 = NULL

    if(!is.null(cont.data$fixed))
    {
        if (is.element(family, known.likelihoodmodels.with.ntheta.0) && length(cont.data$fixed)!=0)
            stop("The length of fixed in `control.data` has to be 0")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.1) && length(cont.data$fixed)!=1)
            stop("The length of fixed in `control.data` has to be 1")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.2) && length(cont.data$fixed)!=2)
            stop("The length of fixed in control.data has to be 2")
    }
    else
    {
        if (is.element(family, known.likelihoodmodels.with.ntheta.0))
            cont.data$fixed = NULL
        else if (is.element(family, known.likelihoodmodels.with.ntheta.1))
            cont.data$fixed = 0
        else if(is.element(family, known.likelihoodmodels.with.ntheta.2))
            cont.data$fixed = rep(0,2)
        else
            stop("Should not happen: 1*")
    }
    
    if(!is.null(cont.data$initial))
    {
        if (is.element(family, known.likelihoodmodels.with.ntheta.0) && length(cont.data$initial) != 0)
            stop("The length of initial in control.data has to be 0")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.1) && length(cont.data$initial) != 1)
            stop("The length of initial in control.data has to be 1")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.2) && length(cont.data$initial) != 2)
            stop("The length of initial in control.data has to be 2")
    }
    if(!is.null(cont.data$param))
    {
        if (is.element(family, known.likelihoodmodels.with.ntheta.0) && length(cont.data$param) != 0)
            stop("The length of param in control.data has to be 0")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.1) && length(cont.data$param) != 2)
            stop("The length of param in control.data has to be 2")
        else if (is.element(family, known.likelihoodmodels.with.ntheta.2) && length(cont.data$param) != 4)
            stop("The length of param in control.data has to be 4")
    }
    if(debug) 
    {
        print("prepare data section")
    }
    ##....then create the new section
    inla.data.section(file=file,family=family,file.data=file.data,control=cont.data)

    ##create the PREDICTOR section. if necessary create a file with
    ##the offset for all likelihood
    if(!is.null(offset.sum))
    {
        if(sum(is.na(offset.sum))>0) 
            stop("No NA values allowed in the offset vector!")
        os = cbind(ind,offset.sum)
        offset.file = tempfile(tmpdir=data.dir)
        file.create(offset.file)
        write(t(os),ncolumns=2,file=offset.file,append=FALSE)
        offset.file = gsub(data.dir, "$DATADIR", offset.file, fixed=TRUE)
    }
    else
        offset.file = NULL
    inla.predictor.section(file=file,n=tot.data,predictor.spec=cont.pred,offset.file=offset.file)

    if (!is.null(cont.pred$compute) && cont.pred$compute)
        all.labels = c(all.labels,"predictor")


    ##FIXED EFFECTS
    if(gp$n.fix>0)
    {
        nc = ncol(gp$model.matrix)
        labels = colnames(gp$model.matrix)
        for(i in 1:nc)
        {
            if (debug)
                cat("write label[", labels[i],"]\n")

            fixed.eff=cbind(ind,as.numeric(gp$model.matrix[,i]))
            ##remove lines with NA
            fixed.eff = fixed.eff[!is.na(fixed.eff[,2]),]
            file.fixed=tempfile(tmpdir=data.dir)
            file.create(file.fixed)
            write(t(fixed.eff),ncolumns=ncol(fixed.eff),file=file.fixed,append=FALSE)
            file.fixed = gsub(data.dir, "$DATADIR", file.fixed, fixed=TRUE)
            inla.linear.section(file=file,file.fixed=file.fixed,label=labels[i],
                                results.dir=paste("fixed.effect",inla.num(i),sep=""),
                                control=cont.fixed,only.hyperparam=only.hyperparam)
        }
    }

    ##RANDOM EFFECT OR LINEAR EFFCT DEFINED VIA f()
    nr=gp$n.random
    n.weights=0
    j=0
    extra.fixed=0
    
    if(nr>0)
    {
        rf$formula =gp$randf
        rf = eval(rf, parent.frame())    

        if(gp$n.weights>0)
        {
            wf$formula = gp$weightf
            wf = eval(wf, parent.frame())
        }
        else 
            wf = NULL
        
        name.random.dir=c()
        if(nr!=(ncol(rf)-1)) stop("SOMETHING STRANGE")
               
        location = list()
        covariate = list()

        count.linear = 0
        count.random = 0
        
        for(i in 1:nr)
        {
            if(gp$random.spec[[i]]$model != "linear" && gp$random.spec[[i]]$model != "z")
            {
                ##in this case we have to add a FFIELD section.........
                count.random = count.random+1
                xx=rf[,i+1]
                if(is.factor(xx))
                {
                    location[[i]] = sort(unique(xx))
                    cov = match(xx,location[[i]])-1
                    cov[is.na(cov)] = -1
                    covariate[[i]] = cov
                }
                else
                {
                    if(!is.null(gp$random.spec[[i]]$values))
                    {
                        location[[i]] = sort(gp$random.spec[[i]]$values)
                        cov = match(xx,location[[i]])-1
                        cov[is.na(cov)] = -1
                        covariate[[i]] = cov
                    }
                    else
                    {
                        location[[i]] = sort(unique(xx))
                        cov = match(xx,location[[i]])-1
                        cov[is.na(cov)] = -1
                        covariate[[i]] = cov
                    }
                }

                ##create a location and covariate file
                file.loc=tempfile(tmpdir=data.dir)
                file.create(file.loc)
                write(location[[i]],ncolumns=1,file=file.loc,append=FALSE)
                file.loc = gsub(data.dir, "$DATADIR", file.loc, fixed=TRUE)
                
                file.cov=tempfile(tmpdir=data.dir)
                file.create(file.cov)
                write(t(cbind(ind,covariate[[i]])),ncolumns=2,file=file.cov,append=FALSE)
                file.cov = gsub(data.dir, "$DATADIR", file.cov, fixed=TRUE)

                ##and in case a file for the extraconstraint
                if(!is.null(gp$random.spec[[i]]$extraconstr))
                {
                    A=gp$random.spec[[i]]$extraconstr$A
                    e=gp$random.spec[[i]]$extraconstr$e
                    
                    if(ncol(A)!=length(location[[i]])) 
                        stop("Ncol in matrix A(extraconstraind) does not correspont to the length of f")
                    file.extraconstr=tempfile(tmpdir=data.dir)
                    file.create(file.extraconstr)
                    write(c(as.vector(t(A)),e),ncolumns=1,file=file.extraconstr,append=FALSE)
                    file.extraconstr = gsub(data.dir, "$DATADIR", file.extraconstr, fixed=TRUE)
                }
                else
                    file.extraconstr = NULL
                
                ##....also if necessary a file for the weights
                if(!is.null(gp$random.spec[[i]]$weights))
                {
                    ##set possible NA to -1
                    www = wf[,n.weights+2]
                    if(sum(is.na(www))!=0)www[is.na(www)] = -1
                    
                    ##create a file for the weights
                    file.weights=tempfile(tmpdir=data.dir)
                    file.create(file.weights)
                    write(t(cbind(ind,www)),ncolumns=2,file=file.weights,append=FALSE)
                    file.weights = gsub(data.dir, "$DATADIR", file.weights, fixed=TRUE)

                    n.weights = n.weights+1
                }
                ##create a FFIELD section
                inla.ffield.section(file=file,file.loc=file.loc,file.cov=file.cov,file.extraconstr=file.extraconstr,
                                    file.weights=file.weights,n=length(location[[i]]),
                                    random.spec=gp$random.spec[[i]],
                                    results.dir=paste("random.effect",inla.num(count.random),sep=""),
                                    only.hyperparam=only.hyperparam, data.dir=data.dir)
            }
            else if(gp$random.spec[[i]]$model == "linear")
            {
                ##....while here we have to add a LINEAR section
                count.linear = count.linear+1
                xx=rf[,i+1]
                file.linear = tempfile(tmpdir=data.dir)
                file.create(file.linear)
                write(t(cbind(ind,xx)),ncolumns=2,file=file.linear,append=FALSE)
                file.linear = gsub(data.dir, "$DATADIR", file.linear, fixed=TRUE)

                cont = list(cdf=gp$random.spec[[i]]$cdf,
                    prec=gp$random.spec[[i]]$prec,mean=gp$random.spec[[i]]$mean)
                
                inla.linear.section(file=file,file.fixed=file.linear,label=gp$random.spec[[i]]$term,
                                    results.dir=paste("fixed.effect",inla.num(gp$n.fix+count.linear),sep=""),
                                    control=cont, only.hyperparam=only.hyperparam)
            }
            else
            {
                ## model == "z"
                if (dim(gp$random.spec[[i]]$Z)[1] != tot.data)
                    stop(paste("Number of data is", tot.data, "but dimension of Z is", dim(gp$random.spec[[i]]$Z)))

                inla.z.section(file=file, random.spec = gp$random.spec[[i]], data.dir = data.dir,
                               results.dir = results.dir, only.hyperparam = only.hyperparam)
            }
        }
    }
    
    if(n.weights != gp$n.weights) 
    {
        stop("strange with weights in the covariate!!")
    }

    ##CREATE INLA SECTION 
    inla.inla.section(file=file,inla.spec=cont.inla)

    ##create mode section
    ##control model
    cont.mode = inla.set.control.mode.default()
    cont.mode[(namc = names(control.mode))] = control.mode
    inla.mode.section(file=file,cont.mode)


    ######################ADD for BYM
    if(!is.null(midFunction) & class(midFunction)=="function"){

      struct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("besag",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
      unstruct=strsplit(gsub("^[[:space:]]*f[[:space:]]*\\([[:space:]]*", "", grep("iid",strsplit(as.character(formula)[3],"\\+")[[1]],value=T)), ",")[[1]][1]
      midFunction(inla.dir=inla.dir,n=length(data[,struct]),struct=struct,unstruct=unstruct)
    }




    if (debug)
    {
        cat("...done\n")
        cat("Run inla...")
    }
    num.t = ifelse(is.numeric(num.threads), paste(" -t ", num.threads, " ", sep=""), "")
    if(.Platform$OS.type == "unix")
    {
        if(is.null(inla.call))
            inla.call="inla"

        if (verbose==TRUE)
            echoc=system(paste(inla.call," -b -v ", num.t, file))
        else
            echoc=system(paste(inla.call," -b  ", num.t, file))
    }
    else if(.Platform$OS.type == "windows")
    {
        if(is.null(inla.call))
            inla.call="inla.exe"
        if (verbose==TRUE)
            echoc=try(system(paste(inla.call," -b -v ", num.t, file)),silent=TRUE)
        else
            echoc=try(system(paste(inla.call," -b  ", num.t, file)),silent=TRUE)
        ##this we need because on windows it may fail on 'exit'
        echoc = 0
    }
    else
        stop("This should not happen: 4*")

    if (debug)
        cat("..done\n")

    if(echoc==0)
    {
        ret = inla.collect.results(results.dir,control.results=cont.res, debug=debug,
            only.hyperparam=only.hyperparam)

        if(family=="stochvol_nig" && cont.compute$cpo==TRUE)
            print("WARNING: PIT values not computed for the stochvol.nig family!")

        ret$control.compute=cont.compute
        ret$control.predictor=cont.pred
        ret$control.inla=cont.inla
        ret$control.data=cont.data
        ret$call=call
        ret$family=family
        ret$data=data
        ret$offset=offset
        ret$Ntrials=Ntrials
        ret$E=E
        ret$formula=formula
        ret$control.fixed=control.fixed
        ret$inla.call = inla.call
        ret$num.threads = num.threads
        ret$disable.check = disable.check
        ret$model.matrix = gp$model.matrix
        class(ret) = "inla"
    }

    if(debug) cat("clean up\n")
    if(!keep) unlink(inla.dir,recursive=TRUE)
    ret$results.dir = results.dir
    ret$dir = inla.dir
    return (ret)
}
inla.check.inla.call = function(inla.call = ifelse(.Platform$OS.type == "windows", "inla.exe", "inla"))
{
    ##
    ## Signal an error if the inla-program cannot be started
    ##

    ok = TRUE
    if(.Platform$OS.type == "windows")
    {
        ret = try(system(paste(inla.call, " -ping"), intern = TRUE, ignore.stderr = TRUE,
            wait = TRUE, input = NULL, show.output.on.console = FALSE, minimized = TRUE,
            invisible = TRUE),silent = TRUE)
        if (class(ret) == "try-error")
            ok = FALSE
    }
    else
    {
        ret = try(system(paste(inla.call, " -ping"), intern = TRUE, ignore.stderr = TRUE,
            wait = TRUE, input = NULL),silent = TRUE)
        ok =  (length(grep("ALIVE", ret)) > 0 || length(grep("unknown option", ret)) > 0)
    }
    if (!ok)
        stop(paste("\n\n***ERROR*** Cannot find the inla-program[", inla.call, "].\n",
                   "            Please modify the [inla.call] argument\n\n", sep=""))

    return (ok)
}

        
