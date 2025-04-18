
allVarsP = function(formula) {
  # return vector of variable names in the formula
  allterms = rownames(attributes(terms(formula))$factors)
  if(!length(allterms)) {
    # might be intercepty only but there may be an offset
    # add a junk variable so that the terms function produces a factors
    allterms =  setdiff(
      rownames(attributes(terms(
        update.formula(formula, .~ .+ junk)
      ))$factors),
      'junk')
  }
  
  firstTerm = as.character(formula)
  firstTerm = trimws(firstTerm[-c(1, length(firstTerm))])
  allterms = setdiff(allterms, firstTerm)
  
  # look for values=1:stuff in inla formula and replace with seq
  
  if(length(allterms)) {
    inlaValuesPattern = 'values[[:space:]]+?=[[:space:]]+?'
    haveInlaValues = grep(inlaValuesPattern, allterms, value=TRUE)
    notInlaValues = grep(inlaValuesPattern, allterms, 
                         value=TRUE, invert=TRUE)
    allterms = c(
      unique(unlist(strsplit(notInlaValues, ":"))),
      haveInlaValues
    )
  }
  allterms = gsub("[[:space:]]", "", allterms)
  # remove offset( or factor(
  alltermsPlain = gsub("^[[:alpha:]]+\\(|\\)$|[,].*", "", allterms)
  attributes(alltermsPlain)$orig = allterms
  alltermsPlain
}

gm.dataRaster = function(
    formula,
    data, grid=data,
    covariates=NULL,
    buffer=0){

  if(abs(diff(res(grid)))>0.000001 )
    warning("data is not on a square grid")
  
  cellsBoth = cellsBuffer(grid, buffer)			
  cellsSmall = cellsBoth$small
  
  # find factors
  
  allterms = allVarsP(formula)
  
  alltermsFull = attributes(allterms)$orig
  
  # find factors
  
  theFactors = grep("^factor", alltermsFull, value=TRUE)
  theFactors = gsub("^factor\\(|\\)$", "", theFactors)
  
  termsInF = grep("^f[(]", alltermsFull, value=TRUE)
  termsInF = gsub("^f[(]|[,].*|[[:space:]]", "", termsInF)
  
  covFactors = NULL
  for(D in names(covariates)) {
    if(any(is.factor(covariates[[D]])))
      covFactors = c(D, covFactors)
  }
  dataFactors = NULL
  for(D in names(data)) {
    if(is.factor(data[[D]]))
      dataFactors = c(D, dataFactors)
  }
  
  inModel = gsub("^[[:alnum:]]+[(]|[,].*|[[:space:]]",
                 "",  alltermsFull)
  inModel = gsub(
    "(,([[:alnum:]]|=|[[:space:]])+)?\\)$",#+?[[:space:]]?\\)[[:space:]]?$",
    "", inModel)
  
  theFormulaTerms= terms(formula)
  
  
  theOffset = rownames(attr(theFormulaTerms, 'factors'))[attr(theFormulaTerms, 'offset')]
  offsetToLogOrig = grep(
    "^offset\\([[:print:]]+,log=TRUE\\)$", 
    gsub("[[:space:]]+", "", theOffset), value=TRUE)
  

  
  if(length(offsetToLogOrig)) {
    names(offsetToLogOrig) = gsub(
      "^[[:space:]]?offset\\(|,[[:space:]]?log[[:space:]]?=[[:space:]]?TRUE[[:space:]]?\\)[[:space:]]?$",
      '', offsetToLogOrig
    )
  }
  
  offsetNotLogged = grep(
    "^offset\\([[:print:]]+,log=TRUE\\)$", 
    gsub("[[:space:]]+", "", theOffset), 
    invert=TRUE, value=TRUE)
  
  Sfactor = c(
    dataFactors,
    covFactors,
    theFactors
  )
  covFactors = intersect(Sfactor,names(covariates))
  dataFactors = intersect(Sfactor,names(data))
  
  inModel = intersect(inModel, names(covariates))
  
  
  if(length(inModel)) {
    
    if(length(grep("SpatRaster", class(covariates)))) {
      covariates = covariates[[inModel]]
    } else {
      covariates = covariates[inModel]
    }
    
    dataFactors = intersect(Sfactor, names(data))
    
    notInData = setdiff(names(covariates), names(data))
    
    
    rmethod = rep("bilinear", length(names(covariates)))
    names(rmethod) = names(covariates)
    rmethod[covFactors] = "near"
    
    
    notLogOffset = ! names(covariates) %in% names(offsetToLogOrig)
    
    if(any(notLogOffset)){
      # have some cvariates, which are not offsets
      if(length(grep("SpatRaster", class(covariates)))) {
        
        covariatesForStack = covariates[[which(notLogOffset)]]
        
        covariatesForStackData = covariates[[notInData]]
        
      } else {
        # covariates excluding offsets, for prediction
        covariatesForStack = covariates[notLogOffset]
        # covariates including offsets, for model
        covariatesForStackData = covariates[notInData]
      }
      
      # covariate raster at resolution of model (cellsSmall, no buffer)
      covariatesStack = stackRasterList(
        covariatesForStack,
        cellsSmall, method=rmethod)
      
      covariatesStack = c(cellsSmall, covariatesStack)
      
      if(length(covariatesForStackData)) {
        covData = stackRasterList(
          covariatesForStackData, 
          data, method=rmethod)
      } else {
        covData = NULL
      }
      
      
    } else { # only log offset
      covariatesStack = cellsSmall # for predictions
      covData = NULL # for model
    }
    
    covData = covData[[setdiff(names(covData), names(offsetToLogOrig))]]
    for(D in names(offsetToLogOrig)) {
      # loop through offsets which should
      # be aggregated before taking logs
      offsetToLog = covariates[[D]]
      
      toCrop = union(
        project(
          ext(covariatesStack),
          crs(covariatesStack), 
          crs(offsetToLog)
        ),
        project(ext(data), crs(data),
                crs(offsetToLog)
        )
      )
      
      
      offsetToLogCrop = crop(
        offsetToLog, 
        toCrop
      )
      
      offsetToLogCrop = project(
        offsetToLogCrop,
        y=crs(covariatesStack),
        method='near')
      
      toAggregatePredictions = floor(min(res(covariatesStack)/res(offsetToLogCrop)))
      toAggregateData = floor(min(res(data)/res(offsetToLogCrop)))
      
      # aggregate for predictions
      if(any(toAggregatePredictions > 1)){
        offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregatePredictions, 
                                   fun=sum, na.rm=TRUE)
      } else {
        toAggregatePredictions = 1
        offsetToLogAgg = offsetToLogCrop
      }
      offsetToLogAgg = project(offsetToLogAgg, covariatesStack)
      
      offsetToLogAgg = classify(
        offsetToLogAgg, 
        t(c(-Inf,0,NA)) 
      )
      
      offsetToLogLogged = log(offsetToLogAgg) - sum(log(rep_len(toAggregatePredictions,2)))
      names(offsetToLogLogged) = paste('log',D,sep='')
      covariatesStack = c(covariatesStack, offsetToLogLogged)
      
      # aggregate differently for model fitting
      if(toAggregateData > 1 ){
          if(toAggregateData != toAggregatePredictions) {
            offsetToLogAgg = aggregate(offsetToLogCrop, fact=toAggregateData, fun=sum)
          }
      } else {
        offsetToLogAgg = offsetToLogCrop
      }
      offsetToLogAgg = project(offsetToLogAgg, covData)
      offsetToLogAgg = classify(
        offsetToLogAgg, 
        t(c(-Inf,0,NA)) 
      )
      
      offsetToLogLogged = log(offsetToLogAgg) + 
          sum(log(res(covariatesStack))) -
          sum(log(res(offsetToLogCrop)))
      names(offsetToLogLogged) = paste('log',D,sep='')
      covData = c(covData, offsetToLogLogged)

    } # loop D offsets to log
    
    if(length(offsetToLogOrig)) {
      theNewOffset = paste0("offset(log", names(offsetToLogOrig), ")")
      
      formulaOrig = formula
      formula = stats::reformulate(
        c( attr(theFormulaTerms, 'term.labels'), offsetNotLogged, theNewOffset), 
        response=rownames(attr(theFormulaTerms, 'factors'))[attr(theFormulaTerms, 'response')], 
        intercept=attr(theFormulaTerms, 'intercept'), 
        env=environment(formula))
    }
    

    covariatesSP = as.points(covariatesStack)
    covariatesDF = values(covariatesSP)
    
    data = c(data, covData)			
    
    
  } else { # if length(inmodel)
    # allcovariates are in data
    covariatesDF = data.frame()
  }
  
  
  
  if(any(res(data)>1.25*res(cellsSmall)))
    warning("data is coarser than grid")
  
  data = c(data, resample(cellsSmall, data, method='near'))	
  
  
  dataSP = suppressWarnings(as.points(data))
  dataDF = values(dataSP)
  
  # get rid of rows with missing response if lgcp with count response
  
  if(names(dataDF)[1] == 'count')
    dataDF = dataDF[!is.na(dataDF$count), ]
  
  # redo factors
  # loop through spatial covariates which are factors
  for(D in intersect(Sfactor, names(covariatesDF))) {
    theTable = sort(table(dataDF[[D]]), decreasing=TRUE)
    theLevels = levels(covariates[[D]])[[1]]
    if(identical(theLevels, "")) {
      theLabels = paste("l", names(theTable),sep="")
    } else {
      theLabels = theLevels[
        match(as.integer(names(theTable)), theLevels$ID)
        ,"Category"]
    }
    dataDF[[D]] = factor(dataDF[[D]], levels=as.integer(names(theTable)),
                         labels=theLabels)			
    covariatesDF[[D]] = factor(covariatesDF[[D]], levels=as.integer(names(theTable)),
                               labels=theLabels)			
    
  }
  

  list(
    data=dataDF,
    grid=cellsSmall,
    covariates=covariatesDF,
    formula = formula
  )
}


#############
# data is a SpatialPointsDataFrame
#############

gm.dataSpatial = function(
    formula, data,  grid, 
    covariates, 
    buffer=0) {
  
  if(missing(covariates)) covariates = list()
  
  
  
  
  # check response variable is in data
  if(!all.vars(formula)[1] %in% names(data)){
    warning(paste(
      'response variable',
      all.vars(formula)[1],
      'not found in data'
    ))
  }
  
  alltermsPlain = allVarsP(formula)
  
  # find factors
  allterms = attributes(alltermsPlain)$orig
  # remove covariates not in the model
  keepCovariates = intersect(alltermsPlain, names(covariates))
  if(is.list(covariates)) {
    covariates = covariates[keepCovariates]
  } else {
    if(length(keepCovariates)) {
      covariates = covariates[[keepCovariates]]	
    } else {
      covariates = list()
    }
  }
  
  # check for missing CRS
  
  if(!nchar(crs(data)) | !nchar(crs(grid)) ) {
    if(nchar(crs(data))) {
      warning("assigning crs of grid to data")
      crs(grid) = crs(data)
    } else if(nchar(crs(grid))) {
      crs(grid) = crs(data)
      warning("assigning crs of data to grid")
    } else if(length(covariates)) {
      if(nchar(crs(covariates[[1]]))) {
        warning("assigning crs of first covariate to data")
        crs(grid) = crs(data) = crs(covariates[[1]])
      } else {
        warning("no crs supplied")
      }
    } else {
      warning("no crs supplied")
    }
  }
  
  theFactors = grep("^factor", allterms, value=T)
  theFactors = gsub("^factor\\(|\\)$", "", theFactors)
  
  
  termsInF = grep("^f[(]", allterms, value=TRUE)
  termsInF = gsub("^f[(]|[,].*|[[:space:]]", "", termsInF)
  
  
  covFactors = NULL
  for(D in names(covariates)) {
    if(any(is.factor(covariates[[D]])))
      covFactors = c(D, covFactors)
  }
  
  
  Sfactors = c(
    names(data)[unlist(lapply(values(data), is.factor))],
    covFactors,
    theFactors
  )
  Sfactors = unique(Sfactors)
  covFactors = intersect(Sfactors,names(covariates))
  
  cantFind = setdiff(Sfactors, c(names(data), names(covariates)))
  if(length(cantFind))
    warning("can't find variables", cantFind)
  
  # the grid
  cellsBoth = cellsBuffer(grid, buffer)
  cellsSmall = cellsBoth$small
  
  # 
  if(length(names(covariates))) {
    
    dataFactors = intersect(Sfactors, names(data))
    
    rmethod = rep("bilinear", length(names(covariates)))
    names(rmethod) = names(covariates)
    rmethod[covFactors] = "near"
    rmethod[intersect(names(covariates), termsInF)] = "near"
    
    covariatesStack = stackRasterList(
      covariates, 
      template=cellsSmall, 
      method=rmethod)
    covariatesStack2 = c(cellsSmall, covariatesStack)
    covariatesSP = suppressWarnings(as.points(covariatesStack2))
    covariatesDF = values(covariatesSP)
  } else { # else no covariates
    covariatesDF = data.frame()
  }
  
  # loop through covariates which aren't in data, extract it from `covariates`
  for(D in setdiff(alltermsPlain, names(data))){
    if(is.null(covariates[[D]]))
      warning("cant find covariate '", D, "' in covariates or data")
    
    if(any(class(covariates[[D]]) == 'SpatRaster')) {
      extractHere = terra::extract(covariates[[D]], 
                                   project(data, crs(covariates[[D]])), ID=FALSE,
                                   method = rmethod[D])
    } else {
      extractHere = terra::extract(covariates[[D]], 
                                   project(data, crs(covariates[[D]])))
      extractHere = extractHere[match(1:length(data), extractHere[,'id.y']), 2]
    }
    
    
    if(is.data.frame(extractHere)) {
      if(nrow(extractHere) != nrow(data)) {warning("mismatch in extracted covariates and data")}
      extractHere = extractHere[,grep("^ID$|^id.y", names(extractHere), invert=TRUE), drop=FALSE]
      terra::values(data)[[D]] = extractHere[,1]
    } else {
      data[[D]] = extractHere
    }
  } # D loop extracting covariates
  
  
  # reproject data to grid
  if(!identical(crs(cellsSmall), crs(data))) {
    data = project(data, crs(cellsSmall))
  }
  data$space = suppressWarnings(terra::extract(cellsSmall, data, 
                                               ID=FALSE, mat=FALSE, dataframe=FALSE, method = 'near'))
  
  # check for strange things where no cell value is found
  theNA = which(is.na(data$space))
  if(length(theNA)) {
    cellHere = cellFromRowCol(cellsSmall, 
                              rowFromY(cellsSmall, crds(data)[theNA, 2]) , 
                              colFromX(cellsSmall, crds(data)[theNA, 1]) 
    ) 
    data$space[theNA] = values(cellsSmall)[cellHere, 'space']
  }
  
  # loop through spatial covariates which are factors
  for(D in intersect(Sfactors, names(covariatesDF))) {
    theLevels = levels(covariates[[D]])[[1]]
    idCol = grep("^id$", names(theLevels), ignore.case=TRUE, value=TRUE)[1]
    if(is.na(idCol)) idCol = 1
    if(D %in% names(theLevels)) {
      labelCol = D
    } else {
      labelCol = grep("^category$|^label$", 
                      names(theLevels), ignore.case=TRUE, value=TRUE)[1]
    }
    if(is.na(labelCol)) labelCol = 2
    
    dataD = unlist(terra::values(data)[[D]])
    if(is.factor(dataD)) {
      
      # give covariatesDF factor levels from data
      if(all(levels(dataD) %in% theLevels[,labelCol])) {
        # match factor levels in data to 
        # factor levels in raster
        theTable = table(dataD)
        if(theTable[1]==0) {
          warning("no data in baseline level ", D)
          theTable = sort(theTable, decreasing=TRUE)
        }
        levelsHave = names(theTable)[theTable > 0]
        covariatesDF[[D]] = factor(
          as.character(covariatesDF[[D]]),
          levels = levelsHave
        )
        terra::values(data)[[D]] = factor(as.character(dataD), levels = levelsHave)
      } else { 
        # levels in data can't be found in raster levels
        # ignore raster levels
        covariatesDF[[D]] = factor(
          covariatesDF[[D]],
          levels = 1:nlevels(data[[D]]),
          labels = levels(data[[D]])
        )
        
      }
      
    } else { # data[[D]] isn't a factor
      
      # choose baseline category 
      
      
      theTable = sort(table(dataD), decreasing=TRUE)
      theTable = theTable[theTable > 0]
      if(is.null(theLevels)) {
        theLabels = paste("l", names(theTable),sep="")
      } else {
        
        if(all(names(theTable) %in% theLevels[,labelCol])) {
          # convert table names to numeric
          # code must work for data where data[[D]] is numeric
          names(theTable) = theLevels[
            match(names(theTable), theLevels[,labelCol])
            , idCol]
        }
        
        
        theLabels = as.character(theLevels[
          match(names(theTable), as.character(theLevels[,idCol])), 
          labelCol
        ])
        if(is.numeric(dataD)) {
          levelsD = as.integer(names(theTable))
        } else {
          levelsD = theLabels
        }
        if(any(is.na(theLabels))) {
          warning(
            'missing labels in covariate raster ', 
            D, ' level ',
            names(theTable)[is.na(theLabels)][1])
          theLabels[is.na(theLabels)] = 
            names(theTable)[is.na(theLabels)]
        }
      } # end else (not is null thelevels)
      
      # re-factor data with new baseline category
      terra::values(data)[[D]] = factor(
        dataD, 
        levels=levelsD,
        labels=theLabels)     
      covariatesDF[[D]] = factor(
        as.character(covariatesDF[[D]]), 
        levels=theLabels,
        labels=theLabels)     
    } # end refactor
  } # end loop D trhoguh factors
  
  
  list(
    data=data,
    grid = cellsSmall,
    covariates=covariatesDF
  )
}
