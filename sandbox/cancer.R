
#cancerRatesCanada = function(area = "canada",
#    year=2000, sex=c("M", "F"), site="Lung") {
#  
#
#  cpath = system.file('extdata', package='diseasemapping')
#  zfile = file.path(cpath, 'canCancer.zip')
#  zTempFile = file.path(tempfile(), 'canCancer.zip')
#  
#  # if neither file exists
#  if(!file.exists(zfile) & !file.exists(zTempFile)){
#    
#    # if we can write to extdata
#    if(file.access(cpath, 2)>=0){
#      dfile = zfile
#    } else {
#      dfile = zTempFile
#    }
#    message('downloading cancer data to', dfile)
#    download.file(canZip,dfile)
#  }
#  
#
#  cpath = tempdir()
#  
#  allFiles = utils::unzip(dfile, list=TRUE)$Name
#  utils::unzip(zfile, exdir=cpath, files=allFiles)
#
#  if(FALSE){
#  library('RMySQL')
#  
#  cClasses=c('numeric',rep('character',5),'NULL','NULL', 'numeric')
#  
#  smallData = read.table(file.path(cpath,allFiles), 
#      sep=',', header=TRUE,
#      stringsAsFactors=FALSE, 
#      na.strings='..', comment.char="",
#      nrows=20)
#
#  
#  library(RPostgreSQL)
#  drv = dbDriver("PostgreSQL")
#  con <- dbConnect(drv, dbname = "cancer")
#  dbWriteTable(con, name="canada", 
#      value=smallData,		
#      overwrite = TRUE)
#  
#  
#  cClasses=unlist(lapply(smallData, function(qq){
#            
#            dbDataType(dbDriver('MySQL'), qq)
#            
#          }))
#  
#  
#  con <- dbConnect(dbDriver("MySQL"), dbname = "cancer")
#  dbWriteTable(con, name="canada", 
#      value=allFiles, 
#      field.types=cClasses,		
#      overwrite = TRUE)
#  }
#  
#  x = NULL
#  Nrow = 200000
#  notDone = TRUE
#  
#  con=file(file.path(cpath, allFiles), 'r')
#
#  header = unlist(scan(con, sep=',', 
#          what=list('c','c','c','c','c','c','c','c', 'c'),
#          na.strings='..', comment.char="",
#          nmax=1))
#
#  while(notDone){
#    cat('.')
#  cData = read.table(con, sep=',',
#      stringsAsFactors=FALSE, 
#      colClasses=cClasses,
#      na.strings='..', comment.char="",
#      nrows=Nrow, col.names=header)
#  if(nrow(cData)<Nrow) notDone=FALSE
#  
#  toKeep = which(
#      cData$Ref_Date %in% year & cData$Value > 0 &
#          substr(cData$UNIT, 1, 2) =='Ne'
#  )
#  x = rbind(x,cData[toKeep,])
#  }
#   close(con)
#   
#   
#  # cancer types file
#  cData$cancers = utils::read.table(
#      file.path(cpath, 
#          grep('^cancer', allFiles, ignore.case=TRUE, value=TRUE)),
#      header=TRUE, sep="\t"
#      )
#  for(D in 1:ncol(cData$cancers))
#        cData$cancers[[D]] = gsub("^[[:space:]]+|[[:space:]]+$", '', cData$cancers[[D]])
#
#  # registry
#  cData$registry = utils::read.fwf(
#      file.path(cpath, 
#          grep('^registry', allFiles, ignore.case=TRUE, value=TRUE)),
#      widths=c(12,5,100),
##      fileEncoding='ISO_8859-15',
##      encoding='UTF-8',
#      stringsAsFactors=FALSE
#  )
#  for(D in 1:ncol(cData$registry))
#    cData$registry[[D]] = gsub("^[[:space:]]+|[[:space:]]+$", '', cData$registry[[D]])
#  names(cData$registry) = c('registry','extra', 'label')
#
#  # find country
#  Sarea = paste(area, collapse='|')
#  SregistryWeak = grep(Sarea, cData$registry$label, ignore.case=TRUE)
#  Sregistry = grep(paste(Sarea, '([[:space:]]?)$', sep=''), 
#      cData$registry$label, 
#      ignore.case=TRUE)
#  if(length(Sregistry)<length(area)){
#    Sregistry = SregistryWeak
#  }
#  
#  thisReg = cData$registry[Sregistry, ]
#  Sregistry = thisReg[,'registry']
#  names(Sregistry) = thisReg[,'label']
#  
#  if(!length(Sregistry)){
#    warning("can't find registry", area)
#  }
#  
#  # find cancer site
#
#  siteRegexp = paste(site, collapse='|')
#
#  Ssite = grep(paste(
#          '^[[:space:]]?(', siteRegexp,')[[:space:]]\\(C', sep=''),
#      cData$cancers$label, ignore.case=TRUE)
#  SsiteWeak = grep(siteRegexp, cData$cancers$label, ignore.case=TRUE)
#  
#  if(length(Ssite) < length(site))
#    Ssite = SsiteWeak
#  
#  if(!length(Ssite)){
#    warning('cant find cancer site', site)
#  }
#  
#  thisSite = cData$cancers[Ssite,]
#  Ssite = thisSite[,'cancer']
#  names(Ssite) = gsub("^[[:space:]]", "", thisSite[,'label'])
#  
#  
#  
#  # load cancer data
#  allCfiles = as.character(utils::unzip(zfile, list=TRUE)$Name)
#  cFiles = grep(paste(Sregistry, collapse='|'), 
#      allCfiles, 
#      value=TRUE)
#  utils::unzip(zfile, exdir=cpath, files=cFiles)
#  for(D in cFiles){
#    cData[[D]] = read.table(file.path(cpath, D), sep=',')
#    names(cData[[D]]) = c('sex','cancer','age','cases','pop')
#    cData[[D]]$age = 5*(cData[[D]]$age -1)
#    cData[[D]]$sex = factor(cData[[D]]$sex, levels=c(1,2), labels=c('M','F'))
#    cData[[D]] = cData[[D]][cData[[D]]$cancer %in% Ssite,]
#  }
#  
#  
#  
#  counts = cData$cases[
#      cData$cases$REGISTRY %in% Sregistry &
#          cData$cases$CANCER %in% Ssite,]
#
#  Sages = grep('^N([[:digit:]]|_unk)', names(counts), value=TRUE)
#
#  counts=reshape(counts, 
#      direction='long',
#      varying = list(
#          Sages
#      ),
#      v.names='count',
#      timevar='age',
#      times=gsub("^[[:alpha:]]", "", Sages),
#      idvar=c('REGISTRY','CANCER', 'SEX','ETHNIC_GROUP')
#      )
#      
#  
#  pop = cData$pop[
#      cData$pop$REGISTRY %in% Sregistry,]
#
#Sages = grep('^P([[:digit:]]|_unk)', names(pop), value=TRUE)
#
#pop=reshape(pop, 
#    direction='long',
#    varying = list(
#        Sages
#    ),
#    v.names='population',
#    timevar='age',
#    times=gsub("^[[:alpha:]]", "", Sages), 
#    idvar=c('REGISTRY', 'SEX','ETHNIC_GROUP')
#)
#
#x = merge(counts, pop, by=c('age','REGISTRY','SEX','ETHNIC_GROUP'))
#
#
#x$sex = factor(x$SEX, levels=1:2, labels=c('M','F'))
#x$rate = x$count / x$population
#x$group = paste(as.character(x$sex), x$age, sep='')
#
#x = x[!is.nan(x$rate),]
#
#x = x[,c('age','sex','group','REGISTRY','ETHNIC_GROUP', 'CANCER','rate')]
#
#x = reshape(x, direction='wide', 
#    timevar='CANCER', 
#    idvar=c('REGISTRY','age','sex','group','ETHNIC_GROUP')
#)
#
#x = x[as.character(x$sex) %in% sex,]
#
#
#Scancers = grep("^rate\\.[[:digit:]]+$", names(x), value=TRUE)
#names(Scancers) = Scancers
#Scancers = gsub("rate\\.", '', Scancers)
#
#cNames = names(Ssite)[match(Scancers, as.character(Ssite))]
#
#cNames = gsub('[[:punct:]]','',cNames)
#cNames = gsub('[[:space:]]','_',cNames)
#
#names(x)[match(names(Scancers),names(x) )] = cNames
#
#if(!any(duplicated(x$group))){
#  rownames(x) = as.character(x$group)
#  x = drop(as.matrix(x[,cNames]))
#}
#
#
#
#attributes(rates)$area = Sregistry
#attributes(rates)$year = year
#
#
# }
