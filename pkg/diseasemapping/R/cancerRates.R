ci5zip = 'http://ci5.iarc.fr/CI5plus/old/CI5plus.zip'
canZip = 'http://www20.statcan.gc.ca/tables-tableaux/cansim/csv/01030550-eng.zip'

cancerRates = function(area = "canada",
   year=2000, sex=c("M", "F"), site="Lung") {

  sexes = c("M"=1, "F"=2)[sex]



  areaCodes =
  c("Canada"=1240099,
  "Newfoundland"=1240899,
"Prince Edward Island"=1240799,
"Nova Scotia"=1240699,
"Ontario"=1241199,
"Manitoba"=1240399,
"Saskatchewan"=1241399,
"Alberta"=1240199,
"British Columbia"=1240299,
"New Zeland"=5540099,
"Sweden"=7520099,
"Slovenia"=7050099,
"Slovakia"=7030099,
Norway=5780099,
Latvia=4280099,
Lithuania=4400099,
Iceland=3520099,
Finland=2460099,
Estonia=2330099,
Denmark=2080099,
"Czech Republic"=2030099,
"Costa Rica"=1880099,
USA=8400199,
Iowa=8400899,
"New Mexico"=8401399
)
if(is.character(area)) {
area = areaCodes[grep(paste("^", area[1], sep=""), 
				names(areaCodes), ignore.case=TRUE)]
} 


result = list()

rates=NULL
for(Dsex in names(sexes)) {
fs<-paste("http://ci5.iarc.fr/CI5plus/old/Table4r.asp?registry=",area,(paste("&period=",year,sep="",collapse="")),"&sex=", sexes[Dsex],"&window=1&text=1&stat=0&submit=Execute",sep="")
tempn = scan(fs, what="a", quiet=TRUE)
theurl=(paste("http://ci5.iarc.fr/", 
					gsub("^HREF=", "", grep("href=/data", 
									iconv(tempn,to="UTF-8"), value=TRUE, 
									ignore.case=TRUE)), sep=""))
theurl = url(theurl)
forAttribute = scan(theurl, what="a", sep="\t", nmax=1, quiet=TRUE)
result[[Dsex]]=read.table(theurl, header=TRUE,skip=1, 
		fill=TRUE, sep="\t", as.is=TRUE)
}

iarcSite=NULL
for (Dsex in names(sexes)) {
 cancerTypes =  result[[Dsex]][,1]
  
 cancerTypes = gsub("[[:space:]]+$", "", cancerTypes)

siterow = grep(site, cancerTypes,ignore.case=TRUE)
if(length(siterow) > 1 )    {
  warning(paste("matched ", paste(cancerTypes[siterow], collapse=","), "\n from cancer types", paste(cancerTypes, collapse=","), 
  "\n using", cancerTypes[siterow[1]]))
  siterow=siterow[1]
  }
if(length(siterow) ) {
  x=result[[Dsex]][siterow,]
iarcSite =gsub("[[:space:]]+$", "", x[1]) 
  x = x[,-c(1,grep("Total|Unk|ICD|CR|ASR", 
						  names(x), ignore.case=TRUE))]
  x = as.vector(x)
  names(x) = gsub("^X|\\.$", "", names(x))
  rates[[Dsex]]=x
  }
}
rates = unlist(rates)/100000
names(rates) = gsub("\\.", "_", names(rates))

attributes(rates)$site = iarcSite
forAttribute = strsplit(forAttribute, "\\(")[[1]]
attributes(rates)$area = gsub("[[:space:]]+$", "", forAttribute[1])
attributes(rates)$year = gsub("\\)$", "", forAttribute[length(forAttribute)])

return(rates)

}

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
