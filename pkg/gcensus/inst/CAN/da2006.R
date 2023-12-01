Dyear = 2006
Dyear = 2001
Dlevel = 5
basedir = '/store/census'
rawdir = file.path(basedir, 'raw')

# only ivt files before 2001
# 1996 https://www12.statcan.gc.ca/English/census96/data/tables/Rp-eng.cfm?LANG=E&APATH=3&DETAIL=1&DIM=0&FL=A&FREE=1&GC=0&GID=0&GK=0&GRP=1&PID=1032&PRID=0&PTYPE=89103&S=0&SHOWALL=No&SUB=0&Temporal=2006&THEME=1&VID=0&VNAMEE=&VNAMEF=
# 1991 https://www12.statcan.gc.ca/English/census91/data/tables/Rp-eng.cfm?LANG=E&APATH=3&DETAIL=1&DIM=0&FL=A&FREE=1&GC=0&GID=0&GK=0&GRP=1&PID=86&PRID=0&PTYPE=4&S=0&SHOWALL=No&SUB=0&Temporal=1991&THEME=101&VID=0&VNAMEE=&VNAMEF=

if(Dyear == 2001) {
	zUrl = 'https://www12.statcan.gc.ca/open-gc-ouvert/2001/95F0300XCB2001001.ZIP'
} else {
# population
zUrl = paste0('https://www12.statcan.gc.ca/census-recensement/',
	Dyear, '/dp-pd/prof/rel/OpenDataDownload.cfm?PID=89109')
}
zFile = file.path(rawdir, paste('da',Dyear,'pop.zip', sep=''))

if(!file.exists(zFile))
	download.file(zUrl, zFile, method='libcurl')

cFile = unzip(zFile, exdir=rawdir, 
		unzip='/usr/bin/unzip')

cFile = unzip(zFile, list=TRUE, 
		unzip='/usr/bin/unzip')

xx = readsdmx::read_sdmx(file.path(rawdir,cFile[1,'Name']))
#xx2 = readsdmx::read_sdmx(file.path(rawdir,cFile[2,'Name']))
table(nchar(xx$GEO))
xx = xx[nchar(xx$GEO)== 11, ]

SageNum = seq(0, 85, by=5)
Sage = c('total', SageNum)
Ssex = c('m','f')

SageSex = expand.grid(age=Sage, sex=Ssex)
SageSex = paste0(SageSex[,'sex'], '_', SageSex[,'age'])

Svars = c('total', 'totalRound',SageSex)


if(Dyear == 2001) ) {
	xx$Sex = factor(xx$DIM1, levels = 1:3, labels = c('total','m','f'))	
	xx$ObsValue = as.numeric(as.character(xx$ObsValue))
	forLabels = XML::readHTMLTable(
		RCurl::getURL('https://www12.statcan.gc.ca/English/census01/products/standard/themes/Rp-eng.cfm?LANG=E&APATH=3&DETAIL=1&DIM=0&FL=A&FREE=1&GC=0&GID=0&GK=0&GRP=1&PID=55434&PRID=0&PTYPE=55430&S=0&SHOWALL=No&SUB=0&Temporal=2006&THEME=37&VID=0&VNAMEE=&VNAMEF=')
	)$tabulation
	xx$AgeGroup = factor(xx$DIM0, levels = 1:nrow(forLabels), labels = forLabels$V1)


	xxSub = xx[grep('[-]|[+]', xx$AgeGroup), ]
	xxSub = xxSub[grep("Total", xxSub$AgeGroup, invert=TRUE), ]
	xxSub = xxSub[grep("total", xxSub$Sex, invert=TRUE), ]
	xxSub$AgeGroup = factor(as.character(xxSub$AgeGroup))
	xxSub$Age1 = as.numeric(gsub("([-]|[+]).*", "", xxSub$AgeGroup))

	xxSub$AgeCut = cut(xxSub$Age1, c(SageNum-0.1, 200), labels = SageNum)
	xxSub2 = aggregate(xxSub[,'ObsValue', drop=FALSE], xxSub[,c('GEO','Sex', 'AgeCut')], sum, na.rm=TRUE)
	xxSub2$Age = as.numeric(as.character(xxSub2$AgeCut))

	xxSub2$var = paste0(xxSub2$Sex, '_', xxSub2$Age)
	xx = xxSub2
} else {



xx$var = Svars[as.numeric(as.character(xx$A06_SexAge49_D1))]
}

xx2 = reshape2::dcast(xx, GEO ~ var, value.var='ObsValue')
names(xx2) = gsub("^GEO$", "idOrig", names(xx2))
xx2$id1 = substr(xx2$id, 1,2)
xx2$id = xx2$id5 = paste0(substr(xx2$idOrig, 1,4), substr(xx2$idOrig, 8,11))

xx2 = xx2[,intersect(c('id','id1','id5', 'idOrig',  Svars), names(xx2))]


foreign::write.dbf(xx2, 
		file.path(basedir, 'CAN', Dyear, Dlevel,"pop.dbf")
)

toCheck = foreign::read.dbf(file.path(basedir, 'CAN',Dyear, Dlevel, 'pop.dbf'))
toCheck[1:2,]

# check with shapefile

library('terra')
xSpatial = vect(file.path(basedir, 'CAN', Dyear, Dlevel, 'map.shp'))
xSpatial2 = merge(xSpatial, xx2, all.x=TRUE)




