
zFile = 'fsaPop2006.zip'
if(!file.exists(zFile))
	download.file(
			'http://www12.statcan.gc.ca/open-gc-ouvert/2006/97-551-XCB2006010.ZIP',
			zFile, mode='wb')

pFile = 'fsaPop2006.RData'

if(!'rsdmx' %in% installed.packages()[,'Package'])
	devtools::install_github("opensdmx/rsdmx")

unzip(zFile, exdir='.')
sFile = grep("Generic", unzip(zFile, list=TRUE)$Name, value=TRUE)

library(rsdmx)
fsaSdmx =  readSDMX(file=sFile, isURL=FALSE)
otherFile = grep("Generic", unzip(zFile, list=TRUE)$Name, 
		value=TRUE, invert=TRUE)
otherSdmx =  readSDMX(file=otherFile, isURL=FALSE)
codeList = otherSdmx@codelists
sapply(slot(codeList, "codelists"), function(x) slot(x, "id"))
sexCode = as.data.frame(
		slot(otherSdmx, "codelists"), 
		codelistId = c("CL_SEX")) 
ageCode = as.data.frame(
		slot(otherSdmx, "codelists"), 
		codelistId = c("CL_AGE")) 

if(!file.exists(pFile)) {
	fsaDf = as.data.frame(fsaSdmx)
	save(fsaDf, file=pFile)
}

load(pFile)


fsaDf$sex = factor(fsaDf$Sex,
		levels = sexCode$id, 
		labels = substr(gsub("[[:space:]]", "", sexCode$label.en), 
				1,1))
fsaDf$age = as.character(factor(fsaDf$Age,
				levels = ageCode$id, labels = 
						gsub("[[:punct:]]|[[:space:]]|years?|Age", "",
								ageCode$label.en)))
fsaDf$age = gsub("to", "_", fsaDf$age)
fsaDf$age = gsub("andover", "plus", fsaDf$age)
fsaDf$age = gsub("Total", "all", fsaDf$age)
fsaDf$ageLow = as.numeric(gsub("_[[:digit:]]+$", "",  fsaDf$age))
fsaDf$ageHigh = as.numeric(gsub("^[[:digit:]]+_", "",  fsaDf$age))
fsaDf$ageDif = fsaDf$ageHigh - fsaDf$ageLow

toKeep = union(which(
				fsaDf$ageDif==4
		),
		grep("^all|^100", fsaDf$age)		
)
fsaDf = fsaDf[sort(toKeep),]

fsaDf$group = paste(
		fsaDf$sex, fsaDf$age, sep=''
)
fsaDf$id1.1 = gsub("^[[:digit:]]{2}", "", fsaDf$GEO)

fsaPop = reshape(
	fsaDf, v.names='obsValue',
	idvar = 'id1.1',
	timevar = 'group',
	direction = 'wide',
	drop = grep(
			"^group$|^id[[:digit:]]|obsValue",
			names(fsaDf), invert=TRUE, value=TRUE
			)
)
names(fsaPop) = gsub("^obsValue\\.", "", names(fsaPop))

# fsa boundaries

zFile = 'fsaBorder2006.zip'
if(!file.exists(zFile))
	download.file(
			'http://data.library.utoronto.ca/dataut/cc06/georef/cbf/arcinfo/gfsa000b06a_e.zip',
			zFile, mode='wb')

unzip(zFile, exdir='.')
sFile = grep("\\.shp$", unzip(zFile, list=TRUE)$Name, value=TRUE)
library(rgdal)
fsa = readOGR(".", gsub('\\.shp$','',sFile),stringsAsFactors=FALSE)


fsa$name0='Canada'
fsa$id0='CAN'

fsa$id1= fsa$PRUID
# fsa$name1=fsa$PRNAME
fsa$id1.1 = as.character(fsa$CENSUS_FSA)

fsa@data = fsa@data[,grep("^(name|id)[[:digit:]]", 
				names(fsa), value=TRUE)]

fsa = spTransform(fsa, CRS("+init=epsg:4326"))

toAdd = fsa@data[
		match(
				as.character(fsaPop$id1.1),
				as.character(fsa$id1.1)
		),
]
fsaPopFull = merge(toAdd, fsaPop, by='id1.1')		

if(FALSE) {
	fsa2006 = fsa
	fsa2006@data = fsaPopFull[
			match(fsa2006$id1.1, fsaPopFull$id1.1),
			]
	save(fsa20061, file='fsa2006.RData', compress='xz')
}

library(mapmisc)
fsaBg = openmap(fsa)

map.new(fsa, buffer=-c(1,1,0.5,5))
plot(fsaBg, add=TRUE)
plot(fsa, add=TRUE, border='red')

if(FALSE){
	theDir = "/store/census/CAN/2006/1.1"
	system(paste("mkdir", theDir))
	png(file.path(theDir,"map.png"))
	map.new(fsa, buffer=-c(1,1,0.5,5))
	plot(fsaBg, add=TRUE)
	plot(fsa, add=TRUE, border='red')
	dev.off()
	
	writeOGR(
			fsa,
			theDir,
			'map',
			driver= "ESRI Shapefile"
			)
	
	foreign::write.dbf(
			fsaPopFull,
			file.path(
					theDir, 'population.dbf'
					)
			)			
	
	
}





