
zFile = 'fsaPop.zip'
if(!file.exists(zFile))
	download.file(
			'http://www12.statcan.gc.ca/open-gc-ouvert/2011/98-311-XCB2011022.ZIP',
			zFile, mode='wb')

pFile = 'fsaPop.RData'

if(!'rsdmx' %in% installed.packages()[,'Package'])
	devtools::install_github("opensdmx/rsdmx")

unzip(zFile, exdir='.')
sFile = grep("Generic", unzip(zFile, list=TRUE)$Name, value=TRUE)

library(rsdmx)
fsaSdmx =  readSDMX(file=sFile, isURL=FALSE)

if(!file.exists(pFile)) {
	fsaDf = as.data.frame(fsaSdmx)
	save(fsaDf, file=pFile)
}

load(pFile)

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

fsaDf$sex = factor(fsaDf$Sex,
		levels = sexCode$id, 
		labels = substr(gsub("[[:space:]]", "", sexCode$label.en), 
				1,1))
fsaDf$age = as.character(factor(fsaDf$AGE,
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

zFile = 'fsaBorder.zip'
if(!file.exists(zFile))
	download.file(
			'http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gfsa000b11a_e.zip',
			zFile, mode='wb')

unzip(zFile, exdir='.')
sFile = grep("\\.shp$", unzip(zFile, list=TRUE)$Name, value=TRUE)
library(rgdal)
fsa = readOGR(".", gsub('\\.shp$','',sFile),stringsAsFactors=FALSE)


fsa$name0='Canada'
fsa$id0='CAN'

fsa$name1=fsa$PRNAME
fsa$id1= fsa$PRUID
fsa$id1.1 = as.character(fsa$CFSAUID)

fsa@data = fsa@data[,grep("^(name|id)[[:digit:]]", 
				names(fsa), value=TRUE)]

fsa = spTransform(fsa, CRS("+init=epsg:4326"))

toAdd = fsa@data[
		match(
				as.character(fsaPop$id1.1),
				as.character(fsa$id1.1)
		),
]
fsaPopFull = cbind(toAdd, fsaPop)		

if(FALSE) {
	fsa2011 = fsa
	fsa2011@data = fsaPopFull[
			match(fsa2011$id1.1, fsaPopFull$id1.1),
			]
	save(fsa2011, file='fsa2011.RData', compress='xz')
}

library(mapmisc)
fsaBg = openmap(fsa)

map.new(fsa, buffer=-c(1,1,0.5,5))
plot(fsaBg, add=TRUE)
plot(fsa, add=TRUE, border='red')

if(FALSE){
	theDir = "/store/census/CAN/2011/1.1"
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






if(FALSE) {
	# are fsa's nested within CD's?
	
	# cd boundaries
	
	zFile = 'cdBorder.zip'
	if(!file.exists(zFile))
		download.file(
				'http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gcd_000a11a_e.zip',
				zFile, mode='wb')
	unzip(zFile, exdir='.')
	sFile = grep("\\.shp$", unzip(zFile, list=TRUE)$Name, value=TRUE)
	library(rgdal)
	cd = readOGR(".", gsub('\\.shp$','',sFile),stringsAsFactors=FALSE)
	
	
	library(rgeos)
	fsaS = gSimplify(fsa, tol=0.0002)
	c(length(fsaS), length(fsa))
	cdS = gSimplify(cd, tol=0.001)
	c(length(cdS), length(cd))
	fsacd = gIntersection(fsaS, cdS, 
			byid=c(TRUE, TRUE),
			drop_lower_td=TRUE)
	length(fsa)
	length(fsacd)
	fsacd = SpatialPointsDataFrame(
			fsacd,
			data.frame(
					index1 = as.numeric(gsub(" [[:digit:]]+$", "", names(fsacd))),
					index2 = as.numeric(gsub("^[[:digit:]]+ ", "", names(fsacd)))
			)
	)
	fsacd$cd = cd$CDUID[
			match(
					as.character(fsacd$index2),
					rownames(cd@data)
			)
	]		
	fsacd$fsa = fsa$CFSAUID[
			match(
					as.character(fsacd$index1),
					rownames(fsa@data)
			)
	]		
	
	themax = names(which.max(table(fsacd$fsa)))
	themax = fsa[fsa$CFSAUID == themax,]
	library(mapmisc)
	bgmap = openmap(themax)
	
	
	map.new(themax,
			buffer=-0.01
	)
	plot(bgmap, add=TRUE)
	plot(themax, add=TRUE, col='orange', lwd=4, lty=0)
	plot(fsa, add=TRUE, border='red', lwd=1)
	plot(cd, add=TRUE, border='blue', lty=3)
	text(
			SpatialPoints(fsa),
			fsa$CFSAUID, cex=0.6
	)
	# no, not nested.  not disjoint either		
}


