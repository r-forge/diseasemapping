# 2001 FSA's

# boundaries
#zUrl = 'http://data.library.utoronto.ca/dataut/cc01/georef/spatial/arcinfo/cbf/national/fsa_cbf_en.exe'
#zFile = 'fsa2001arc.exe'

zUrl = 'http://data.library.utoronto.ca/dataut/cc01/georef/spatial/mapinfo/cbf/national/fsa_cbf_en.exe'
zFile = 'fsa2001mapinfo.exe'
if(!file.exists(zFile))
	download.file(zUrl, zFile)

unzip(zFile, exdir='.')
unzip(zFile, list=TRUE)
mifFile = grep('MIF$',unzip(zFile, list=TRUE)$Name, value=TRUE)

library('rgdal')
ogrInfo('.',mifFile)
system(paste('ogr2ogr -f "ESRI Shapefile" -a_srs EPSG:4326',
  			'fsa2001mapinfo', mifFile))

fsa2001 = readOGR("fsa2001mapinfo",'FSA_CBF_EN')

fsa2001$id1.1 = as.character(fsa2001$CENSUS_FSA)
fsa2001$id1 = as.character(fsa2001$PR)
fsa2001$id0 = 'CAN'
fsa2001$name0 = 'Canada'

fsa = fsa2001[,grep("^(id|name)[[:digit:]]", names(fsa2001))]

# populations
zFile = 'fsaPop2001.zip'
if(!file.exists(zFile))
	download.file(
			'http://www12.statcan.gc.ca/open-gc-ouvert/2001/95F0300XCB2001005.ZIP',
			zFile, mode='wb')

pFile = 'fsaPop2001.RData'



unzip(zFile, exdir='.')
sFile = grep("Generic", unzip(zFile, list=TRUE)$Name, value=TRUE)

library(rsdmx)
fsaSdmx =  readSDMX(file=sFile, isURL=FALSE)

if(!file.exists(pFile)) {
	fsaDf = as.data.frame(fsaSdmx)
	save(fsaDf, file=pFile)
}


otherFile = grep("Generic", unzip(zFile, list=TRUE)$Name, 
		value=TRUE, invert=TRUE)
otherSdmx =  readSDMX(file=otherFile, isURL=FALSE)
codeList = otherSdmx@codelists
sapply(slot(codeList, "codelists"), function(x) slot(x, "id"))
sexCode = as.data.frame(
		slot(otherSdmx, "codelists"), 
		codelistId = c("CL_DIM1")) 
ageCode = as.data.frame(
		slot(otherSdmx, "codelists"), 
		codelistId = c("CL_DIM")) 

fsaCode = as.data.frame(
		slot(otherSdmx, "codelists"), 
		codelistId = c("CL_GEO")) 


fsaDf$sex = factor(fsaDf$DIM1,
		levels = sexCode$id, 
		labels = substr(gsub("[[:space:]]", "", sexCode$label.en), 
				1,1))
fsaDf$age = as.character(
			factor(fsaDf$DIM0,
			levels = ageCode$id, 
			labels = gsub("[[:space:]]", "", ageCode$label.en))
)

fsaDf$id1.1 = as.character(
		factor(fsaDf$GEO,
				levels = fsaCode$id, 
				labels = gsub("[[:space:]]", "", fsaCode$label.default))
)

fsaDf$age = gsub("Total.+Age$", "all", fsaDf$age)

fsaDf$age = gsub("-", "_", fsaDf$age)
fsaDf$age = gsub("\\+$", "plus", fsaDf$age)
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


toAdd = fsa@data[
		match(
				as.character(fsaPop$id1.1),
				as.character(fsa$id1.1)
		),
]
fsaPopFull = cbind(toAdd, fsaPop)		


if(FALSE) {
	fsa2001 = fsa
	fsa2001@data = fsaPopFull[
			match(fsa2001$id1.1, fsaPopFull$id1.1),
	]
	save(fsa2001, file='fsa2001.RData', compress='xz')
}

library(mapmisc)
fsaBg = openmap(fsa)

if(FALSE){
	theDir = "/store/census/CAN/2001/1.1"
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

fsaPopSub = 	fsaPopFull[grep('^[A-Z][[:digit:]][A-Z]$',fsaPopFull$id1.1),]
	foreign::write.dbf(
fsaPopSub,
			file.path(
					theDir, 'population.dbf'
			)
	)			
	
	
}



# 1996

zFile = 'fsaBorder1996.zip'
if(!file.exists(zFile))
	download.file(
			#'http://geo2.scholarsportal.info/proxy.html?http:__maps2.scholarsportal.info/files/zips/DLI/DLI_1996_Census_CBF_Eng_Nat_fsa.zip',
			'http://maps2.scholarsportal.info/files/zips/DLI/DLI_1996_Census_CBF_Eng_Nat_fsa.zip',
			zFile, mode='wb')

unzip(zFile, exdir='.')
unzip(zFile, list=TRUE)
zFile2 = unzip(zFile, list=TRUE)$Name[
		which.max(unzip(zFile, list=TRUE)$Length)]
unzip(zFile2, exdir='.')
unzip(zFile2, list=TRUE)
mFile = grep('TAB$',
		unzip(zFile2, list=TRUE)$Name, 
		value=TRUE)


library('rgdal')
ogrInfo('.',mFile)
system(paste('ogr2ogr -f "ESRI Shapefile" -a_srs EPSG:4326',
  			'fsa1996mapinfo', mFile))
fsa1996 = readOGR('fsa1996mapinfo', 'GFSA000B')
fsa1996 = spTransform(fsa1996, CRS("+init=epsg:4326"))

prov = raster::getData("GADM", country='CAN', level=1)
prov = spTransform(prov, fsa1996@proj4string)
fsaP = SpatialPoints(fsa1996, proj4string=fsa1996@proj4string)
provId = over(fsaP, prov)

fsa1996$id1 = as.character(provId$CCA_1)
fsa1996$name1 = as.character(provId$NAME_1)
fsa1996$id0 = 'CAN'
fsa1996$name0 = 'Canada'


zFile = 'fsaPop1996.zip'
if(!file.exists(zFile))
	download.file(
			'http://www12.statcan.ca/open-gc-ouvert/1996/95F0184XDB.ZIP',
			zFile)



pFile = 'fsaPop1996.RData'


unzip(zFile, exdir='.')
sFile = grep("Generic", unzip(zFile, list=TRUE)$Name, value=TRUE)

library(rsdmx)
fsaSdmx =  readSDMX(file=sFile, isURL=FALSE)

if(!file.exists(pFile)) {
	fsaDf = as.data.frame(fsaSdmx)
	save(fsaDf, file=pFile)
}




library(mapmisc)
fsaBg = openmap(fsa1996)

if(FALSE){
	theDir = "/store/census/CAN/1996/1.1"
	system(paste("mkdir", theDir))
	png(file.path(theDir,"map.png"))
	map.new(fsa1996, buffer=-c(1,1,0.5,5))
	plot(fsaBg, add=TRUE)
	plot(fsa1996, add=TRUE, border='red')
	dev.off()
	
	writeOGR(
			fsa1996,
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

