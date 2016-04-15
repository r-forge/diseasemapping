
library(rgdal)
library("mapmisc")
rawDir = "/store/census/raw"

# zip 2000 2010
zUrl = 'https://data2.nhgis.org/extracts/92245/19/nhgis0019_shape.zip'

pzUrl = 'https://data2.nhgis.org/extracts/92245/19/nhgis0019_csv.zip'

zFile = file.path(rawDir,basename(zUrl))
pzFile = file.path(rawDir,basename(pzUrl))
if(!file.exists(zFile))
	download.file(zUrl, zFile)
if(!file.exists(pzFile))
	download.file(pzUrl, pzFile)

fNames = unzip(zFile, list=TRUE)$Name
fNames

pfNames = grep("csv$", unzip(pzFile, list=TRUE)$Name, value=TRUE)
pfNames


for(D in pfNames) {
	dFile = file.path(rawDir, D)
	if(!file.exists(dFile)) {
		unzip(pzFile, exdir=rawDir)
	}
	pop = read.table(dFile, header=TRUE, 
		sep=',',stringsAsFactors=FALSE,
		skip=1)	
	names(pop) = gsub("(em)?ale\\.+|years|", "", names(pop))
	names(pop) = gsub("\\.+(and|to)\\.+", "_", names(pop))
	names(pop) = gsub("_over\\.?$", "plus", names(pop))
	names(pop) = gsub("\\.$", "", names(pop))
	names(pop) = gsub("Under\\.", "0_", names(pop))

	idCol = grep("5.digit.zip.code", names(pop), ignore.case=TRUE, value=TRUE)
	pop$id0 = 'USA'
	pop$id1 = paste("USA", substr(pop[[idCol]],1,2), sep='')
	pop$id = pop$id1.2 = paste("USA", pop[[idCol]], sep='')
	pop = pop[,grep("^(M[[:digit:]]|F[[:digit:]]|id)", names(pop))]
	
	
	Dyear = gsub("_zcta.csv$", "", D)
	Dyear = gsub(".*_", "", Dyear)
	
	
	foreign::write.dbf(
			pop,
			file.path("/store/census/USA", Dyear, "1.2", "pop.dbf")
			)
	
}

for(D in fNames) {
	if(!file.exists(D)) {
		unzip(zFile, exdir=rawDir)
	}
	Dnames = unzip(D, list=TRUE)$Name	
	if(!all(file.exists(Dnames))){
		unzip(D, exdir=rawDir)	
	}
	x = readOGR(rawDir, 
		gsub("\\.shp$", "", 
			grep("shp$", Dnames, value=TRUE)[1]
		)
	)
	x$id0 = 'USA'
	zcol = grep("^ZCTA", names(x), value=TRUE)
	x$id1 = paste("USA", substr(x[[zcol]],1,2), sep='')
	x$id = x$id1.2 = paste("USA", x[[zcol]], sep='')

	x = spTransform(x, CRS("+init=epsg:4326"))
	
	xDf = x@data[,grep("^id", names(x))]
	
	x2 = x
	x2@data = xDf
	
	Dyear = gsub(".*_", "", D)
	Dyear = gsub("\\.zip$", "", Dyear)
	Dpath = file.path("/store/census/USA", Dyear, "1.2")
	writeOGR(x2, Dpath, 'map', driver="ESRI Shapefile", 
		overwrite_layer = file.exists(file.path(Dpath, "map.shp"))
	)
	
	bgMap = openmap(x2)
	png(file.path(Dpath, "map.png"), height=600,width=800)
	map.new(x2)
	plot(bgMap,add=TRUE)
	plot(x2, add=TRUE)
	dev.off()
	
}


