
Dyear = 2021
Dlevel = 5
basedir = '/store/census'
rawdir = file.path(basedir, 'raw')

# population
zUrl = 'https://www150.statcan.gc.ca/n1/en/tbl/csv/98100023-eng.zip?st=epPACv5w'

zFile = file.path(rawdir, paste('da',Dyear,'pop.zip', sep=''))

if(!file.exists(zFile))
	download.file(zUrl, zFile, method='libcurl')

cFile = unzip(zFile, exdir=rawdir, 
		unzip='/usr/bin/unzip')

cFile = unzip(zFile, list=TRUE, 
		unzip='/usr/bin/unzip')

Sfiles = grep("CSV$", cFile$Name, value=TRUE, ignore.case=TRUE)
Sfiles = grep("Metadata|-DQ[.]", Sfiles, value=TRUE, invert=TRUE, ignore.case=TRUE)


	daPopOrig = daPop = read.csv(file.path(rawdir,Sfiles[1]), stringsAsFactors=FALSE)
	names(daPop) = gsub("Age.*", "AgeGroup",  names(daPop))
	names(daPop) = gsub(".*Total.*", "Total",  names(daPop))
	names(daPop) = gsub(".*Men.*", "m",  names(daPop))
	names(daPop) = gsub(".*Women.*", "f",  names(daPop))


	daPop$age = as.integer(gsub(" .*", "", daPop$AgeGroup))
	daPop = daPop[nchar(daPop$DGUID) == 17, ]

	daPop$other = daPop$Total - daPop$m - daPop$f
	daPop=daPop[,grep("Symbol", names(daPop), invert=TRUE)]

	theOthers = daPop[grep("Total", daPop$AgeGroup), ]

	daPop = daPop[grep(" to | and ", daPop$AgeGroup), ]
	daPop = daPop[grep("^[89]. to|^100|15 to 64|65 years and over|0 to 14", daPop$AgeGroup, invert=TRUE), ]


	daPop2 = reshape2::melt(daPop[,c('DGUID','age','m', 'f')], id.var = c('DGUID','age'))
	daPopWide = reshape2::dcast(daPop2, DGUID ~ variable + age)
	daPopWide$other = theOthers[match(daPopWide$DGUID, theOthers$DGUID), 'other']
	daPopWide$other = pmax(0, daPopWide$other, na.rm=TRUE)


theMap = terra::vect(file.path(basedir, 'CAN', Dyear, Dlevel, 'map.shp'))
daPopWide$id5 = values(theMap)[match(daPopWide$DGUID, theMap$DGUID), 'id5']
daPopWide = daPopWide[!is.na(daPopWide$id5), ]

foreign::write.dbf(daPopWide, 
		file.path(basedir, 'CAN', Dyear, Dlevel,"pop")
)



