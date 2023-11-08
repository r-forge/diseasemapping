

basedir = '/store/census'
rawdir = file.path(basedir, 'raw')
dir.create(rawdir)
dir.create(file.path(basedir, 'CAN'))




library("terra")

data('worldMap', package='mapmisc');worldMap = unwrap(worldMap)

canadaLL = geodata::gadm(country='CAN', level=0, path=rawdir)
canada = terra::project(canadaLL,terra::crs("epsg:3347"))

canadaBgLL = mapmisc::openmap(canadaLL)

canadaBg = mapmisc::openmap(canada, zoom=1)



Sid = 		c(
		id1 = 'pr',
		id1.2 = 'cma',
		id2 = 'cd',
		id3 = 'csd',
		id4 = 'ct',
		id5 = 'da',
		id6 = 'db'
)

SidEx = paste(Sid,
		c("_", "")[nchar(Sid)-1],
		sep=''
)

SidNames = list(
		'1' = c(
				id1 = 'PRUID',
				name1 = 'PRENAME',
				fr1 = 'PRFNAME'
		),
		'1.2' = c(
				id1.2 = 'CMAPUID',
				name1.2 = 'CMANAME',
				cmatype = 'CMATYPE'
		),
		'2' = c(
				id2 = 'CDUID',
				name2 = 'CDNAME',
				cdtype = 'CDTYPE'
		),
		'3' = c(
				id3 = 'CSDUID',
				name3 = 'CSDNAME',
				type3 = 'CSDTYPE',
				sactype = 'SACTYPE'
		),
		'4' = c(
				id4 = 'CTUID'
		),
		'5' = c(
				id5 = 'DAUID'
				),
		'6' = c(
				id6 = 'DBUID'
				)
)



Syear = c(2001, 2006, 2011, 2016, 2021)  
Sprefix = c('2001' = 'g', '2006' = 'g', '2011' = 'g', '2016' = 'l', '2021'='l')
SinsideFolder = c('2001' = '', '2006' = '', '2011' ='' , '2016' = '2016/', '2021'='')
bndList = list()
#for(Dyear in Syear) {
for(Dyear in 2001) {
candir = file.path(basedir, 'CAN', Dyear)

dir.create(candir)

Surl = paste(
		'http://www12.statcan.gc.ca/census-recensement/', 2011, 
		'/geo/bound-limit/files-fichiers/',
#		Dyear, '/',
		SinsideFolder[as.character(Dyear)],
		Sprefix[as.character(Dyear)],
		SidEx, '000b', substr(Dyear, 3,4),
		'a_e.zip', sep=''
)

if(Dyear == 2001) {
	Surl = gsub("gdb_000b01a", "gdb_000a01a", Surl)
}

if(Dyear == 2021) {
	Surl = gsub("bound-limit", "sip-pis/boundary-limites", Surl)
}

Zfile = file.path(rawdir, basename(Surl))

#for(Dlevel in 1:length(Zfile)) {
#if(!file.exists(Zfile[Dlevel])) 
#	download.file(Surl[Dlevel], Zfile[Dlevel])
#}
Sfile = Pmisc::downloadIfOld(
		Surl, path=rawdir,
		#file=Zfile, 
		verbose=TRUE, exdir=rawdir
)



SshpFile = grep("[.]shp$", Sfile, value=TRUE)
SshpFile = gsub("[.]shp$", "", basename(SshpFile))

Se00files =  grep("[.]e00$", Sfile, value=TRUE)
Se00files = Se00files[nchar(Se00files)!= max(nchar(Se00files))]



bndList[[Dyear]] = list()
for(Dlevel in names(Sid)) {
	cat(Dlevel, ' ')
	DlevelId = as.numeric(gsub("id", "", Dlevel))
	leveldir = file.path(candir, DlevelId)
	dir.create(leveldir, showWarnings=FALSE)
	
if(!length(SshpFile)) {

	De00file =  grep(
			paste("[a-z]", Sid[Dlevel], sep=''),
			Se00files, value=TRUE
	)
	bnd1 = terra::vect(De00file, layer='ARC')
	bndLabels = terra::vect(De00file, layer='LAB')
	if(FALSE) {
		mapmisc::map.new(extent())
		plot(bnd1, add=TRUE)
	}

	values(bnd1) =  values(bndLabels)[bnd1$LPOLY_,grep("NAME|UID", names(bndLabels))]

	bnd = aggregate(bnd1, by=grep(paste0(Sid[Dlevel],'UID'), names(bnd1), ignore.case=TRUE, value=TRUE))
	
} else {


	DshpFile = grep(
			paste("[a-z]", Sid[Dlevel], sep=''),
			SshpFile, value=TRUE
	)
	


	bnd = terra::vect(file.path(rawdir, paste0(DshpFile, '.shp')))
}
	bnd = terra::project(bnd, terra::crs("epsg:3347"))
	
	bndDf = data.frame(
			id0 = rep('CAN', length(bnd)),
			name0 = 'Canada'
	)
	
	for(DlevelNames in names(SidNames)) {
		
		idHere = SidNames[[as.character(DlevelNames)]]
		
		for(Dname in names(idHere)) {
			if(idHere[Dname] %in% names(bnd))
				bndDf[[Dname]] = values(bnd)[[idHere[Dname]]]
		}
		
	}
	bndS = terra::simplifyGeom(bnd, tol=50)
#	rownames(bndDf) = names(bndS)
	for(D in names(bndDf)) bndS[[D]] = bndDf[[D]]
	bndList[[Dyear]][[Dlevel]] = bndS	
	# write map
	
	png(file.path(leveldir, "map.png"))
	mapmisc::map.new(canada, buffer=-c(0,0.5,0.5,5))
	plot(canadaBg, add=TRUE)
	plot(bndS, add=TRUE, border='red')
	dev.off()
	
	# write shapefile

	writeVector(
			bndS,
			file.path(leveldir, 'map.shp'),
			filetype= "ESRI Shapefile", overwrite=TRUE
	)
	
	
	

}

}
# population
Dyear = 2016

# canada, province, cd, csd da.  my DA is 5, ct is 4
SlevelsRecode = c('0' = 0, '1' =1, '2' = 2, '3'=3, '4' = 5)



rawDirDyear = file.path(rawdir, Dyear)
dir.create(rawDirDyear)
popFile = Pmisc::downloadIfOld(
'https://www12.statcan.gc.ca/census-recensement/2016/dp-pd/dt-td/CompDataDownload.cfm?LANG=E&PID=109525&OFT=CSV',
path = rawDirDyear,
file=file.path(rawDirDyear,'pop.zip'),
exdir = rawDirDyear)
popFile = grep("data.csv$", popFile, value=TRUE)
popHere = read.csv(popFile)
popHere$level = SlevelsRecode[as.character(popHere$GEO_LEVEL)]
popHere = popHere[grep("^[[:digit:]]+$|Average|Under|^0 to 14|^9|^100|65 years and over|85 to 89|15 to 64", popHere$DIM..Age, invert=TRUE), ]
popHere$age = gsub("[[:space:]]+.*", "", popHere$DIM..Age)
colnames(popHere) = gsub(".*(Male|Female)", "\\1", colnames(popHere))
colnames(popHere) = gsub("^GEO_CODE..POR.*", "id", colnames(popHere))
popLong = popHere[, c('level','id', 'GEO_NAME','age','Male','Female')]
popLong2 = reshape2::melt(popLong, id.var = c('level','id','GEO_NAME','age'))
popWide = reshape2::dcast(popLong2, level + id + GEO_NAME~ variable + age, value.var='value')

popSplit = split(popWide[,grep("level", colnames(popWide), invert=TRUE)], popWide$level)
for(D in names(popSplit)) {
	popSplit[[D]][[paste0("id", D)]] = popSplit[[D]]$id
}
for(D in setdiff(names(popSplit), '0')) {
	idHere = paste0('id', D)
	popSplit[[D]] = popSplit[[D]][
		match(bndList[[Dyear]][[idHere]][[idHere]], popSplit[[D]][[idHere]]),]
}
for(D in setdiff(names(popSplit), '0')) {
foreign::write.dbf(
	popSplit[[D]],
	file.path('/store','census','CAN',Dyear, D, 'pop.dbf')
)
	}


	Dyear = 2011
zUrl = 'https://www12.statcan.gc.ca/census-recensement/2011/dp-pd/prof/details/download-telecharger/comprehensive/comp_download.cfm?CTLG=98-316-XWE2011001&FMT=CSV1501&Lang=E&Tab=1&Geo1=PR&Code1=01&Geo2=PR&Code2=01&Data=Count&SearchText=&SearchType=Begins&SearchPR=01&B1=All&Custom=&TABID=1'
zFile = file.path(rawdir, paste('da',Dyear,'pop.zip', sep=''))

if(!file.exists(zFile))
	download.file(zUrl, zFile, method='libcurl')

cFile = unzip(zFile, exdir=rawdir, 
		unzip='/usr/bin/unzip')

cFile = unzip(zFile, list=TRUE, 
		unzip='/usr/bin/unzip')

Sfiles = grep("CSV$", cFile$Name, value=TRUE, ignore.case=TRUE)
Sfiles = grep("Metadata|-DQ[.]", Sfiles, value=TRUE, invert=TRUE, ignore.case=TRUE)

cat('\n')
daPopR = list()
for(D in Sfiles) {
	cat(D, '\n')
	daPop = read.csv(file.path(rawdir,D), stringsAsFactors=FALSE)
	daPop = daPop[grep("Age", daPop$Topic),]
	daPop = daPop[grep("years|Total", daPop$Characteristic),]
	daPop$age = gsub("[[:space:]]|years|population by age groups", "", daPop$Characteristic)
	daPop$age = tolower(gsub("andover", "plus", daPop$age))
	daPop = daPop[grep("^[[:digit:]]+$", daPop$age, invert=TRUE),]
	
	daPop = daPop[,c('Geo_Code','Prov_name','age','Male','Female')]
	names(daPop) = gsub("(em)?ale$", "", names(daPop))
	names(daPop) = gsub("Geo_Code", "id", names(daPop))
	names(daPop) = tolower(gsub("Prov_name", "name1", names(daPop)))
	
	daPopR[[D]] = reshape(
			daPop,
			direction = 'wide',
			idvar = c('id','name1'),
			timevar = 'age',
			sep='_'
			)
}
cat('\n')

daPopAll = do.call(rbind, daPopR)
daPopAll$id0 = 'CAN'
daPopAll$id1 = substr(daPopAll$id,1,2)
daPopAll$id = daPopAll$id5 = as.character(daPopAll$id)


idCols = grep("^id|^name", names(daPopAll), value=TRUE)
daPopAll = daPopAll[,c(idCols, setdiff(names(daPopAll), idCols))]
rownames(daPopAll) = daPopAll$id
daPopAll= daPopAll[match(
	bndList[[Dyear]]$id5$id5,
	daPopAll$id5),]

colnames(daPopAll) = gsub("(to[[:digit:]]|plus).*", "", colnames(daPopAll))
colnames(daPopAll) = gsub("_", ".", colnames(daPopAll))


foreign::write.dbf(
	daPopAll,
	file.path('/store','census','CAN',Dyear, '5', 'pop.dbf')
)
