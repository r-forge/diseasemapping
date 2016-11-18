
Dyear = 2016
basedir = '/store/census'
rawdir = file.path(basedir, 'raw')
candir = file.path(basedir, 'CAN', Dyear)

dir.create(candir)

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

Surl = paste(
		'http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/l',
		SidEx, '000b16a_e.zip', sep=''
)
Zfile = file.path(rawdir, basename(Surl))

Sfile = Pmisc::downloadIfOld(
		Surl, Zfile, verbose=TRUE, exdir=rawdir
)


SshpFile = grep("[.]shp$", Sfile, value=TRUE)
SshpFile = gsub("[.]shp$", "", basename(SshpFile))

library("rgdal")
canada = raster::getData("GADM", country='CAN', level=0, path=rawdir)
canada = spTransform(canada, mapmisc::crsLL)
canadaBg = mapmisc::openmap(canada)

for(Dlevel in names(Sid)) {
	cat(Dlevel, ' ')
	DlevelId = as.numeric(gsub("id", "", Dlevel))
	leveldir = file.path(candir, DlevelId)
	dir.create(leveldir)
	
	DshpFile = grep(
			paste("l", Sid[Dlevel], sep=''),
			SshpFile, value=TRUE
	)
	
	bnd = readOGR(rawdir, DshpFile, stringsAsFactors=FALSE)
	bnd = spTransform(bnd, mapmisc::crsLL)
	
	bndDf = data.frame(
			id0 = rep('CAN', length(bnd)),
			name0 = 'Canada'
	)
	
	for(DlevelNames in names(SidNames)) {
		
		idHere = SidNames[[as.character(DlevelNames)]]
		
		for(Dname in names(idHere)) {
			if(idHere[Dname] %in% names(bnd))
				bndDf[[Dname]] = bnd@data[[idHere[Dname]]]
		}
		
	}
	bnd@data = bndDf

	# write map
	
	png(file.path(leveldir, "map.png"))
	mapmisc::map.new(canada, buffer=-c(0,0.5,0.5,5))
	plot(canadaBg, add=TRUE)
	plot(bnd, add=TRUE, border='red')
	dev.off()
	
	# write shapefile
	
	writeOGR(
			bnd,
			leveldir,
			'map',
			driver= "ESRI Shapefile", overwrite=TRUE
	)
	
	
	
}


