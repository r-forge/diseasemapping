
Dyear = 2011
Dlevel = 5
basedir = '/store/census'
rawdir = file.path(basedir, 'raw')

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

daPopAll[1:4,-(10:30)]

foreign::write.dbf(daPopAll, 
		file.path(basedir, 'CAN', Dyear, Dlevel,"population.dbf")
)


database = c(dbname='gcensus', user='gcensus', host='localhost',  port=5433)

con <- do.call(
		DBI::dbConnect,
		c(
				list(drv=RPostgreSQL::PostgreSQL()),
				as.list(database)
		)
)

tableName = paste("can", Dyear, ".l",Dlevel,"pop", sep="")

if(DBI::dbExistsTable(con,tableName)) {
	DBI::dbRemoveTable(con,tableName)
}


DBI::dbWriteTable(con,
		tableName, 
		daPopAll, row.names = FALSE)

DBI::dbDisconnect(con)
