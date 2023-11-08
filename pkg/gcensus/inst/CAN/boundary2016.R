

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


reversePoly = function(x) {
  theDf = values(x)
  toSwap = theDf$TNODE_
  theDf$TNODE_ = theDf$FNODE_
  theDf$FNODE_ = toSwap
  vect(crds(x)[nrow(crds(x)):1, ], atts=theDf, crs=crs(x), type='lines')
}


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
	bnd1a = terra::vect(De00file, layer='ARC')
	
	
#	(toronto = values(bndLabels)[bndLabels$CDUID == 3520, ])
	
#	plot(bnd1a, xlim = -79.4 + c(-1,1), ylim = 43.6 + c(-1,1))
	
	bndLabels = terra::vect(De00file, layer='LAB')
	
	IDvar = grep(paste0(Sid[Dlevel], 'UID$'), names(bndLabels), value=TRUE, ignore.case=TRUE)
	
	bndLabels1 = split(values(bndLabels), bndLabels[[IDvar]])
	
	verbose = FALSE
	bndPoly = list()
	for(D in 1:length(bndLabels1)) {
	  if(verbose) cat(D, " ")
	  xx = bndLabels1[[D]] #lapply(bndLabels1, function(xx) {
	  res = bnd1a[bnd1a$RPOLY_ %in% xx$PolyId | bnd1a$LPOLY_ %in% xx$PolyId,]
	  res$isRight = res$RPOLY_ %in% xx$PolyId
	  res$isLeft = res$LPOLY_ %in% xx$PolyId
	  isBoth = which(res$isRight & res$isLeft)
	  if(length(isBoth)) {
	    # get rid, it's dividing a region in two
	    res = res[-isBoth, ]
	  }
	  res$isSame = res$FNODE_ == res$TNODE_
	  res=res[order(!res$isSame, res$FNODE_),]
	  res$polyId = NA
	  resOutPoints = res[res$isSame, ]
	  resOutPoints$polyId = 1:nrow(resOutPoints)
	  resRemaining = res[!res$isSame, ]
	  tableOfNodes= sort(table(c(resRemaining$FNODE_, resRemaining$TNODE_)))
	  theFours = as.numeric(names(tableOfNodes[tableOfNodes==4]))
	  
	  
	  #	  theOnes = as.integer(names(which(table(unlist(values(resRemaining)[!isTheFours,c("FNODE_", 'TNODE_')]))==1)))
	  
	  #isTheOnes = resRemaining$FNODE_ %in% theOnes | resRemaining$TNODE_ %in% theOnes
	  #isTheOnes = isTheOnes & !isTheFours
	  
#	  resIsOnes = resRemaining[isTheOnes, ]
	  resRemainingOrig = resRemaining
	  if(length(theFours)) {
	    isTheFours = mapply(function(fourHere) {
	      which(resRemaining$FNODE_ %in% fourHere | resRemaining$TNODE_ %in% fourHere)
	    }, fourHere = theFours, SIMPLIFY=TRUE)

	  theFoursReserve = isTheFours#[-1, , drop=FALSE]

	  resRemaining = resRemainingOrig[-as.vector(theFoursReserve), ]
	  resTheFours = resRemainingOrig[as.vector(theFoursReserve), ]
	  theFoursReserveId = apply(values(resTheFours)[,c('FNODE_', 'TNODE_')], 1, setdiff, y=theFours)
	  theFoursIsT = resTheFours$TNODE_ %in% theFoursReserveId
	} # length(theFours)
	  
	  
	  theTable = table(unlist(values(resRemaining)[,c("FNODE_", 'TNODE_')]))
	  theOnes = as.numeric(names(which(theTable==1)))

	  resRemaining$oneNeighbour = resRemaining$FNODE_ %in% theOnes | resRemaining$TNODE_ %in% theOnes
	  
	  resRemaining$oneNeighbour = resRemaining$oneNeighbour - 5*(resRemaining$FNODE_ %in% theFours) - 5*(resRemaining$TNODE_ %in% theFours)
	  
	  resRemaining=resRemaining[order(!resRemaining$oneNeighbour, resRemaining$FNODE_),]
	  
	  Dpolygon = nrow(resOutPoints)+1
	  resOut = resRemaining[1,]
	  resOut$polyId = Dpolygon
	  
	  # reverse if can't find a 'to' node but do have a 'from'
	  if(nrow(resRemaining)) {
	  if(!( resOut$TNODE_ %in% resRemaining$FNODE_) ){
  	  resOut = reversePoly(resOut)
	  }
	  }
	  
	  resRemaining = resRemaining[-1, ]
	  Sline  = seq(from=1, len=nrow(resRemaining), by=1)
	  for(Dline in Sline) {
	    inTnode = which(resRemaining$TNODE_ == resOut$TNODE_[nrow(resOut)])
	    TnodeHere = resRemaining$TNODE_[inTnode]
	    inFnode = which(resRemaining$FNODE_ == resOut$TNODE_[nrow(resOut)])
	    FnodeHere = resRemaining$FNODE_[inFnode]
	    inC = c(inFnode, inTnode)
	    if(length(inC) > 1) {
	      warning("found two polygons", Dline)
	      # check if can close a loop
	      firstInThisPolyid = which.min(resOut$polyId == Dpolygon)
	      toMatch = values(resOut)[firstInThisPolyid, 'FNODE_']
	      toMatch1 = values(resRemaining)[inFnode,'TNODE_'] == toMatch
	      toMatch2 = values(resRemaining)[inTnode,'FNODE_'] == toMatch
        if(any(toMatch1)) {
          inC = inFnode[toMatch1][1]
        } else if(any(toMatch2)) {
          inC = inTnode[toMatch2][1]
        }
	      sameLeft = resRemaining$isLeft[inC] == resOut$isLeft[[nrow(resOut)]]
	      if(any(sameLeft)) inC = inC[sameLeft]
	      inC = inC[1]
	    }
	    if(!length(inC) ) {
	      # new polygon
	      Dpolygon = Dpolygon + 1
	      inFnode = inC = 1
	      newRes = resRemaining[1, ]
	      nodesRemaing = unlist(values(resRemaining)[-1,c('FNODE_', 'TNODE_')])
	      matchedTnode =  newRes$TNODE_ %in% nodesRemaing
	      matchedFnode = newRes$TNODE_ %in% nodesRemaing
	                                          
	      if(matchedFnode & (!matchedTnode)){
	        newRes = reversePoly(newRes)
	      }
	    } else {
  	    newRes = resRemaining[inC, ]
	    }
	    
	    newRes$polyId = Dpolygon
	    resRemaining = resRemaining[-inC, ]
	    if(length(inTnode)) { # reverse
	      newRes = reversePoly(newRes)
	    }
	    resOut = vect(c(resOut, newRes))
	  } # Dline
	  
	  if(length(theFours)) {
	    # put nodes with four connections back in
      for(Dfour in 1:nrow(resTheFours)) {
        newRes = resTheFours[Dfour,]
        matchF = which(resOut$FNODE_ == theFoursReserveId[Dfour])
        matchT = which(resOut$TNODE_ == theFoursReserveId[Dfour])

        newRes$polyId = resOut$polyId[c(matchF, matchT)[1]]
          
        if( !(
          (length(matchF) & (theFoursIsT[Dfour]) ) 
            | 
          (length(matchT) & (!theFoursIsT[Dfour]) ) 
        )
          ){
          # reverse
          newRes = reversePoly(newRes)
        }
        if(length(matchT)) { # goes after
          resOut = vect(c(
            resOut[seq(from=1, to=matchT[1])],
            newRes,
            resOut[seq(from=matchT[1]+1, by=1, len=nrow(resOut)-matchT[1])]
          ))
        } else if(length(matchF)) {  # before
          resOut = vect(c(
            resOut[seq(from=1, by=1, len=matchF[1]-1)],
            newRes,
            resOut[seq(matchF[1], nrow(resOut))]
          ))
        }
      }
	    
	  } # adding back in theFours
	  
	  resOut2 = vect(c(resOutPoints, resOut))
	  resOut3 = as.points(resOut2)
	  resOut4 = split(resOut3, resOut3$polyId)

	  resOut5 = lapply(resOut4, function(xxx) {
	    vect(crds(xxx), atts=values(xxx), crs=crs(xxx), type='polygons')
	  })
	  resOut6 = aggregate(vect(resOut5))	  
	  values(resOut6) = xx[1,grep("NAME|ID", names(xx))]
	  bndPoly[[D]] = resOut6
	}	# D
	if(verbose) cat("\n ")
	
	
	bnd = vect(bndPoly)
  

	if(FALSE) {
		mapmisc::map.new(bndPoly1[grep("Toronto|Essex|Ottawa", bndPoly1$CDNAME), ])
		plot(bndPoly1, add=TRUE)
	}


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
# Dyear = 2016

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
