#' # 2011
#' 
#' ## Da
#' 
#' ### population
#' 
zurl = 'http://dc2.chass.utoronto.ca/grid2/census/uixCoSIQqSq_data.csv.zip'
hurl = 'http://dc2.chass.utoronto.ca/grid2/census/uixCoSIQqSq_header.txt'
file = basename(zurl)
hfile = basename(hurl)
if(!file.exists(file))
	download.file(zurl, file)
if(!file.exists(hfile))
	download.file(hurl, hfile)
unzip(file, list=TRUE)
unzip(file)

dat = read.csv(unzip(file, list=TRUE)$Name, sep=',', header=TRUE)
header = read.table(hfile, skip=1, sep='-', fill=TRUE, as.is=TRUE, col.names=c('a','b','c'))
header[nchar(header$c)==0 ,'c'] = header[nchar(header$c) ==0 ,'b'] 

cnames = gsub("Total.*groups", "", header[,'c'])
names(cnames) = trimws(header[,'a']) 
cnames = gsub("Females / ", "f_", cnames)
cnames = gsub("Males / ", "m_", cnames)
cnames = gsub(" to ", "to", cnames)
cnames = gsub("years", "", cnames)
cnames = gsub("[[:space:]]+and over", "plus", cnames)
cnames = gsub("^[[:space:]]?(m|f)_;", "", cnames)
cnames = gsub("; (Fem|M)ales$", "", cnames)
cnames[grep("Both sexes", cnames)] = 'Total'
toJunk = grep(";|Population", cnames)
cnames[toJunk] = paste("junk", seq(1, length(toJunk)))
cnames = trimws(cnames)
cnames = gsub("%", "pct", cnames)
cnames = gsub("CD code", "id2", cnames)
cnames = gsub("CD name", "name2", cnames)
cnames = gsub("Province code", "id1", cnames)
cnames = gsub("Province name", "name1", cnames)
cnames = gsub("GEO UID", "id5", cnames)
cnames = gsub("[[:space:]]", "", cnames)

names(dat) = cnames[names(dat)]

dat = dat[, grep('junk', names(dat), invert=TRUE)]
dat = dat[dat$DAname > 0,]

foreign::write.dbf(dat, "/store/census/CAN/2011/5/population.dbf")


#' ## CT
#' 
#' ### population
#' 
zurl = 'http://dc2.chass.utoronto.ca/grid2/census/NLR63uPrSKK_data.csv.zip'
hurl = 'http://dc2.chass.utoronto.ca/grid2/census/NLR63uPrSKK_header.txt '
file = basename(zurl)
hfile = basename(hurl)
if(!file.exists(file))
	download.file(zurl, file)
if(!file.exists(hfile))
	download.file(hurl, hfile)
unzip(file, list=TRUE)
unzip(file)

dat = read.csv(unzip(file, list=TRUE)$Name, sep=',', header=TRUE)
header = read.table(hfile, skip=1, sep='-', fill=TRUE, as.is=TRUE, col.names=c('a','b','c'))
header[nchar(header$c)==0 ,'c'] = header[nchar(header$c) ==0 ,'b'] 

cnames = gsub("Total.*groups", "", header[,'c'])
names(cnames) = trimws(header[,'a']) 
cnames = gsub("Females / ", "f_", cnames)
cnames = gsub("Males / ", "m_", cnames)
cnames = gsub(" to ", "to", cnames)
cnames = gsub("years", "", cnames)
cnames = gsub("[[:space:]]+and over", "plus", cnames)
cnames = gsub("^[[:space:]]?(m|f)_;", "", cnames)
cnames = gsub("; (Fem|M)ales$", "", cnames)
cnames[grep("Both sexes", cnames)] = 'Total'
toJunk = grep(";|Population", cnames)
cnames[toJunk] = paste("junk", seq(1, length(toJunk)))
cnames = trimws(cnames)
cnames = gsub("%", "pct", cnames)
cnames = gsub("CD code", "id2", cnames)
cnames = gsub("CD name", "name2", cnames)
cnames = gsub("Province code", "id1", cnames)
cnames = gsub("Province name", "name1", cnames)
cnames = gsub("GEO UID", "id", cnames)
cnames = gsub("[[:space:]]", "", cnames)

names(dat) = cnames[names(dat)]

dat = dat[, grep('junk', names(dat), invert=TRUE)]

foreign::write.dbf(dat, "/tmp/population.dbf")
foreign::write.dbf(dat, "/store/census/CAN/2011/4/population.dbf")



#' ## CSD
#' 
#' ### population
#' 
zurl = 'http://dc2.chass.utoronto.ca/grid2/census/ziNvk2EIJECAB6kcp_data.csv.zip'
hurl = 'http://dc2.chass.utoronto.ca/grid2/census/ziNvk2EIJECAB6kcp_header.txt '
file = basename(zurl)
hfile = basename(hurl)
if(!file.exists(file))
	download.file(zurl, file)
if(!file.exists(hfile))
	download.file(hurl, hfile)
unzip(file, list=TRUE)
unzip(file)

dat = read.csv(unzip(file, list=TRUE)$Name, sep=',', header=TRUE)
header = read.table(hfile, skip=1, sep='-', fill=TRUE, as.is=TRUE, col.names=c('a','b','c'))
header[nchar(header$c)==0 ,'c'] = header[nchar(header$c) ==0 ,'b'] 

cnames = gsub("Total.*groups", "", header[,'c'])
names(cnames) = trimws(header[,'a']) 
cnames = gsub("Females / ", "f_", cnames)
cnames = gsub("Males / ", "m_", cnames)
cnames = gsub(" to ", "to", cnames)
cnames = gsub("years", "", cnames)
cnames = gsub("[[:space:]]+and over", "plus", cnames)
cnames = gsub("^[[:space:]]?(m|f)_;", "", cnames)
cnames = gsub("; (Fem|M)ales$", "", cnames)
cnames[grep("Both sexes", cnames)] = 'Total'
toJunk = grep(";|Population", cnames)
cnames[toJunk] = paste("junk", seq(1, length(toJunk)))
cnames = trimws(cnames)
cnames = gsub("%", "pct", cnames)
cnames = gsub("CD code", "id2", cnames)
cnames = gsub("CD name", "name2", cnames)
cnames = gsub("CSD name", "name", cnames)
cnames = gsub("Province code", "id1", cnames)
cnames = gsub("Province name", "name1", cnames)
cnames = gsub("GEO UID", "id", cnames)
cnames = gsub("[[:space:]]", "", cnames)

names(dat) = cnames[names(dat)]

dat = dat[, grep('junk', names(dat), invert=TRUE)]

foreign::write.dbf(dat, "/tmp/population.dbf")
foreign::write.dbf(dat, "/store/census/CAN/2011/3/population.dbf")




#' # 1996
#' 
#' ## Da
#' 
#' ### population
#' 
zurl = 'http://dc2.chass.utoronto.ca/grid2/census/9q6nt3OdfNQPRNd_data.csv'
hurl = 'http://dc2.chass.utoronto.ca/grid2/census/9q6nt3OdfNQPRNd_header.txt'
file = basename(zurl)
hfile = basename(hurl)
if(!file.exists(file))
	download.file(zurl, file)
if(!file.exists(hfile))
	download.file(hurl, hfile)

dat = read.csv(file, sep=',', header=TRUE)
header = scan(hfile, skip=1, what='a', sep='\n')
names(header) = gsub(" - .*$", "", header)
cnames = gsub("^COL[[:digit:]]+ - ", "", header)

cnames = gsub("Female, total / ", "f_", cnames)
cnames = gsub("Male, total / ", "m_", cnames)
cnames = gsub("-", "to", cnames)
cnames = gsub("\\+", "plus", cnames)
cnames = gsub("^[[:space:]]?(m|f)_;", "", cnames)
cnames = trimws(cnames)
cnames = gsub("Province code", "id1", cnames)
cnames = gsub("Province name", "name1", cnames)
cnames[grep("^EAuid", cnames)] = "id5"
cnames = gsub("[[:space:]]", "", cnames)
cnames[grep("Male,total", cnames)] = 'm_total'
cnames[grep("Female,total", cnames)] = 'm_total'
cnames[grep("^Totalpopulation", cnames)] = 'a_total'


names(dat) = cnames[names(dat)]
dat = dat[dat$EnumerationArea > 0,]
dat = dat[,grep("^(m|f|id|a|name)", names(dat), value=TRUE)]
dat = cbind(id = dat$id5, dat)
foreign::write.dbf(dat, "/store/census/CAN/1996/5/population.dbf")

#' ## CT

zurl = 'http://dc2.chass.utoronto.ca/grid2/census/1455674033_LVAE76OrkGTbqILJ.dbf'

hurl = 'http://dc2.chass.utoronto.ca/grid2/census/1455674033_LVAE76OrkGTbqILJ.txt'
file = basename(zurl)
hfile = basename(hurl)
if(!file.exists(file))
	download.file(zurl, file)
if(!file.exists(hfile))
	download.file(hurl, hfile)
dat = foreign::read.dbf(file)
header = read.table(hfile, skip=0, sep=',', stringsAsFactors=FALSE)
cnames = header[,2]
names(cnames) = header[,1]

mPos = grep("^Male", cnames)
fPos = grep("^Female", cnames)

numericPos = grep("^[[:digit:]]", cnames)
cnames[numericPos[numericPos < fPos]] = paste("m_", 
		cnames[numericPos[numericPos < fPos]], sep=""
		)
		cnames[numericPos[numericPos > fPos]] = paste("f_", 
				cnames[numericPos[numericPos > fPos]], sep=""
		)
cnames = gsub("^CTUID.*", "id4", cnames)
cnames = gsub("Male total", "m_total", cnames)
cnames = gsub("Feale total", "f_total", cnames)
cnames = gsub("Female total", "f_total", cnames)
cnames = gsub("^Population 1996.*", "a_total", cnames)
cnames = gsub("\\$", "", cnames)
cnames = gsub("[[:space:]]", "", cnames)
cnames = gsub("-", "to", cnames)

names(dat) = cnames[names(dat)]
dat = dat[dat$CTCode != '0000',]
dat = dat[,grep("^(m|f|id|a|name)", names(dat), value=TRUE)]

#prov = readOGR("/store/census/CAN/1996/1", "boundary")
bnd = readOGR("/store/census/CAN/1996/4", "boundary")

dat = cbind(id = dat$id4, 
		id1 = bnd$id1[match(dat$id4, bnd$id)],
		dat)
foreign::write.dbf(dat, "/store/census/CAN/1996/4/population.dbf")
