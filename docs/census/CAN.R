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