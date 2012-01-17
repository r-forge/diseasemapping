cancerRates = function(area = "canada",
   year=2000, sex=c("M", "F"), site="Lung") {

  sexes = c("M"=1, "F"=2)[sex]



  areaCodes =
  c("Canada"="0099",
  "Newfoundland"=0899,
"Prince Edward Island"=0799,
"Nova Scotia"=0699,
"Ontario"=1199,
"Manitoba"=0399,
"Saskatchewan"=1399,
"Alberta"=0199,
"British Columbia"=0299)
# add more countries!

area = areaCodes[grep(paste("^", area[1], sep=""), names(areaCodes), ignore.case=TRUE)]


result = list()

  rates=NULL
  for(Dsex in names(sexes)) {
fs<-paste("http://ci5.iarc.fr/CI5plus/Table4r.asp?registry=124",area,(paste("&period=",year,sep="",collapse="")),"&sex=", sexes[Dsex],"&window=1&text=1&stat=0&submit=Execute",sep="")
tempn = scan(fs, what="a", quiet=T)
theurl=(paste("http://ci5.iarc.fr/", gsub("^HREF=", "", grep("href=/data", tempn, value=T, ignore.case=T)), sep=""))
result[[Dsex]]=read.table(url(theurl), header=T,skip=1, fill=T, sep="\t", as.is=T)
forAttribute = scan(url(theurl), what="a", sep="\t", nmax=1, quiet=T)

  }


for (Dsex in names(sexes)) {
 cancerTypes =  result[[Dsex]][,1]
  
 cancerTypes = gsub("[[:space:]]+$", "", cancerTypes)

siterow = grep(site, cancerTypes,ignore.case=T)
if(length(siterow) > 1 )    {
  warning(paste("matched ", paste(cancerTypes[siterow], collapse=","), "\n from cancer types", paste(cancerTypes, collapse=","), 
  "\n using", cancerTypes[siterow[1]]))
  siterow=siterow[1]
  }
if(length(siterow) ) {
  x=result[[Dsex]][siterow,]
iarcSite =gsub("[[:space:]]+$", "", x[1]) 
  x = x[,-c(1,grep("Total|Unk|ICD|CR|ASR", names(x), ignore.case=T))]
  x = as.vector(x)
  names(x) = gsub("^X|\\.$", "", names(x))
  rates[[Dsex]]=x
  }
}
rates = unlist(rates)/100000
names(rates) = gsub("\\.", "_", names(rates))

attributes(rates)$site = iarcSite
forAttribute = strsplit(forAttribute, "\\(")[[1]]
attributes(rates)$area = gsub("[[:space:]]+$", "", forAttribute[1])
attributes(rates)$year = gsub("\\)$", "", forAttribute[2])

return(rates)

}


