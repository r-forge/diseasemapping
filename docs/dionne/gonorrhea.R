library(maptools)
library(spdep)
popdata <-  readShapePoly(fn="NC_gonorrhea_yearly_and_quarterly_with_covariates")


popAdjacency = poly2nb(popdata, as.character(popdata[["CENTRACTID"]]))
load("adds2.Rdata")


adds2$CenTractID = as.character(adds2$CenTractID)
popdata$CENTRACTID = as.character(popdata$CENTRACTID)

newdata <- merge(popdata, adds2, by.x="CENTRACTID", by.y="CenTractID", all.x=TRUE)


newdata$RUCA = factor(newdata$RUCA)
newdata$pctNotOwned = 1-newdata$PCT_OWNER

newdata$MFlog = log(newdata$MF_RATIO)
newdata$MFlog[newdata$MFlog > 1] = 1


# reshape the newdata to long format

ColPop = grep("^(Pop[[:digit:]]+$)", names(newdata), value = T, ignore.case = T)
# get rid of 2004 and 2005q1
yearsToGetRidOf = grep("^CasesY04Q[[:digit:]]$|^CasesY05Q1$", names(newdata), value = T, ignore.case = T)


ColCase = grep("^CasesY0[[:digit:]]Q[[:digit:]]$", names(newdata), value = T, ignore.case = T)
years = substr(ColCase, 7, 8)
# get rid of 2004 and 2005q1
yearsToGetRidOf = grep("^CasesY04Q[[:digit:]]$|^CasesY05Q1$", names(newdata), value = T, ignore.case = T)
ColCase = ColCase[! ColCase %in% yearsToGetRidOf]

years = gsub("^CasesY|Q[[:digit:]]$", "", ColCase)

ColPop = paste("Pop20", years, sep="")

poplong = reshape(newdata,  varying=list(ColPop, ColCase), direction="long",
	 timevar="time", times = substr(ColCase, 7, 10), v.names=c("Population", "Cases"))

poplong$logPop = log(poplong$Population) - log(4)

# centre some covariates
for(Dcentre in c("PCT__30K", "FHH_CHILD", "PCT_RENTER", "PCT__HS", "PCT_UNEMP", "PCTBLACK"))
  poplong[,Dcentre] = poplong[,Dcentre] - mean(poplong[,Dcentre])
     

library(glmmBUGS)


GONQragged <- glmmBUGS(Cases + logPop ~ time +  MFlog + 
    PCTBLACK+ FHH_CHILD+ PCT_NO_PLM +
    PCT__HS + RUCA * (PCT_RENTER +PCT__30K) + PCT_UNEMP, 
                        data=poplong,
                      effects = c("CENTRACTID", "time"), spatial=popAdjacency, spatialEffect= "CENTRACTID", family="poisson")

save(GONQragged, file="GONQragged.RData")

load("GONQragged.RData")
startingValues = GONQragged$startingValues
source("getInits.R") 
library(R2WinBUGS)


GONQResult = bugs(GONQragged$ragged, getInits, 
	parameters.to.save = names(getInits()), 
	model.file="model.bug", n.chain=3, n.iter=30000, n.burnin=2000, n.thin=100, program="winbugs",useWINE = TRUE, debug=T,
bugs.directory="C:/Program Files/WinBUGS64/WinBUGS14/")   
save(GONQResult, GONQragged, file="GONQResult.RData")


 gonpar = restoreParams(GONQResult, GONQragged$ragged)
gonsummary = summaryChain(gonpar)
library(diseasemapping)
pop2 = mergeBugsData(popdata, gonsummary, "CENTRACTID") 


# code to reformat the time bit

# split the Rtime bit by region and time
thenames = dimnames(gonsummary$Rtime)[[1]]
thenames = matrix(unlist(strsplit(thenames, "/")), ncol=2, byrow=T, dimnames=list(NULL, c("region", "time")))

Stime = unique(thenames[,"time"])
theregionnames = unique(thenames[,"region"])
summaryNames = dimnames(gonsummary$Rtime)[[2]]
Rtime = array(NA, 
	c(length(theregionnames), length(Stime), length(summaryNames)),
	dimnames = list(theregionnames, Stime, summaryNames) )
for(Dtime in Stime) {

  tomerge = thenames[,"time"]==Dtime

	Rtime[thenames[tomerge, "region"], Dtime, ] = 
	 gonsummary$Rtime[tomerge, ]
}

RtimeMean = Rtime[,,"mean"]
colnames(RtimeMean) = paste("mean", colnames(RtimeMean), sep="")

ncFinal = pop2
ncFinal@data = merge(pop2@data, RtimeMean, by.x="CENTRACTID", 
  by.y="row.names", all.x=T) 
save(ncFinal, file="ncResultsSpaceTime.RData")
writePolyShape(ncFinal, "ncResultsSpaceTime")
