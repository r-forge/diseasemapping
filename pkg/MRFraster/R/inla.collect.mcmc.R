inla.collect.mcmc <-
function(x,all=T,nBurn=100){

		covs = colnames(x$model.matrix)[-1]
wd.mcmc = x$.args$working.directory
		
dir.mcmc = paste(wd.mcmc, "/", "results.files", sep = "")
dd = dir(dir.mcmc)

#dd = grep("fixed|random|hyperparameter", dd, value=T)

if(all){##if read all results
result = summarys=as.list(dd)
names(result)=names(summarys) = dd



for (id in 1:length(dd)){
d.mcmc = paste(dir.mcmc, "/", dd[id], sep = "")
fnm = paste(d.mcmc, "/", "trace.dat", sep = "")
#print(fnm)

if (file.exists(fnm)) {
allS = read.table(fnm)
result[[id]]  =  data.frame(allS[-c(1:nBurn),])
}
if(class(result[[id]])=="data.frame") {
tempData =  data.frame(result[[id]][-(1:nBurn),])
temp = data.frame(t(rbind(rbind(apply(result[[id]],2,summary), apply(result[[id]],2,quantile,c(0.05,0.95))),sd=apply(result[[id]],2,sd))))
rownames(temp)=NULL
#temp = data.frame(t(rbind(rbind(apply(result[[id]],2,summary), apply(result[[id]],2,quantile,c(0.05,0.95))),sd=apply(result[[id]],2,sd))))
summarys[[id]] =temp[,c(4,9,7,8,1,2,3,5,6)]
}
}
}

##read last configuration only
lastConf = as.list(1:3)
names(lastConf) =c("idTable","Theta","x")
lastDIR = paste(dir.mcmc, "/", dd[grep("^last",dd,ignore.case=T)], sep = "")
fnmTable = paste(lastDIR, "/", "idx-table.dat", sep = "")
fnmTheta = paste(lastDIR, "/", "theta.dat", sep = "")
fnmX= paste(lastDIR, "/", "x.dat", sep = "")
if (file.exists(fnmTable)) {lastConf[["idTable"]] = read.table(fnmTable)}
if (file.exists(fnmTheta)) {lastConf[["Theta"]] = read.table(fnmTheta)}
if (file.exists(fnmX)) {lastConf[["x"]] = read.table(fnmX)}

####read offset
offDir = paste(dir.mcmc, "/", dd[grep("^totaloffset",dd,ignore.case=T)], sep = "")
offTable = paste(offDir, "/", "totaloffset.dat", sep = "")
if (file.exists(offTable)) {totalOffset = read.table(offTable)}


#result
if(all){
x = list(result,summarys,lastConf,totalOffset)
names(x) = c("samples","summary","lastConf","offset")


torename = paste('^fixed\\.effect(0+)',seq(from=2, length=length(covs)), '$',sep="" )
names(torename)=covs

torenamehyp = paste('^hyperparameter([[:digit:]]|-|random\\.effect)+parameter', rep(0:1, c(2,2)), rep(c('','-user-scale'), 2), '$' ,  sep='')
names(torenamehyp) = paste( rep(c('log',''),2),rep(c("Precision", "Range"), c(2,2)),sep='')

torename = c(torename, 'intercept'='^fixed\\.effect(0+)1$', torenamehyp )

for(D1 in c('samples','summary')) {
	for(D2 in 1:length(torename)) 
	names(x[[D1]]) = gsub(torename[D2], names(torename)[D2], names(x[[D1]]))
	
}

}

x
}
