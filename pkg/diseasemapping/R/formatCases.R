`formatCases` <- function(casedata, ageBreaks=NULL, years=NULL, aggregate.by=NULL) {

# are there age and sex columns?
haveAgeSex = length(grep("^age$", names(casedata), ignore.cases=T)) &
    length(grep("^sex$", names(casedata), ignore.cases=T)) &

# if not, is there an age_sex_group column, in rif format?  use it to create age and sex
groupvar = grep("^AGE_SEX_GROUP$", names(casedata), value=TRUE, ignore.case=TRUE)
 
if(length(groupvar)){
# create age and sex columns
if(length(grep("_", casedata[[groupvar]], ignore.case=TRUE ))){ 
casedata$sex = substr(casedata[[groupvar]], 1, 1)
ss= getBreaks(casedata[[groupvar]])
casedata$age = ss$age
#casedata$group =ss$newNames
}else{
        casedata$sex = factor(substring(casedata[[groupvar]], 1, 1), levels=c("1", "2"), labels=c("M", "F"))
        casedata$age = 5*as.numeric(substring(casedata[[groupvar]], 2, 3))
#        casedata$group = paste(casedata$sex, casedata$age, sep=".")
}     
}

# are there age and sex columns?
if(!length(groupvar)){
haveSex = grep("sex", names(casedata), value=TRUE, ignore.case=TRUE)
haveAge = grep("age", names(casedata), value=TRUE, ignore.case=TRUE)
      if(length(haveSex)){casedata$sex = substr(casedata[[haveSex]], 1, 1)}
      if(length(haveAge)){
             if(length(grep("_", casedata[[haveAge]], ignore.case=TRUE))) {casedata$age = substr(casedata[[haveAge]], 1, 2)}
                   }
#             casedata$groupd = paste(casedata$sex, casedata$age, sep=".")      
                   }
                  
## if there are underscores in age, take the number before the underscores
##if(there are underscores) {
##casedata$age = as.numeric(grep("(_[[:digit:]]+$|PLUS$)", "", casedata$age))
##}

## if not, is there an age_sex_group column, in rif format?  use it to create age and sex

# use the cut function on age

if(!is.null(ageBreaks)) {
 casedata$cutAge = as.character(cut(as.numeric(as.character(casedata$age)), ageBreaks$breaks, right=F))
}

 # aggregate, if necessary
if(!is.null(aggregate.by)) {
   popa = aggregate(casedata$cases, list(age=casedata$cutAge), sum)
   names(popa)[names(popa)=="x"] = "CASES"
   casedata <- popa
}

attributes(casedata)$breaks = ageBreaks
casedata
}
