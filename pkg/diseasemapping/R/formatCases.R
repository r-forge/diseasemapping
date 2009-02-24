`formatCases` <- function(casedata, ageBreaks=NULL, years=NULL, aggregate.by=NULL) {

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

# this can probably go (?)
 
 groupvar = grep("^AGE_SEX_GROUP$", names(casedata), value=TRUE, ignore.case=TRUE)
 if(length(groupvar) & ageBreaks$mustAggregate ==TRUE)  {
 # there is a group variable, use it to create age and sex
  # get rid of NA's
  casedata = casedata[!is.na(casedata[[groupvar]]),]
  # change NA's in cases to zero
 casedata$cases[which(is.na(casedata$cases))]= 0
   casedata$cutAge = cut(as.numeric(as.character(casedata$age)), ageBreaks$breaks, right=F)
   popa = aggregate(casedata$cases, list(age=casedata$cutAge), sum)
   names(popa)[names(popa)=="x"] = "CASES"
   casedata <- popa
}
 # check for underscores
if(length(groupvar)& ageBreaks$mustAggregate == FALSE)  {
 # there is a group variable, use it to create age and sex
  # get rid of NA's
  casedata = casedata[!is.na(casedata[[groupvar]]),]
 underscores = as.logical(length(grep("_", casedata[[groupvar]])))
 if(underscores) {
    n <- regexpr("_", casedata[[groupvar]], fixed = TRUE)
    casedata$age = cut(as.integer(substr(casedata[[groupvar]], n-2, 100)), breaks=ageBreaks$breaks)
    casedata$sex = factor(substr(casedata[[groupvar]], 1, 1))
    casedata$AGE_NUM= as.numeric(substring(casedata$age, 0, 2))

 } else {
 # assume we have RIF format age groups
 casedata$sex = factor(substring(casedata[[groupvar]], 1, 1),
    levels=c("1", "2"), labels=c("M", "F"))
 ageterm = 5*as.numeric(substring(casedata[[groupvar]], 2, 3))

 casedata$age = factor(paste(ageterm, ageterm+4, sep="_"))
 }
 }

casedata
}

