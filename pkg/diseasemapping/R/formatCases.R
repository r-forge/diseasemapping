`formatCases` <- function(casedata) {


 groupvar = grep("^AGE_SEX_GROUP$", names(casedata), value=TRUE,
    ignore.case=TRUE)
 if(length(groupvar))  {
 # there is a group variable, use it to create age and sex
  # get rid of NA's
  casedata = casedata[!is.na(casedata[[groupvar]]),]

 
 # check for underscores
 underscores = as.logical(length(grep("_", casedata[[groupvar]])))
 if(underscores) {
    n <- regexpr("_", casedata[[groupvar]], fixed = TRUE)
    casedata$age = factor(substr(casedata[[groupvar]], n-2, 100))
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

