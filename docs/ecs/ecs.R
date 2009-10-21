
library(ECS)
data(controls)
data(cases)

names(controls) = gsub("^smoked$", "smk", names(controls))
controls$cancer = 0
cases$cancer = as.integer(gsub("[[:space:]]", "", cases$cancer))
ECS = rbind(cases[,names(controls)], controls)

ECS$cancer = factor(ECS$cancer, levels = c(0, 2, 5, 6), labels = c("control","stomach","lung","breast"))

names(ECS) = gsub("^smk$", "smk_ever", names(ECS))

# get rid of spaces
for(D in c("sex","smk_ever", "smk_now")) {
  ECS[[D]] = gsub("[[:space:]]", "", ECS[[D]])
}
# change never smokers to not smoking now
ECS$smk_now[ECS$smk_ever == "N"] = "N"

# change Y,N to T,F
for(D in c("smk_ever", "smk_now"))
  ECS[[D]] = ECS[[D]]=="Y"

ECS$smoke = "never"
ECS$smoke[ECS$smk_ever] = "former"
ECS$smoke[ECS$smk_now] = "current"
ECS$smoke = factor(ECS$smoke, levels=c("never","former","current")   )

# integers
for(D in c("diagage","smk_age","smk_avg","smk_quit","packyrs"))
  ECS[[D]] = as.integer(ECS[[D]])

save(ECS, file="ECS.RData")

