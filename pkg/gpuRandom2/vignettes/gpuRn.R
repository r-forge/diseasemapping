## ----date----------------------------------------------------------------
date()

## ----packages, results='hide'--------------------------------------------
library('gpuR')
# set context to the second gpu
setContext( grep('gpu', listContexts()$device_type) [2]    )


## ----setType-------------------------------------------------------------
theType = c('float','double')[1+gpuInfo()$double_support]
theType

