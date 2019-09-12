## ----packages, results='hide'--------------------------------------------
library('gpuR')

## ----memoryAvailable, echo=TRUE------------------------------------------
gpuInfo()$deviceName
gpuInfo()$maxAllocatableMem/(1024^3)

